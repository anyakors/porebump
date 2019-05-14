import numpy as np
from bokeh.client import push_session
from bokeh.io import curdoc
from bokeh.plotting import figure
from bokeh.layouts import column, widgetbox
from bokeh.palettes import Spectral11
from bokeh.models.widgets import TextInput
import matplotlib.pyplot as plt
import h5py
import threshold_utils as thu
import os


files = os.listdir('/var/lib/MinKNOW/data/reads/20181213_0929_RNA11/fast5/pass/0')
i = 0

for file in files:
    if file.endswith('.fast5'):
        f = h5py.File('/var/lib/MinKNOW/data/reads/20181213_0929_RNA11/fast5/pass/0/{}'.format(file), 'r')
        read_no = list(f['Raw']['Reads'].keys())
        signal = f['Raw']['Reads'][read_no[0]]['Signal'][()]


#y = np.loadtxt('/Users/iamqoqao/workspace/minion/training_sets/truncated/rna3/ch_4_7', encoding='latin1', delimiter="\n")
f = h5py.File('/Users/iamqoqao/workspace/minion/data/20181213_0929_RNA11/fast5/pass/0/mookse_Veriton_X4650G_20181213_FAH54070_MN21778_sequencing_run_RNA11_27121_read_2037_ch_242_strand.fast5', 'r')
y = [f['Raw']['Reads'][key]['Signal'][()] for key in f['Raw']['Reads'].keys()]
y = np.array(y[0])
x = np.arange(0, len(y))
filename = 'ch_2_20'

lag = 30
threshold = 5
influence = 0.0

result = thu.thresholding_algo(y, lag=lag, threshold=threshold, influence=influence)
avg = result["avgFilter"]
std = result["stdFilter"]
signals = result["signals"]
x2 = np.arange(0, len(signals))

frequency = thu.peak_frequency(result["signals"])
#freq_threshold = frequency*20
freq_threshold = 2000

drops = thu.drop_occurence(signals, freq_threshold)
x3 = np.arange(0, len(drops))

# figure() function auto-adds the figure to curdoc()
p1 = figure(title='rna9 {}'.format(filename), x_range=(min(x), max(x)), y_range=(min(y), max(y)), plot_width=1000, plot_height=200)
r1 = p1.line(x=[], y=[], line_width=1, color='navy')
#mypalette=Spectral11[0:3]
#r1 = p1.multi_line(xs=[], ys=[], line_color=mypalette)

p2 = figure(x_range=(min(x), max(x)), y_range=(-2, 2), plot_width=1000, plot_height=200)
r2 = p2.step(x=[], y=[], line_width=1, color='green')

p3 = figure(x_range=(min(x), max(x)), y_range=(-2, 2), plot_width=1000, plot_height=200)
r3 = p3.step(x=[], y=[], line_width=1, color='red')
# open a session to keep our local document in sync with server
session = push_session(curdoc())
i = 0
ds1 = r1.data_source
ds2 = r2.data_source
ds3 = r3.data_source

def update():
    global i
    ds1.data['x'].append(x[i])
    ds1.data['y'].append(y[i])
    #ds1.data['xs'].append([x[i], x[i], x[i]])
    #ds1.data['ys'].append([y[i], avg[i]+threshold*std[i], avg[i]-threshold*std[i]])
    ds1.trigger('data', ds1.data, ds1.data)

    ds2.data['x'].append(x2[i])
    ds2.data['y'].append(signals[i])
    ds2.trigger('data', ds2.data, ds2.data)

    ds3.data['x'].append(x3[i])
    ds3.data['y'].append(drops[i])
    ds3.trigger('data', ds3.data, ds3.data)

    if i < len(x)-1:
        i += 1
    else:
        i = 0
        ds1.data['x'] = []
        ds1.data['y'] = []

        ds2.data['x'] = []
        ds2.data['y'] = []

        ds3.data['x'] = []
        ds3.data['y'] = []

text_input = TextInput(value="default", title="Label:")

doc = curdoc()
doc.add_periodic_callback(update, 1)

session.show(column(widgetbox(text_input), p1, p2, p3)) # open the document in a browser
session.loop_until_closed() # run forever