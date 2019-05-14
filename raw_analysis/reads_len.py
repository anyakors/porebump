#extracts random strand events from a given bulk file 

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import UnivariateSpline
import argparse
import statistics
import seaborn as sns

sns.set()

def find_read(name):
    if 'Read_' in name:
        return name

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
    help="path to folder with fast5 files")
#ap.add_argument("-c", "--channel", required=True,
#    help="channel number")
ap.add_argument("-o", "--output", required=True,
    help="path to savedir")
ap.add_argument("-p", "--plot", required=True,
    help="path to plot file savedir")
args = vars(ap.parse_args())

fast5s = os.listdir(args["input"])

strand_events = []

for file in fast5s:

    if file.endswith("fast5"):

        print(file)
        fast5 = h5py.File(os.path.join(args["input"], file), 'r')

        for key in fast5['Raw']['Reads'].keys():
            signal = fast5['Raw']['Reads'][key]['Signal'][()]
            signal = signal.astype(np.int16) 
            key_str = fast5['Raw']['Reads'].visit(find_read)

        strand_events.append(len(signal))


print("After", len(strand_events), "strands, mean =", statistics.mean(strand_events), "std =", statistics.stdev(strand_events))

np.savetxt(args["output"], strand_events, delimiter=',')

sns_plt = sns.distplot(strand_events)
#sns_plt.set_xlim(0, 200000)
sns_plt.figure.savefig(args["plot"])

# the part for plotting an overlay
#for instance in overlay:
#    plt.plot(np.arange(0,len(instance)), instance, alpha=0.5)

plt.show()