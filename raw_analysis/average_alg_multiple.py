# Implementation of algorithm from https://stackoverflow.com/a/22640362/6029703
import numpy as np
import h5py
import os
import argparse

def thresholding_algo(y, lag):
    min_length = 5000
    signals = np.zeros(len(y))
    blockage = np.zeros(len(y))
    avg = [0]*len(y)
    std = [0]*len(y)
    avg[lag - 1] = np.mean(y[0:lag])
    std[lag - 1] = np.std(y[0:lag])
    current = 0
    first_occurence = 0
    for i in range(lag, len(y)):
        avg[i] = np.mean(y[(i-lag+1):i+1])
        std[i] = np.std(y[(i-lag+1):i+1])
        if abs(avg[i]-avg[i-100]) < 0.015*avg[i]:
            signals[i] = 1
            current += 1
            if current>min_length:
                blockage[i-min_length:i] = 1
                if first_occurence==0:
                    first_occurence = i-min_length
                    print("std:", std[i])
        else:
            current = 0

    return dict(signals = np.asarray(signals),
                blockage = np.asarray(blockage), 
                avg = np.asarray(avg),
                first_occurence = first_occurence)

def locate_blk(fast5):
    lag = 400

    for key in fast5['Raw']['Reads'].keys():
        y = fast5['Raw']['Reads'][key]['Signal'][()]
        y = y.astype(np.int16) 
    result = thresholding_algo(y, lag=lag)

    return result["first_occurence"]

def find_read(name):
    if 'Read_' in name:
        return name


ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True,
    help="path to fast5 files")
ap.add_argument("-o", "--output", required=True,
    help="path to savedir")
args = vars(ap.parse_args())

files = os.listdir(args["input"])
k = 0

for file in files:

    if file.endswith("fast5"):

        fast5 = h5py.File(os.path.join(args["input"], file), 'r')
        k += 1
        CUT = locate_blk(fast5)

        for key in fast5['Raw']['Reads'].keys():
            signal = fast5['Raw']['Reads'][key]['Signal'][()]
            signal = signal.astype(np.int16) 
            key_str = fast5['Raw']['Reads'].visit(find_read)

        if CUT>5000:
            print(file)
            with h5py.File(os.path.join(args["output"], file), "w") as f:

                file_version = fast5['/'].attrs['file_version']
                f['/'].attrs['file_version'] = file_version

                fast5.copy('PreviousReadInfo', f)

                raw = f.create_group('Raw')
                raw_reads = raw.create_group('Reads')
                print(fast5.get('Raw/Reads/'+key_str))
                fast5.copy('Raw/Reads/'+key_str, raw_reads)
                length = f['Raw/Reads/'+key_str].attrs['duration']
                #length = np.uint32(length - CUT)
                length = np.uint32(CUT)
                f['Raw/Reads/'+key_str].attrs['duration'] = length
                del f['Raw/Reads/'+key_str+'/Signal']
                signal_ = signal[:CUT]
                if len(signal_)<201536:
                    chunks = signal_.shape
                else:
                    chunks = (201536,)

                signal_ = raw_reads.create_dataset(key_str+'/Signal', data=signal_, compression="gzip", compression_opts=1, \
                    shape=signal_.shape, dtype=signal_.dtype, maxshape=signal_.shape, chunks=chunks)

                uniq = f.create_group('UniqueGlobalKey')
                fast5.copy('UniqueGlobalKey/channel_id', uniq)
                fast5.copy('UniqueGlobalKey/context_tags', uniq)
                fast5.copy('UniqueGlobalKey/tracking_id', uniq)

        if k>4000:
            break
