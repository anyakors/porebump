import numpy as np

def thresholding_algo(y, lag, threshold, influence):
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)-1):
        if abs(y[i] - avgFilter[i-1]) > threshold*stdFilter[i-1] and abs(y[i+1] - avgFilter[i-1]) > threshold*stdFilter[i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence)*filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))


def peak_frequency(signals):
    k = 0
    last = 0
    periods = []
    frequency = 20
    for i in range(0,len(signals)):
        if signals[i]!=last:
            periods.append(i-k)
            k = i
            last = signals[i]
    if len(periods)>5:
        frequency = np.mean(periods)
    return frequency


def drop_occurence(signals, freq_threshold):
    k = 0 #last change
    last = 0 #last value
    drops = [0]*len(signals)
    for i in range(0,len(signals)):
        if signals[i]==last:
            current = i-k #length of constant value
            if current>freq_threshold:
                drops[i-current:i] = [1]*current
        else:
            current = 0
            last = signals[i]
            k = i
    return drops


def continuous_drop(array):
    k = np.argmax(array) #last change
    drop_number = 0
    current = 0
    for i in range(0,len(array)):
        if array[i]==1:
            current += 1 #length of constant value
            if current==1:
                drop_number += 1
        else:
            current = 0
    return drop_number


def first_drop(array):
    k = np.argmax(array) #last change
    current = 0
    drop_index = 0
    for i in range(0,len(array)):
        if array[i]==1:
            drop_index = i
            break
    return drop_index


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()