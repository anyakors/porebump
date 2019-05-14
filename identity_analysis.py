import os
import h5py
import numpy as np
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

#blast putput -6 as input '6 qseqid sseqid pident qlen length mismatch gapopen gaps qstart qend sstart send evalue bitscore'
ap = argparse.ArgumentParser()
ap.add_argument('file', type=argparse.FileType('r'), nargs='+')
args = ap.parse_args()

pident_all = []

for f in args.file:
	print(f)
	pident, length, aln_len, mismatch, gaps, evalue = np.loadtxt(f, usecols=(2,3,4,5,7,12), unpack=True)
	pident_all.append(pident)

print(len(pident))
print(pident)

total_len = np.sum(length)
total_aln = np.sum(aln_len)
average_aln = np.average(aln_len)
total_mismatch = np.sum(mismatch)
total_gaps = np.sum(gaps)

print("total len - {}, cumulative aln - {}, average aln - {}, total mismatch N - {}, mismatch/100bp - {}, total gaps - {}, gaps/100bp - {}".format(total_len, \
	total_aln, average_aln, total_mismatch, total_mismatch*100/total_aln, total_gaps, total_gaps*100/total_aln))

#sns_plt = sns.distplot(data, rug=True, bins=25)
plt.figure(figsize=(8,6))

for stack in pident_all:
	n, bins, patches = plt.hist(x=stack, bins='auto', 
	                            alpha=0.6, rwidth=0.85)
#sns_plt.set_xlim(min(data), max(data))
#sns_plt.set(xlabel='Identity percent', ylabel='Number of reads')
#sns_plt.figure.savefig("violin_dist.png")
plt.xlabel('Identity percent')
plt.ylabel('Number of reads')
plt.show()