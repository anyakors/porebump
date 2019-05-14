import numpy as np
from itertools import islice

#AUCG 
word = 'AUGCUCAGAG'

def rolling_window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

#word = 'GGUGGG'
def random_alt(word = 'GGAUAAAG'): 
	words = []
	bases = 'AUCG'
	for i in np.arange(len(word)):
		for base in bases:
			if word[i]!=base:
				words.append(word[:i]+base+word[i+1:])
	return words

#fastq = np.genfromtxt('/home/mookse/workspace/MinKNOW/data/reads/20180614_0352_RNA9_mut/fast5/0/bc/all.fastq', dtype='str', delimiter='\n')
#fastq = np.genfromtxt('/home/mookse/workspace/minion_scripts/quad_train/cut_20180614_0352_RNA9_mut/bc/all.fastq', dtype='str', delimiter='\n')
#fastq = np.genfromtxt('/home/mookse/workspace/MinKNOW/data/reads/20180612_0522_RNA3_G4_2/fast5/0/bc/all.fastq', dtype='str', delimiter='\n')
fastq = np.genfromtxt('/home/mookse/workspace/minion_scripts/quad_train/manual_cut_RNA3_0522/bc/all.fastq', dtype='str', delimiter='\n')
words = random_alt()
words.append(word)
match_no = 0

print(len(fastq))

for k in np.arange(1, len(fastq), 4):
	#print(fastq[k])
	#fastq[1], [5], [9] -> 1+4*i	
	for word in words:
		#print("WORD:", word)
		split = rolling_window(fastq[k], len(word))
		i = 0
		for x in split:
			i += 1
			#print("current X:", ''.join(x))
			indicator = (''.join(x) == word)
			#print(indicator)
			if indicator:
				match_no += 1
				print("MATCH! {}, word: {}".format(i,word))

print("Total matches:", match_no)