#usage: python transform4fastx.py yourFastaq.fastq > reformatted.fastq

import sys

inFile = open(sys.argv[1],'r')

header = ''
seq = ''
qual = ''

seqs = False
quals = False
for line in inFile:
    if line[0] == "@":
        if header != '':
            print("@" + header)
            print(seq.upper())
            print("+" + header)
            print(qual)

        header = line[1:].strip()
        seqs = True
        quals = False
        qual = ''
        seq = ''
    elif line[0] == "+":
        seqs = False
        quals = True
    else:
        if quals:
            qual += line.strip()
        if seqs:
            seq += line.strip()

print("@" + header)
print(seq)
print("+" + header)
print(qual)