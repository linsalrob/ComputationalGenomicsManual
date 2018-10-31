import sys
from roblib import sequences


barcodes = set(['CTCGATG', 'ACCAACT', 'TCGCAGG', 'GCTCGAA', 'GGATCAA'])

barF=open('barcodes.fastq', 'w')
seqF=open('sequences.fastq', 'w')
count={}

# Maximum number of sequences to write per tag
#maxper = 10000000000000000
maxper = 100000

files = ['SRR2080423_pass.fastq.gz', 'SRR2080425_pass.fastq.gz', 'SRR2080427_pass.fastq.gz', 'SRR2080434_pass.fastq.gz', 'SRR2080436_pass.fastq.gz']

for f in files:
    for (sid, label, seq, qual) in sequences.stream_fastq(f):
        bc = seq[0:7]
        if bc in barcodes:
            count[bc] = count.get(bc, 0)+1
            if count[bc] <= maxper:
                barF.write("@{}\n{}\n+\n{}\n".format(label, seq[0:7], qual[0:7]))
                seqF.write("@{}\n{}\n+\n{}\n".format(label, seq[7:], qual[7:]))

for c in count:
    print("{}\t{}\n".format(count[c], c))
