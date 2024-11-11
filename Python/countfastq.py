#!/usr/bin/env python3

"""
Count the characters in a fasta file. We summarize the longest and shortest reads
and the N50 of the data set.
"""

import os
import sys
import gzip
import argparse

__author__ = 'Rob Edwards'


def stream_fastq(fqfile):
    """Read a fastq file and provide an iterable of the sequence ID, the
    full header, the sequence, and the quaity scores.

    Note that the sequence ID is the header up until the first space,
    while the header is the whole header.
    """

    if fqfile.endswith('.gz'):
        qin = gzip.open(fqfile, 'rt')
    else:
        qin = open(fqfile, 'r')

    linecounter = 0
    while True:
        header = qin.readline()
        linecounter += 1
        if not header:
            break
        if not header.startswith("@"):
            print(f"The file {fqfile} does not appear to be a four-line fastq file at line {linecounter}", file=sys.stderr)
            sys.exit(-1)
        header = header.strip()
        seqidparts = header.split(' ')
        seqid = seqidparts[0]
        seqid = seqid.replace('@', '')
        seq = qin.readline().strip()
        linecounter += 1
        qualheader = qin.readline()
        if not qualheader.startswith("+"):
            print(f"The file does not appear to be a four-line fastq file at line {linecounter}", file=sys.stderr)
            sys.exit(-1)
        linecounter += 1
        qualscores = qin.readline().strip()
        linecounter += 1
        header = header.replace('@', '', 1)
        if len(qualscores) != len(seq):
            print(f"The sequence and qual scores are not the same length at line {linecounter}", file=sys.stderr)
            sys.exit(-1)
        yield seqid, header, seq, qualscores



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', nargs='+', help='fastq file')
    parser.add_argument('-d', nargs='+', help='directory of fastq files')
    parser.add_argument('-t', help='tab separated summmary of name, total len, shortest, longest, n50, n75', action="store_true")
    parser.add_argument('-s', help="summarize the counts. Useful if you run on a directory", action="store_true")
    args = parser.parse_args()

    if not args.f and not args.d:
        sys.stderr.write(f"FATAL: Please specify either -f or -d or use -h for more help\n")
        sys.exit(1)

    if args.f:
        files = args.f
    else:
        files = []

    if args.d:
        for subdir in args.d:
            for f in os.listdir(subdir):
                if 'fastq' in f or 'fq' in f:
                    files.append(os.path.join(subdir, f))
                else:
                    sys.stderr.write(f"Skipped {os.path.join(subdir, f)}. Does not appear to be fastq\n")

    overall = {'number': 0, 'total': 0, 'shortest':1e6, 'longest': 0}
    for faf in files:
        if not os.path.exists(faf):
            sys.stderr.write(f"FATAL: {faf} not found\n")
            sys.exit(1)

        lens = []
        for (sid, label, seq, qual) in stream_fastq(faf):
            lens.append(len(seq))
        lens.sort()
        length=sum(lens)

        len_so_far = 0
        n50 = None
        n75 = None
        auN = 0
        for i in lens:
            len_so_far += i
            if not n50 and len_so_far >= length * 0.5:
                n50 = i
            if not n75 and len_so_far >= length * 0.75:
                n75 = i
            auN += i**2

        auN /= length

        if args.t:
            print(f"{faf}\t{len(lens):,}\t{length:,}\t{lens[0]:,}\t" \
                  + f"{lens[-1]:,}\t{n50:,}\t{n75:,}\t{int(auN):,}")
        else:
            print(f"""
File name: {faf}
Number of sequences: {len(lens):,}
Total length: {length:,}
Shortest: {lens[0]:,}
Longest: {lens[-1]:,}
N50: {n50:,}
N75: {n75:,}
auN: {int(auN):,}  """
            )
        overall['number'] += len(lens)
        overall['total']  += length
        if lens[0] < overall['shortest']:
            overall['shortest'] = lens[0]
        if lens[-1] > overall['longest']:
            overall['longest'] = lens[-1]

if args.s:
    print(f"""

OVERALL SUMMARY
Number of sequences: {overall['number']:,}
Total length: {overall['total']:,}
Shortest: {overall['shortest']:,}
Longest: {overall['longest']:,}
"""
      )



