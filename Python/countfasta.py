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


def read_fasta(fname: str, whole_id: bool = False, qual: bool = False) -> dict:
    """
    Read a fasta file and return a hash.

    If wholeId is set to false only the first part of the ID
    (upto the first white space) is returned

    :param fname: The file name to read
    :param whole_id: Whether to keep the whole id, or trim to first whitespace (default = all)
    :param qual: these are quality scores (so add a space between lines!)
    :return: dict
    """

    try:
        if fname.endswith('.gz'):
            f = gzip.open(fname, 'rt')
        elif fname.endswith('.lrz'):
            f = subprocess.Popen(['/usr/bin/lrunzip', '-q', '-d', '-f', '-o-', fname], stdout=subprocess.PIPE).stdout
        else:
            f = open(fname, 'r')
    except IOError as e:
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write("Message: \n" + str(e.message) + "\n")
        sys.exit("Unable to open file " + fname)

    seqs = {}
    seq = ""
    seqid = ""
    for line in f:
        line = line.rstrip('\r\n')
        if line.startswith(">"):
            if seqid != "":
                seqs[seqid] = seq
                seq = ""
            seqid = line.replace(">", "", 1)
            if not whole_id and seqid.count(" ") > 0:
                seqids = seqid.split(" ")
                seqid = seqids[0]
        else:
            if qual:
                seq += " " + line
            else:
                seq += line

    seqs[seqid] = seq.strip()
    return seqs



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-f', nargs='+', help='fasta file')
    parser.add_argument('-d', nargs='+', help='directory of fasta files')
    parser.add_argument('-l', help='list the lengths for each sequence (default = not to)', action='store_true')
    parser.add_argument('-m', help='minimum length fo be inclued', type=int, default=0)
    parser.add_argument('-n', help='do NOT print the summary at the end', action='store_true')
    parser.add_argument('-t', help='tab separated output. Fields: [SeqID, # seqs, total bp, shortest, longest, N50, N75, auN]', action='store_true')
    parser.add_argument('-v', help='more output', action='store_true')
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
                files.append(os.path.join(subdir, f))
                if args.v:
                    sys.stderr.write(f"Added {files[-1]}\n")

    overall = {'number': 0, 'total': 0, 'shortest':1e6, 'longest': 0}


    for faf in files:
        if args.v:
            sys.stderr.write(f"Counting sequences in {faf}\n")
        fa=read_fasta(faf)

        if len(fa.keys()) == 1 and list(fa.keys())[0] == '':
            sys.stderr.write(f"No sequences found in {faf}\n")
            sys.exit(0)

        if args.l:
            for i in fa:
                print("{}\t{}".format(i, len(fa[i])))
            print()

        lensall=[len(fa[i]) for i in fa]
        lens = list(filter(lambda x: x > args.m, lensall))
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

        if args.n:
            continue

        if args.t:
            print("\t".join(map(str, [faf, len(lens), length, lens[0], lens[-1], n50, n75, auN])))
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

