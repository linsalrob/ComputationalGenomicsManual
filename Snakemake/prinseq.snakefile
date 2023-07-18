"""

An example snakefile written for use on deepthought.

You may use this snakefile as you wish, but I am not liable if you do
use it!

Rob Edwards, September 2020

"""


import os
import sys


# set this to whatever the name of your directory
# with the reads is. If you are following along with the
# tutorial, you can leave this as fastq
READDIR = 'fastq'

# Note that this example requires an R1 file AND an R2 file
# and that each file should match *_R1.* and *_R2.* 
SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1.{extentions}'))

# just get the first file extension as we don't need to iterate all of them
file_extension = EXTENSIONS[0]


# just check there is something to actually do!
if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write(f"Do you have a directory called {READDIR} with some fastq files in it?\n")
    sys.stderr.write("Do those fastq files have _R1. and _R2.?\n")
    sys.exit()


rule all:
    input:
        expand(os.path.join("prinseq", "{sample}_R1.good.fastq"), sample=SAMPLES)

rule run_prinseq:
    input:
        r1 = os.path.join(READDIR, "{sample}_R1." + file_extension),
        r2 = os.path.join(READDIR, "{sample}_R2." + file_extension)
    output:
        r1 = os.path.join("prinseq", "{sample}_R1.good.fastq"),
        s1 = os.path.join("prinseq", "{sample}_R1.single.fastq"),
        b1 = os.path.join("prinseq", "{sample}_R1.bad.fastq"),
        r2 = os.path.join("prinseq", "{sample}_R2.good.fastq"),
        s2 = os.path.join("prinseq", "{sample}_R2.single.fastq"),
        b2 = os.path.join("prinseq", "{sample}_R2.bad.fastq")
    shell:
        """
        prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
            -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
            -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
            -trim_qual_right 30 -trim_qual_window 10 \
            -out_good {output.r1} -out_single {output.s1} -out_bad {output.b1} \
            -out_good2 {output.r2} -out_single2 {output.s2} -out_bad2 {output.b2} \
            -fastq {input.r1} \
            -fastq2 {input.r2};
        """


