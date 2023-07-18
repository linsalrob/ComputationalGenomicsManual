"""

An example snakefile

You may use this snakefile as you wish, but I am not liable if you do
use it!

Rob Edwards, July 2023

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
        expand(os.path.join("fastp", "{sample}_R1.fastq.gz"), sample=SAMPLES)

rule run_fastp:
    input:
        r1 = os.path.join(READDIR, "{sample}_R1." + file_extension),
        r2 = os.path.join(READDIR, "{sample}_R2." + file_extension),
        ad = "IlluminaAdapters.fa"
    output:
        r1 = os.path.join("fastp", "{sample}_R1." + file_extension),
        r2 = os.path.join("prinseq", "{sample}_R2." + file_extension),
        fp = temporary("fastp.html")
        fj = temporary("fastp.json")
    shell:
        """
        fastp fastp -n 1 -l 100 -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --adapter_fasta {input.ad}
        """


