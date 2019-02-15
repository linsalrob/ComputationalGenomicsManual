# Removing Human (or other) contamination from your metagenome.

If you have sequenced a metagenome, or are analyzing a metagenome, that comes from a human, or is associated with another animal, there is a good chance that you will find DNA from the host contaminating your library. 

Before you analyze the data, you need to remove that contamination so that you don't skew the results.

There are several ways to do this. For example, [deconseq](http://deconseq.sourceforge.net/) is a standalone program that can do this for you. However, if you have a unique or different host, you may need to build your own databases to remove the contamination. Here, we work through the steps of removing contamination. In this example, we use build 38 of the human genome to demonstrate how to do so. If you are using your own reference genome, you can skip the first few steps.

## Download and prepare the human genome

In this example, we're using build 38 of the human genome available from UC Santa Cruz. There are several [different files available](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/) for the human genome, but we are using:

<pre>
hg38.fa.gz - "Soft-masked" assembly sequence in one file.
    Repeats from RepeatMasker and Tandem Repeats Finder (with period
    of 12 or less) are shown in lower case; non-repeating sequence is
    shown in upper case.
</pre>
  
First, we download this file, and verify the download worked correctly:


```bash
curl -lo hg38.fa.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt | grep hg38.fa.gz > md5sum.txt
  md5sum --check md5sum.txt
```

At this point you should see the message:
```
hg38.fa.gz  OK
```

If you do not, your download had an issue!

Next, we uncompress the human genome and build a bowtie index. Here, we are using [pigz](https://zlib.net/pigz/) but you can also use `gunzip` to get the same result (just slower!).

Finally, we use `bowtie2-build` to create an index of that file (note, here I am using 16 threads, but you may need more or less!):

```bash
bowtie2-build --threads 16 hg38.fa humanGenome
```


Now we are ready to screen against the human genome!


## Mapping your reads to the human genome

Lets assume your sequences are in three files

left.fastq – the left reads in a fastq file
right.fastq – the right reads in a file
single.fastq – the singleton reads in a fastq file

run the mapping (using 16 processes):

```
bowtie2 -p 16 -x sharks -1 left.fastq -2 right.fastq -U single.fastq | samtools sort > seqs.sharks.bam
```

## Using samtools to separate contamination and not contamination




 # build the index
bowtie2-build sharks.fasta sharks

Now we have our bamfile we can use samtools to separate the sequences into those sequences that map to the sharks and those sequences that do not map to the sharks. We’ll do this in two separate steps.

Note that for all these steps, we are making use of the samtools flags filters, using either -f or -F (match a flag or do not match a flag) or -G (exclude reads that match). You can learn more about the samtools flags and explore the options with this samtools flags calculator.

Sequences that map to the reference
We can extract all the sequences that map to the shark reference genomes into three files: the left reads, the right reads, and the unmapped reads.

# extract the left reads that match the reference genome
samtools fastq -G 12 -f 65 seqs.sharks.bam > left.shark.fastq
# extract the right reads that match the reference genome
samtools fastq -G 12 -f 129 seqs.sharks.bam > right.shark.fastq
# extract the single reads that match the reference genome
samtools fastq -F 5 seqs.sharks.bam > single.shark.fastq
In the first two commands, we use -G 12 which excludes reads where the read is unmapped and the mate is unmapped, and then we look for reads that are paired and it is the first read in the pair (-f 65) or it is the second read in the pair (-f 129). This combination gives us either the left or right mapped reads.

Next, we look for reads that are not pairs (that would normally be -F 1) and are mapped (that would normally be -F 4), hence we use -F 5 to find the singletons that are mapped to the reference, and are hence shark sequences.

Sequences that do not map to the reference
We apply a similar logic to find the sequences that do not map to the reference genome. Again, we’ll end up with three files, a left reads, right reads, and unmapped reads.

# extract the left reads that do NOT match the reference genome
samtools fastq -f 77 seqs.sharks.bam > left.notshark.fastq
# extract the right reads that do NOT match the reference genome
samtools fastq -f 141 seqs.sharks.bam > right.notshark.fastq
# extract the single reads that do NOT match the reference genome
samtools fastq -f 4 -F 1 seqs.sharks.bam > single.notshark.fastq
The logic here is very similar. A samtools flag of -f 77 means the read is paird, neither the read nor mate is mapped, and the it is the first read in the pair, while a flag of -f 141 means the same except it is the second mate in the pair.

Finally, we use two flags: -f 4 (the read is unmapped) and -F 1 (the read is not paired) to find the single reads that are not mapped.

By combining the flags in the samtools files we can easily separate our sequences into those reads that match the reference genome, and those that do not. The next question is, which ones are you interested in!

