# Removing Human (or other) contamination from your metagenome.

If you have sequenced a metagenome, or are analyzing a metagenome, that comes from a human, or is associated with another animal, there is a good chance that you will find DNA from the host contaminating your library. 

Before you analyze the data, you need to remove that contamination so that you don't skew the results.

There are several ways to do this. For example, [deconseq](http://deconseq.sourceforge.net/) is a standalone program that can do this for you. However, if you have a unique or different host, you may need to build your own databases to remove the contamination. Here, we work through the steps of removing contamination. In this example, we use build 38 of the human genome to demonstrate how to do so. If you are using your own reference genome, you can skip the first few steps.

## Download and prepare the human genome

In this example, we're using build 38 of the human genome available from UC Santa Cruz. There are several [different files available](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/) for the human genome, but we are using `hg38.fa.gz`:

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
unpigz hg38.fa.gz
bowtie2-build --threads 16 hg38.fa humanGenome
```


Now we are ready to screen against the human genome!


## Mapping your reads to the human genome

For this example, we're using the data set [SRR2005706](https://www.ncbi.nlm.nih.gov/sra/SRR2005706). We happen to know that metagenome has human contamination in it, and is a relatively small metagenome so it is perfect to test for this example! (Our new version of [partie](https://github.com/linsalrob/partie) includes a measurement of human sequence contamination!)

First, we download the whole dataset, and split it into three files, a left reads, a right reads, and an unpaired reads:

```bash
fastq-dump --outdir fastq --skip-technical  --readids --read-filter pass --dumpbase --clip --split-3 SRR2005706
```

This creates a directory called `fastq` with three files: `fastq/SRR2005706_pass_1.fastq`, `fastq/SRR2005706_pass_2.fastq`, and `fastq/SRR2005706_pass.fastq`. The file that ends `pass_1.fastq` are the *left* reads, the file that ends `pass_2.fastq` are the *right* reads, and the file that ends `pass.fastq` are the unpaired reads.

We use bowtie2 to map this data to our human genome, using the index file that we created above. Note in this example I am using 16 threads, but you may need to change that to match your compute resources as appropriate:

```bash
bowtie2 --threads 16 -1 fastq/SRR2005706_pass_1.fastq -2 fastq/SRR2005706_pass_2.fastq -U fastq/SRR2005706_pass.fastq  -x humanGenome | samtools sort -o SRR2005706.human.bam 
samtools index SRR2005706.human.bam
```

This makes an indexed `.bam` file of the alignments, and reports how many sequences were aligned:

<pre>
42120 reads; of these:
  34823 (82.68%) were paired; of these:
    33665 (96.67%) aligned concordantly 0 times
    974 (2.80%) aligned concordantly exactly 1 time
    184 (0.53%) aligned concordantly >1 times
    ----
    33665 pairs aligned concordantly 0 times; of these:
      404 (1.20%) aligned discordantly 1 time
    ----
    33261 pairs aligned 0 times concordantly or discordantly; of these:
      66522 mates make up the pairs; of these:
        66363 (99.76%) aligned 0 times
        42 (0.06%) aligned exactly 1 time
        117 (0.18%) aligned >1 times
  7297 (17.32%) were unpaired; of these:
    4116 (56.41%) aligned 0 times
    2004 (27.46%) aligned exactly 1 time
    1177 (16.13%) aligned >1 times
8.40% overall alignment rate
</pre>

## Splitting the reads.

Our goal is to use the [samtools flags](https://broadinstitute.github.io/picard/explain-flags.html) to parse the alignment file into reads that map the human genome and reads that do not map the human genome. However, we also want to preserve our left and right pairs so we can do interesting things with them later.

First, we make a set of human-only reads using a combination of two flags. `-G 12` excludes any reads that match having an unmapped read or an unmapped mate. In other words, we only find mapped reads. `-f 65` says find a read that is paired *and* is the first read in the pair (i.e. the *left* read). So together, we are just finding the left reads that *are* mapped to the human genome.

We repeat this with `-f 129` which means the read is paired and it is the second read in the pair.

Finally, we add `-F 5` which is looking for reads that are paired and unmapped, except we invert that match by using a capital F rather than a lower case f, and so we get reads that are NOT paired and are NOT mapped. Note the subtle difference between the `-G` and `-F` options to `samtools`.

Check out the [samtools flag calculator](https://broadinstitute.github.io/picard/explain-flags.html) for more help parsing these flags.


```bash
samtools fastq -G 12 -f 65 SRR2005706.human.bam > is_human/left.human.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1659 reads
samtools fastq -G 12 -f 129 SRR2005706.human.bam > is_human/right.human.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1659 reads
samtools fastq -F 5 SRR2005706.human.bam > is_human/single.human.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 3181 reads
```


And if we take a look at these sequences, we see the length distributions (note that `countfastq.py` is available from [our git repo](https://github.com/linsalrob/EdwardsLab).

```bash
countfastq.py -f is_human/left.human.fastq 
Number of sequences: 1659
Total length: 165563
Shortest: 80
Longest: 100
N50: 100
N75: 100

countfastq.py -f is_human/right.human.fastq 
Number of sequences: 1659
Total length: 165412
Shortest: 81
Longest: 100
N50: 100
N75: 100

countfastq.py -f is_human/single.human.fastq 
Number of sequences: 3181
Total length: 317353
Shortest: 80
Longest: 100
N50: 100
N75: 100
```

We have a total of 1,659 pairs that map to the human genome. Of these, 1,158 were concordantly mapped. If we want just those, we can add a flag 2 to our selection:

```bash
samtools fastq -G 13 -f 131 SRR2005706.human.bam > is_human/right_concordantly_aligned.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 1158 reads

```

and similarly we can identify the discordantly aligned sequences:

```bash
samtools fastq -G 13 -f 129 -F 2 SRR2005706.human.bam > is_human/right_discordantly_aligned.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 501 reads
```


The 3,181 unpaired reads are reported by bowtie2:
<pre>
    2004 (27.46%) aligned exactly 1 time
    1177 (16.13%) aligned >1 times
</pre>


## Extract the reads that are NOT human

We take essentially the same approach to map the non human reads, this time using `samtools` flag `-f 77` to find a read that is paired, where neither the read nor its mate are mapped, and the read is the first in the pair. Similarly, `-f 141` retrieves the second in the pair that is not mapped, and finally `-f 4` finds unmapped, unpaired reads.

```bash
samtools fastq -f 77 SRR2005706.human.bam > not_human/left.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 33164 reads
samtools fastq -f 141 SRR2005706.human.bam > not_human/right.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 33164 reads
samtools fastq -f 4 -F 1 SRR2005706.human.bam > not_human/single.fastq
[M::bam2fq_mainloop] discarded 0 singletons
[M::bam2fq_mainloop] processed 4116 reads
```

### Read counts


Source | Left and right pairs | Singletons |
--- | --- | --- |
Human sequences | 1,659 | 3,181
Non Human sequences | 33,164 | 4,116
Total (human + not human) | 34,823 | 7,297
&nbsp; | &nbsp; | &nbsp;
Original downloaded dataset | 34,823 | 7,297 |









