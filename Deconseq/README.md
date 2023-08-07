# Removing host genome contamination

We have writen several tools to remove host genomes like [Deconseq](https://deconseq.sourceforge.net/) and written [blog posts](https://edwards.flinders.edu.au/command-line-deconseq/) to describe how to remove host genomes. It's critical that you do this before you analyse your data, because otherwise you will get spurious results when you do your other analyses.

# Step 1. Find your genome

You need to find a host genome that you want to remove. If you are using a non-model organism, take a look at the [NCBI Genomes](https://www.ncbi.nlm.nih.gov/genome) collection that has a lot of curated genomes. 

If you are using the human genome, the NCBI has some specific resources:

There are several [human genome versions](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/) specifically designed for inclusion in pipelines. 

```
A. GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

A gzipped file that contains FASTA format sequences for the following:
1. chromosomes from the GRCh38 Primary Assembly unit.
   Note: the two PAR regions on chrY have been hard-masked with Ns.
   The chromosome Y sequence provided therefore has the same
   coordinates as the GenBank sequence but it is not identical to the
   GenBank sequence. Similarly, duplicate copies of centromeric arrays
   and WGS on chromosomes 5, 14, 19, 21 & 22 have been hard-masked
   with Ns (locations of the unmasked copies are given below).
2. mitochondrial genome from the GRCh38 non-nuclear assembly unit.
3. unlocalized scaffolds from the GRCh38 Primary Assembly unit.
4. unplaced scaffolds from the GRCh38 Primary Assembly unit.
5. Epstein-Barr virus (EBV) sequence
   Note: The EBV sequence is not part of the genome assembly but is
   included in the analysis set as a sink for alignment of reads that
   are often present in sequencing samples.

B. GCA_000001405.15_GRCh38_full_analysis_set.fna.gz

A gzipped file that contains all the same FASTA formatted sequences as
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz, plus:

6. alt-scaffolds from the GRCh38 ALT_REF_LOCI_* assembly units.

C. GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz

A gzipped file that contains all the same FASTA formatted sequences as
GCA_000001405.15_GRCh38_full_analysis_set.fna.gz, plus:

7.  human decoy sequences from hs38d1 (GCA_000786075.2)

D. GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz

A gzipped file that contains all the same FASTA formatted sequences as
GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz, plus:

7.  human decoy sequences from hs38d1 (GCA_000786075.2)
```

We usually use [GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz) which contains everything, but you do need to remember that it contains adenovirus sequences and if you are interested in those they may get removed.

If you want to download it, you can use `wget` to get the file:

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
```


## Use minimap2 and samtools to filter the host sequences

You can use any [HTS mapper](https://academic.oup.com/bioinformatics/article/28/24/3169/245777), like [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [BWA](https://github.com/lh3/bwa), but nowadays we prefer [minimap2](https://github.com/lh3/minimap2) because it is a lot faster.

In this example, we're using the human genome from above. But if you have downloaded a different genome, just replace your fasta file. Note that another advantage of `minimap2` is that it handles `gzip` compressed files natively, you don't need to uncompress them!

```
minimap2 --split-prefix=tmp$$ -a -xsr GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz R1.fastq.gz R2.fastq.gz | samtools view -bh | samtools sort -o output.bam
samtools index output.bam
```

## Filter mapped reads

Now, we use `samtools` flags to filter out the host and non host sequences. The reads that *ARE* mapped are host, while the reads that ARE *NOT* mapped are non-host.

We use samtools to get fastq output from the `.bam` format files.

You can find out what the flags mean using the [samtools flag explainer](https://broadinstitute.github.io/picard/explain-flags.html)

Here is the [samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the description of the columns is on page 6.

### host sequences

```
mkdir host not_host
samtools fastq -F 3588 -f 65 output.bam | gzip -c > host/output_S_R1.fastq.gz
echo "R2 matching host genome:"
samtools fastq -F 3588 -f 129 output.bam | gzip -c > host/output_S_R2.fastq.gz
```

### sequences that are not host

```
samtools fastq -F 3584 -f 77 output.bam  | gzip -c > not_host/output_S_R1.fastq.gz
samtools fastq -F 3584 -f 141 output.bam | gzip -c > not_host/output_S_R2.fastq.gz
samtools fastq -f 4 -F 1 output.bam | gzip -c > not_host/output_S_Singletons.fastq.gz
```




