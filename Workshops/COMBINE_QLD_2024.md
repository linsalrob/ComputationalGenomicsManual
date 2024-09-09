# COMBINE Metagenomics Workshop, 2024

University of Queensland, Queensland, Australia, 2024


# Workshop Location

Room 3.142 [Institute for Molecular Bioscience](https://imb.uq.edu.au/ ), University of Queensland, St Lucia, QLD.

You can find the room here:

![Room 3.142 is behind this secret door!](images/metagenomics_map_uq_2024.png)

# Workshop Schedule

_Please note:_ The workshop schedule is subject to change depending on how quickly or slowly we progress.

## Day 1: Monday September 9<sup>th</sup>

Time | Topic 
--- | --- 
0900-0915 | Welcome and Introductions. 
0915-1000 | [Introduction to metagenomics, Linux, and Bash](https://linsalrob.github.io/ComputationalGenomicsManual/Linux/)
1000-1030 | Bash hands on and practice
1030-1100 | Morning coffee
1100-1200 | Introduction to methods for identifying species GTDB, SingleM 
1200-1330 | Lunch 
1330-1400 | Hands-on with SingleM 
1400-1500 | Introduction to binning 
1500-1530 | Afternoon tea 
1530-1600 | Microbial binning 
1600-1700 | Invited Presentation:  Professor Elizabeth Dinsdale


## Day 2: Tuesday September 10<sup>th</sup>

Time | Topic 
-- | --
0900-0930 | Recap of Day 1 
0930-1000 | Introduction to methods for functional analysis 
1000-1030 | SUPER-FOCUS hands-on 
1030-1100 | Morning coffee 
1100-1200 | Viral identification using Hecatomb 
1200-1230 | Hecatomb hands-on 
1230-1330 | Lunch 
1400-1500 | Hands-on data visualisation 
1500-1530 | Afternoon tea 
1530-1645 | Hands-on data visualisation 
1645-1700 | Wrap up and summary 


# Download

If you are using a MS Windows machine, please download and install [MobaXterm](https://mobaxterm.mobatek.net/) before we start. If you are using a Mac, you are good to go!

# Metagenomics

We are going to jump right in with metagenomics, but [here is a brief introduction](https://linsalrob.github.io/ComputationalGenomicsManual/Metagenomics/) if you want to read something while Rob is talking.


We have created servers for you with all the software and data that you will need for these excercises. 

There are two machines that you can use, if you don't have access to a server:

```
IP Addresses:
1: 35.184.190.79
2: 34.138.12.48
3: 
4: 
5: 
```


# Usernames and passwords

Access them [here](https://docs.google.com/spreadsheets/d/1RDjBFT1KdO_G_2bJw-GZv5NJDeGVclAy22Cftg_ujyk/edit?usp=sharing)

# Learning BASH

If you need more helo with bash, you can follow the Pony!

We have installed [PonyLinux](https://github.com/NCGAS/PonyLinux). You only need two commands to get started.

First, log into the accounts provided for you and then follow the commands:

Next type `cd PonyLinux` and press `enter` (or on some computers it is called `return`).

Now, type `./Ponylinux.sh` and press `enter` (or `return`).



# Installing  software

Our first excercise is installing software using mamba. 


Install `conda`, `fastp`, `minimap2`, `samtools` using [conda](../Conda/)

We will use all of these programs today.


# Downloading Data

Download the [CF data](../Datasets/CF) fastq files, or if you want, you can use your own fastq files and see what you find!

For this example, I will download:

- [R1 file](https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/Datasets/CF/788707_20180129_S_R1.fastq.gz)
- [R2 file](https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/Datasets/CF/788707_20180129_S_R2.fastq.gz)

If you are using a remote server, you can use `wget` again to download these reads:


```
wget https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/Datasets/CF/788707_20180129_S_R1.fastq.gz
wget https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/Datasets/CF/788707_20180129_S_R2.fastq.gz
```

# Use `fastp` to trim bad sequences and remove the adapters.

We are going to use the [Illumina Adapters](https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/SequenceQC/IlluminaAdapters.fa), and trim out:

1. Sequences that are less 100 bp
2. Sequences that contain 1 N
3. Trim the adapters off the 3' and 5' ends of the sequences

Do you remember how to download the Illumina Adapters? The URL is `https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/SequenceQC/IlluminaAdapters.fa`

Once you have downloaded the adapters, we can use this command:

```bash
mkdir fastp
fastp -n 1 -l 100 -i 788707_20180129_S_R1.fastq.gz -I 788707_20180129_S_R2.fastq.gz -o fastp/788707_20180129_S_R1.fastq.gz -O fastp/788707_20180129_S_R2.fastq.gz --adapter_fasta IlluminaAdapters.fa
```

When `fastp` runs, you will get an HTML output file called [fastp.html](fastp_788707_20180129.html). This shows some statistics about the run.


# Filter out the human and non-human sequences.

The [NCBI](http://www.ncbi.nlm.nih.gov/) has several [human genome versions](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/) specifically designed for inclusion in pipelines like this. 

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

For this work, we are going to use [GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz) which contains everything.

If you want to download it, you can use `wget` like we have done before!

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
```

## Use minimap2 and samtools to filter the human sequences


```
minimap2 --split-prefix=tmp$$ -a -xsr GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz fastp/788707_20180129_S_R1.fastq.gz fastp/788707_20180129_S_R2.fastq.gz | samtools view -bh | samtools sort -o 788707_20180129.bam
samtools index 788707_20180129.bam
```

Here is the [samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the description of the columns is on page 6.

Now, we use `samtools` flags to filter out the human and not human sequences. You can find out what the flags mean using the [samtools flag explainer](https://broadinstitute.github.io/picard/explain-flags.html)

### human only sequences

```
mkdir human not_human
samtools fastq -F 3588 -f 65 788707_20180129.bam | gzip -c > human/788707_20180129_S_R1.fastq.gz
echo "R2 matching human genome:"
samtools fastq -F 3588 -f 129 788707_20180129.bam | gzip -c > human/788707_20180129_S_R2.fastq.gz
```

### sequences that are not human

```
samtools fastq -F 3584 -f 77 788707_20180129.bam  | gzip -c > not_human/788707_20180129_S_R1.fastq.gz
samtools fastq -F 3584 -f 141 788707_20180129.bam | gzip -c > not_human/788707_20180129_S_R2.fastq.gz
samtools fastq -f 4 -F 1 788707_20180129.bam | gzip -c > not_human/788707_20180129_S_Singletons.fastq.gz
```


# Using snakemake

In the above example, we started with two fastq files, removed adapter sequences, mapped them to the human genome, and then separated out the human and not human sequences.

We can combine all of that into a single `snakemake` file, and it will do all of the steps for us.

See the [Snakemake](../Snakemake) section for details on how to run these two commands in a single pipeline.



