# COMBINE Metagenomics Workshop, 2024

Pawsey, November 11<sup>th</sup> and 12<sup>th</sup>, 2024.


# Workshop Location

The workshop will be held at the [Pawsey Supercomputing Centre](https://maps.app.goo.gl/FCV7nFaRGJLYEvbu8), 1 Bryce Avenue, Kensington WA 6151. The workshop is for registered attendees only.

# Workshop Schedule

_Please note:_ The workshop schedule is subject to change depending on how quickly or slowly we progress.

## Day 1: Monday November 11<sup>th</sup>

Time | Topic 
--- | --- 
0900-0915 | Welcome and Introductions. 
0915-1000 | Introduction to metagenomics, Linux, and Bash 
1000-1030 | Bash hands on and practice: using conda
1030-1100 | Morning coffee
1100-1200 | Downloading data and filtering human genomes using minimap2
1200-1300 | Lunch 
1300-1330 | Introduction to methods for identifying species GTDB, SingleM 
1330-1400 | Hands-on with SingleM 
1400-1500 | Introduction to binning 
1500-1530 | Afternoon tea 
1530-1600 | Microbial binning 
1600-1700 | Hands on with binning


## Day 2: Tuesday November 12<sup>th</sup>

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

Here are some machines that you can use, if you don't have access to a server:

```
IP Addresses:
1: 
2: 
3: 
4: 
5: 
6: 
7: 
```


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

All the data we are going to use in the workshop is present on the servers in `/storage/data/cf_data`

# Use `fastp` to trim bad sequences and remove the adapters.

We are going to use the [Illumina Adapters](https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/SequenceQC/IlluminaAdapters.fa), and trim out:

1. Sequences that are less 100 bp
2. Sequences that contain 1 N
3. Trim the adapters off the 3' and 5' ends of the sequences

Do you remember how to download the Illumina Adapters? The URL is `https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/SequenceQC/IlluminaAdapters.fa`

Once you have downloaded the adapters, we can use this command:

```bash
mkdir fastp
fastp -n 1 -l 100 -i /storage/data/cf_data/reads/788707_20180129_S_R1.fastq.gz -I /storage/data/cf_data/reads/788707_20180129_S_R2.fastq.gz -o fastp/788707_20180129_S_R1.fastq.gz -O fastp/788707_20180129_S_R2.fastq.gz --adapter_fasta IlluminaAdapters.fa
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

We have made this data available for you in `/storage/data/human`

## Install minimap and samtools

This also installs the conda channels for you in the right order!
```
conda config --add channels bioconda
conda config --add channels conda-forge
mamba create -n minimap2 minimap2 samtools
```

## Use minimap2 and samtools to filter the human sequences


```
mkdir -p bam/
minimap2 --split-prefix=tmp$$ -a -xsr /storage/data/human/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz fastp/788707_20180129_S_R1.fastq.gz fastp/788707_20180129_S_R2.fastq.gz | samtools view -bh | samtools sort -o bam/788707_20180129.bam
samtools index bam/788707_20180129.bam
```

Here is the [samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the description of the columns is on page 6.

Now, we use `samtools` flags to filter out the human and not human sequences. You can find out what the flags mean using the [samtools flag explainer](https://broadinstitute.github.io/picard/explain-flags.html)

### human only sequences

```
mkdir human not_human
samtools fastq -F 3588 -f 65 bam/788707_20180129.bam | gzip -c > human/788707_20180129_S_R1.fastq.gz
echo "R2 matching human genome:"
samtools fastq -F 3588 -f 129 bam/788707_20180129.bam | gzip -c > human/788707_20180129_S_R2.fastq.gz
```

### sequences that are not human

```
samtools fastq -F 3584 -f 77 bam/788707_20180129.bam  | gzip -c > not_human/788707_20180129_S_R1.fastq.gz
samtools fastq -F 3584 -f 141 bam/788707_20180129.bam | gzip -c > not_human/788707_20180129_S_R2.fastq.gz
samtools fastq -f 4 -F 1 bam/788707_20180129.bam | gzip -c > not_human/788707_20180129_S_Singletons.fastq.gz
```


### Using snakemake


In the above example, we started with two fastq files, removed adapter sequences, mapped them to the human genome, and then separated out the human and not human sequences.

We can combine all of that into a single `snakemake` file, and it will do all of the steps for us.

See the [Snakemake](../Snakemake) section for details on how to run these two commands in a single pipeline.


## Assembling the sequences

*Note:* Assembling _may_ take a while, and for the workshops, Rob has already assembled the sequences. We may, however, assemble some of them depending on computational resources!

We will assemble with megahit:

```
megahit -1 not_human/788707_20180129_S_R1.fastq.gz -2 not_human/788707_20180129_S_R2.fastq.gz -o megahit_assembled/788707_20180129_S -t 8
```

This generates a contig file called `final.contigs.fa`.

**Questions:**

 - How long is the longest sequence?
 - How long is the shortest sequence?
 - What are the N<sub>50</sub>, N<sub>75</sub>, and auN?


### Counting the lengths of the sequences

We can use [count_fasta.py](bin/count_fasta.c) to count the lengths of the sequences:

```
count_fasta -p -f megahit_assembled/788707_20180129_S/final.contigs.fa | perl -pe 's/^>(\S+).*\t/$1\t/' > megahit_assembled/788707_20180129_S/final.contigs.lengths.tsv
```

## Mapping the reads to the contigs

For this patient, we sequenced several different samples. You can click on the links to get each of the samples

Sample | R1 or R2 | Number of sequences | Total length | Shortest | Longest | N50 | N75 | auN
-- | -- | -- | -- | -- | -- | -- | -- | --
[788707_20171213_S](../Datasets/CF/788707_20171213_S_R1.fastq.gz) | R1 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20171213_S](../Datasets/CF/788707_20171213_S_R2.fastq.gz) | R2 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20180129_S](../Datasets/CF/788707_20180129_S_R1.fastq.gz) | R1 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20180129_S](../Datasets/CF/788707_20180129_S_R2.fastq.gz) | R2 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20180313_S](../Datasets/CF/788707_20180313_S_R2.fastq.gz) | R2 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20180313_S](../Datasets/CF/788707_20180313_S_R1.fastq.gz) | R1 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20181126_S](../Datasets/CF/788707_20181126_S_R1.fastq.gz) | R1 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300
[788707_20181126_S](../Datasets/CF/788707_20181126_S_R2.fastq.gz) | R2 | 125,000 | 37,500,000 | 300 | 300 | 300 | 300 | 300

Download each of the datasets that we don't already have and put them in the `reads` directory. We are going to map all these reads to the contigs we have just created.

We are going to use `minimap`, like we did beore. However, here is a little bit of code that can run `minimap` on all of the samples!

```
mkdir bam_contigs
for R1 in $(find reads/ -name \*R1\* -printf "%f\n"); do 
	R2=${R1/R1/R2}; 
	BAM=${R1/_R1.fastq.gz/.contigs.bam}; 
	minimap2 --split-prefix=tmp$$ -t 8 -a -xsr  megahit_assembled/788707_20180129_S/final.contigs.fa reads/$R1 reads/$R2 | samtools view -bh | samtools sort -o bam_contigs/$BAM;
done
find bam_contigs -type f -exec samtools index {} \;
```

## Generating a depth profile

Now we can generate a depth profile for each contig in each of the four samples. Before we do this, take a look at the output from `samtools coverage` from _one_ `.bam` file!

```
samtools coverage bam_contigs/788707_20171213_S.contigs.bam | less
```

Now we iterate over all the files and get the first column, the contig name, and the 7<sup>th</sup> column which has the mean depth for that contig.

```
for BAM in $(find bam_contigs -type f -name \*bam -printf "%f\n"); do 
	OUT=${BAM/.contigs.bam/.tsv}; 
	samtools coverage bam_contigs/$BAM | cut -f 1,7 > bam_contigs_tsv/$OUT; 
done
```

We have created an [example Jupyter notebook](Workshop_MAG_demo.ipynb) so you can see some of the commands

## Identifying contigs that belong to the same genome

We are going to move the data to [Google Colab](https://colab.research.google.com/) to analyse the data and identify contigs that co-occur across multiple samples.

# Predicting the species that are in the sample

## SingleM

take a [look at the manual](https://wwood.github.io/singlem/) for detailed `singlem` instructions.


Install singleM with mamba:

```
mamba create -n singlem -c bioconda singlem
mamba activate singlem
```

Before you use it, make sure you add this command: 

`export SINGLEM_METAPACKAGE_PATH=/storage/data/metapackage`


Now run singleM on the CF data:

```
singlem pipe -1 /storage/data/cf_data/CF_Data_R1.fastq.gz -2 /storage/data/cf_data/CF_Data_R2.fastq.gz -p output_profile.tsv
```

