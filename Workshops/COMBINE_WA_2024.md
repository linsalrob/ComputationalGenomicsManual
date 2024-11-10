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
1300-1330 | Introduction to methods for identifying species focus, SingleM 
1330-1500 | Hands-on with SingleM and focus
1500-1530 | Afternoon tea 
1530-1600 | Introduction to methods for functional analysis
1600-1700 | SUPER-FOCUS hands-on


## Day 2: Tuesday November 12<sup>th</sup>

Time | Topic 
-- | --
0900-0930 | Recap of Day 1 
0930-1000 | Introduction to binning
1000-1030 | Megahit assembly
1030-1100 | Morning coffee 
1100-1200 | Binning by hand
1200-1230 | Viral identification using Hecatomb 
1230-1330 | Lunch 
1400-1500 | Hands-on data visualisation 
1500-1530 | Afternoon tea 
1530-1645 | Hands-on data visualisation 
1645-1700 | Wrap up and summary 

# Download

If you are using a MS Windows machine, please download and install [MobaXterm](https://mobaxterm.mobatek.net/) before we start. If you are using a Mac, you are good to go!



# Metagenomics

We are going to jump right in with metagenomics, but [here is a brief introduction](https://linsalrob.github.io/ComputationalGenomicsManual/Metagenomics/) if you want to read something while Rob is talking.


We have created accounts for you on pawsey, and we will share the usernames and passwords with you at the workshop. *These are temporary accounts and will be deleted at the end of the workshop*




# Installing  software

Our first excercise is installing software using mamba. 

Before we begin, we are going to make lives slightly easier for ourselves by making an `alias` or `symbolic link`:

```
ln -s /software/projects/courses01/$USER software
```

This will create a directory called software.

*Important*: When you install mamba, it will ask you for a location. Use `/home/$USER/software/miniforge3` as the location

Install `conda`, `fastp`, `minimap2`, `samtools` using [conda](../Conda/)

We will use all of these programs today.

You can check that they installed by using the command:

```
which fastp
```

If that works, it will tell you!


# Downloading Data

We are going to use the CF data that Rob talked about. To start we are just going to download two files, an R1 and an R2 file to work with:

```
mkdir fastq
cd fastq
curl -LO https://github.com/linsalrob/ComputationalGenomicsManual/raw/refs/heads/master/Datasets/CF/788707_20180129_S_R1.fastq.gz
curl -LO https://github.com/linsalrob/ComputationalGenomicsManual/raw/refs/heads/master/Datasets/CF/788707_20180129_S_R2.fastq.gz
cd
ls
```

# Use `fastp` to trim bad sequences and remove the adapters.

We are going to use the [Illumina Adapters](https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/SequenceQC/IlluminaAdapters.fa), and trim out:

1. Sequences that are less 100 bp
2. Sequences that contain 1 N
3. Trim the adapters off the 3' and 5' ends of the sequences

Do you remember how to download the Illumina Adapters? The URL is `https://github.com/linsalrob/ComputationalGenomicsManual/raw/master/SequenceQC/IlluminaAdapters.fa`

Once you have downloaded the adapters, we are going to make a slurm script to run the command on the cluster

Use `nano` (or `vi` or `emacs`) to edit a file, and copy this text:

```bash
#!/bin/bash
#SBATCH --job-name=fastp
#SBATCH -o fastp-%j.out
#SBATCH -e fastp-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB

mkdir fastp
fastp -n 1 -l 100 -i fastq/788707_20180129_S_R1.fastq.gz -I fastq/788707_20180129_S_R2.fastq.gz -o fastp/788707_20180129_S_R1.fastq.gz -O fastp/788707_20180129_S_R2.fastq.gz --adapter_fasta IlluminaAdapters.fa
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


## Use minimap2 and samtools to filter the human sequences

We are going to make another slurm script, called `mapping.slurm`:

```
nano mapping.slurm
```

And copy and paste these contents:

```
#!/bin/bash
#SBATCH --job-name=mapping
#SBATCH -o mapping-%j.out
#SBATCH -e mapping-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB


mkdir -p bam/
minimap2 --split-prefix=tmp$$ -t 16 -a -xsr /scratch/courses01/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz fastp/788707_20180129_S_R1.fastq.gz fastp/788707_20180129_S_R2.fastq.gz \
	| samtools view -bh | samtools sort -o bam/788707_20180129.bam
samtools index bam/788707_20180129.bam
```

Here is the [samtools specification](https://samtools.github.io/hts-specs/SAMv1.pdf), and the description of the columns is on page 6.

Now, we use `samtools` flags to filter out the human and not human sequences. You can find out what the flags mean using the [samtools flag explainer](https://broadinstitute.github.io/picard/explain-flags.html)


Again, edit a file, this time called `filtering.slurm`

```
nano filtering.slurm
```

And paste these contents:

```
#!/bin/bash
#SBATCH --job-name=filtering
#SBATCH -o filtering-%j.out
#SBATCH -e filtering-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

mkdir human not_human
samtools fastq -F 3588 -f 65 bam/788707_20180129.bam | gzip -c > human/788707_20180129_S_R1.fastq.gz
samtools fastq -F 3588 -f 129 bam/788707_20180129.bam | gzip -c > human/788707_20180129_S_R2.fastq.gz

# sequences that are not human

samtools fastq -F 3584 -f 77 bam/788707_20180129.bam  | gzip -c > not_human/788707_20180129_S_R1.fastq.gz
samtools fastq -F 3584 -f 141 bam/788707_20180129.bam | gzip -c > not_human/788707_20180129_S_R2.fastq.gz
```


### Using snakemake


In the above example, we started with two fastq files, removed adapter sequences, mapped them to the human genome, and then separated out the human and not human sequences.

We can combine all of that into a single `snakemake` file, and it will do all of the steps for us.

See the [Snakemake](../Snakemake) section for details on how to run these two commands in a single pipeline.


# Read based annotations

In metagenomics, there are two fundemental approaches: read-based annotations and assembly based approaches. We are going to start with read based annotations.



# Predicting the species that are in the sample

## SingleM

take a [look at the manual](https://wwood.github.io/singlem/) for detailed `singlem` instructions.


Install singleM with mamba. _Note:_ Here, we introduce named `mamba` environments. What is the advantage of creating a named environment?

```
mamba create -n singlem -c bioconda singlem
mamba activate singlem
```

After installing singlem, you will get a warning from krona. DO NOT run the `ktUpdate.sh` script. Instead, create a new symlink like so:

```
rm -rf /software/projects/courses01/$USER/miniforge3/envs/singlem/opt/krona/taxonomy
ln -s /scratch/courses01/krona/taxonomy /software/projects/courses01/$USER/miniforge3/envs/singlem/opt/krona/taxonomy
```

Next, before you use singlem, make sure you add this command: 

```
export SINGLEM_METAPACKAGE_PATH='/scratch/courses01/singlem/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'
```


Now run singleM on the CF data:

```
#!/bin/bash
#SBATCH --job-name=singlem
#SBATCH -o singlem-%j.out
#SBATCH -e singlem-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

eval "$(conda shell.bash hook)"
conda activate singlem
export SINGLEM_METAPACKAGE_PATH='/scratch/courses01/singlem/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'
singlem pipe -1 not_human/788707_20180129_S_R1.fastq.gz -2 not_human/788707_20180129_S_R2.fastq.gz -p output_profile.tsv --taxonomic-profile-krona krona.html --threads 16
```

## Unzip the data

Note: Before we carry on, both `focus` and `super-focus` require that we unzip the data.

```
cd not_human
find . -name \*gz -exec gunzip {} \;
cd ..
```


## Focus

Another way to identify the species present is to use [FOCUS](https://github.com/metageni/FOCUS)

We create a mamba environment just for focus:

```
mamba create -n focus -c bioconda focus
```

Now, we need to unpack the database. Here is a trick, since we don't know exactly where the database is:

```
FOCUS=$(find software/miniforge3/envs/ -name db.zip -printf "%h\n")
unzip $FOCUS/db.zip -d $FOCUS
```

This should create the directory `software/miniforge3/envs/focus/lib/python3.13/site-packages/focus_app/db/` with two files inside of it.

Now we can run focus on our data:

```
#!/bin/bash
#SBATCH --job-name=focus
#SBATCH -o focus-%j.out
#SBATCH -e focus-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

eval "$(conda shell.bash hook)"
conda activate focus

focus -q not_human/ -o focus -t 16
```


# SUPER-FOCUS

We are going to assess the functions using [SUPER-FOCUS](https://github.com/metageni/SUPER-FOCUS)

We are going to make _another_ mamba environment for super-focus:

```
mamba create -n superfocus -c bioconda super-focus mmseqs2
```

Now we can run super-focus on our data. _Note_: Superfocus creates a _lot_ of data, and you will likely get an error if you just output the results to your home directory. In this command, we put the results somewhere else!

```
#!/bin/bash
#SBATCH --job-name=superfocus
#SBATCH -o superfocus-%j.out
#SBATCH -e superfocus-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

eval "$(conda shell.bash hook)"
conda activate superfocus

export SUPERFOCUS_DB=/scratch/courses01/superfocus/
superfocus -q not_human/ -dir /scratch/courses01/$USER/superfocus -a mmseqs -t 16 -db DB_95
```

## Recompressing the files.

Now that we are done with `focus` and `superfocus`, we can recompress the files. _Question:_ Why should we compress the files (or not compress them)?

```
cd not_human
find . -type f -exec gzip {} \;
cd ..
```

## Assembling the sequences

*Note:* Assembling _may_ take a while, and for the workshops, Rob has already assembled the sequences. We may, however, assemble some of them depending on computational resources!

We will assemble with megahit.

We need to create a mamba environment for megahit ... how are you going to do that?


```
#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH -o megahit-%j.out
#SBATCH -e megahit-%j.err
#SBATCH --account=courses01
#SBATCH --time=1-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB

eval "$(conda shell.bash hook)"
conda activate megahit

mkdir -p megahit_assembled/
megahit -1 not_human/788707_20180129_S_R1.fastq.gz -2 not_human/788707_20180129_S_R2.fastq.gz -o megahit_assembled/788707_20180129_S -t 16
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

