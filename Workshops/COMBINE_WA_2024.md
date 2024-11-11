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
ln -s /scratch/courses01/$USER/miniforge3 software/miniforge3
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

We can use [countfasta.py](https://raw.githubusercontent.com/linsalrob/ComputationalGenomicsManual/refs/heads/master/Python/countfasta.py) to count the lengths of the sequences:

```
python countfasta.py -f megahit_assembled/788707_20180129_S/final.contigs.fa
```


## Making Sankey Plots

We can use the [countfasta.py](https://raw.githubusercontent.com/linsalrob/ComputationalGenomicsManual/refs/heads/master/Python/countfasta.py) and the related [countfastq.py](https://raw.githubusercontent.com/linsalrob/ComputationalGenomicsManual/refs/heads/master/Python/countfastq.py) scripts to make beautiful SanKey plots. 

Open a text file, and by counting the sequences and the contigs, create some data that looks like:

```
fastq [1176881906] fastp
fastq [1425302] low quality
fastp [672903626] human
fastp [515964090] not human
not human [429542979] sequence similarity
not human [86421111] unknown
sequence similarity [7100783] Eukaryote
sequence similarity [330717020] Bacteria
sequence similarity [71361] Archaea
sequence similarity [1029932] Virus
sequence similarity [89215792] Multiclass
```

Then, we use the awesome [sankeymatic](https://sankeymatic.com/build/) to make our SanKey plots!


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


## BEFORE WE GO ON!

We are running out of space on `/home`, where we logged in, so now we are going to work on `/scratch`.

Here's the hack to make life easy for you:

```
ln -s /scratch/courses01/$USER scratch
```

Now you see a directory called `scratch`. Moving forwards, we are going to do all our work in there.

```
cd scratch
```

(Note the difference between `/scratch` and `scratch`)


## Cross-assembly

We are going to run a cross assembly on this data to get more contigs. I've staged the data on `/scratch/courses01/cf_data`.

Here is the code that we need to run this assembly

```
ALLR1=""; ALLR2="";

for R1 in $(find /scratch/courses01/cf_data/ -name \*R1\*); do
        R2=${R1/R1/R2};
        ALLR1+="$R1,";
        ALLR2+="$R2,";
done;

ALLR1=$(echo $ALLR1 | sed -e 's/,$//');
ALLR2=$(echo $ALLR2 | sed -e 's/,$//');

megahit -1 $ALLR1 -2 $ALLR2 -o megahit_assembly -t 16

```

**Note:** 
1. Let's walk through this code and see what it does!
2. You need some slurm "directives" before you can start that code running.


## Map the reads to all the contigs

We are going to use `minimap`, like we did beore. However, here is a little bit of code that can run `minimap` on all of the samples!

```
READDIR=/scratch/courses01/cf_data/

mkdir -p bam_contigs
for R1 in $(find $READDIR -name \*R1\* -printf "%f\n"); do
	R2=${R1/R1/R2};
	BAM=${R1/_R1.fastq.gz/.contigs.bam};
	minimap2 --split-prefix=tmp$$ -t 8 -a -xsr  megahit_assembled/cross_assembly/final.contigs.fa $READDIR/$R1 $READDIR/$R2 | samtools view -bh | samtools sort -o bam_contigs/$BAM;
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
mkdir bam_contigs_tsv
for BAM in $(find bam_contigs -type f -name \*bam -printf "%f\n"); do 
	OUT=${BAM/.contigs.bam/.tsv}; 
	samtools coverage bam_contigs/$BAM | cut -f 1,7 > bam_contigs_tsv/$OUT; 
done
```

We have created an [example Jupyter notebook](Workshop_MAG_demo.ipynb) so you can see some of the commands

## Identifying contigs that belong to the same genome

We are going to move the data to [Google Colab](https://colab.research.google.com/) to analyse the data and identify contigs that co-occur across multiple samples.


---


# Viral Analysis With Hecatomb


We have already run the hecatomb pipeline for you, and the data is in `/storage/data/hecatomb`, but you need to download the data to your computer before we can analyse the files.

There are three files:


- [bigtable](https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/bigtable.tsv.gz?download=)
- [VMR table of viral families and hosts](https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/VMR_MSL39_v1.ascii.tsv.gz?download=)
- [CF Metadata table](https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/CF_Metadata_Table-2023-03-23.tsv.gz?download=)


## Start a Google Colab Notebook

Connect to [Google Colab](https://colab.google.com/)


### Load the data into pandas dataframes:

```
import os
import sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import statsmodels.stats.api as sms
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
```

### Copy the files into your Jupyter notebooks and then find the files

```
os.chdir('Data')
os.listdir()
```

### Read the data into pandas data frames

```
data = pd.read_csv('bigtable.tsv.gz',compression='gzip',header=0,sep='\t')
# without downloading - copy link address from above
# data = pd.read_csv('https://github.com/linsalrob/CF_Data_Analysis/raw/main/hecatomb/bigtable.tsv.gz?download=',compression='gzip',header=0,sep='\t')
metadata = pd.read_csv('CF_Metadata_Table-2023-03-23.tsv.gz',compression='gzip',header=0,sep='\t')
vmr = pd.read_csv('VMR_MSL39_v1.ascii.tsv.gz', compression='gzip',header=0,sep='\t')
```

Take some time to look at the tables. 

How many columns are there in each?
How many rows?
What are the data types in the rows and columns?

### First, filter the bigtable for _just_ the Viruses

```
viruses = data[(data.kingdom == "Viruses")]
```

### Group them by alignment type, length, and family

```
virusesGroup = viruses.groupby(by=['family','alnType','alnlen','pident'], as_index=False).count()
```

### Now use seaborn to make some pretty plots!

First, we start by creating a style that we like.

```
#styling
sizeScatter = 10 * virusesGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(12,8)})

# and now plot the data
g = sns.FacetGrid(virusesGroup, col="family", col_wrap=10)
g.map_dataframe(sns.scatterplot, "alnlen", "pident", alpha=.1, hue="alnType", sizes=(100,500), size=sizeScatter)
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
```


That's a lot of data!

#### Restrict this plot to just E value < 1 x 10<sup>-20</sup>

```
virusesFiltered = viruses[viruses.evalue<1e-20]
virusesGroup = virusesFiltered.groupby(by=['family','alnType','alnlen','pident'], as_index=False).count()
g = sns.FacetGrid(virusesGroup, col="family", col_wrap=10)
g.map_dataframe(sns.scatterplot, "alnlen", "pident", alpha=.1, hue="alnType", sizes=(100,500), size=sizeScatter)
for ax in g.axes.flat:
    ax.tick_params(axis='both', labelleft=True, labelbottom=True)
    ax.axhline(y=75, c='red', linestyle='dashed', label="_horizontal")
    ax.axvline(x=150, c='red', linestyle='dashed', label="_vertical")
    
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
```


### More detailed analysis

There are two proteins that are bogus. We don't know why, so we just ignore. For the next steps, we also _only_ look at the amino acid hits, not the nucleotide hits.

```
data = pd.read_csv('bigtable.tsv.gz',compression='gzip',header=0,sep='\t')

blacklist = ['A0A097ZRK1', 'G0W2I5']
virusesFiltered = data[(data.alnType == "aa") & (data.kingdom == "Viruses") & (~data.targetID.isin(blacklist) & (data.evalue < 1e-20))]
```

### Patient and Sample Date

Now we create new columns for the patient and the sample date

```
virusesFiltered[['patient', 'date', 'Sputum or BAL']] = virusesFiltered['sampleID'].str.split('_', expand=True)
```

## Merge the real data and the host data

```
virusesFiltHost = pd.merge(virusesFiltered, vmr[['Family', 'Host source']], left_on="family", right_on="Family", how='left')
```

Take a look at this new table. What have we accomplished? What was done? How did it work?

## Correct more errors

There are two families that are incorrectly annotated in this table. Let's just manually correct them:

```
really_bacterial = ['unclassified Caudoviricetes family', 'unclassified Crassvirales family']
virusesFiltHost.loc[(virusesFiltHost['family'].isin(really_bacterial)), 'Host source'] = 'bacteria'
```

## Group-wise data

Now lets look at some data, and put that into a group. We'll also remember what we've looked 

```
to_remove = []
bacterial_viruses = virusesFiltHost[(virusesFiltHost['Host source'] == 'bacteria')]
to_remove.append('bacteria')
```

Can you make a plot of just the bacterial viruses like we did?

### How many reads map per group, and what do we know about them

```
host_source = "archaea"
rds = virusesFiltHost[(virusesFiltHost['Host source'] == host_source)].shape[0]
sps = virusesFiltHost[(virusesFiltHost['Host source'] == host_source)].species.unique()
fams = virusesFiltHost[(virusesFiltHost['Host source'] == host_source)].family.unique()
print(f"There are {rds} reads that map to {host_source} viruses, and they belong to {len(sps)} species ", end="")
if len(sps) < 5:
    spsstr = "; ".join(sps)
    print(f"({spsstr}) ", end="")
print(f"and {len(fams)} families: {fams}")
to_remove.append(host_source)
```

Can you repeat this for each of (`plants` OR `plants (S)`), `algae`, `protists`, `invertebrates`, `fungi`, 

## Filter out things we don't want

Once we remove the above, what's left?

```
virusesFiltHost = virusesFiltHost[~virusesFiltHost['Host source'].isin(to_remove)]
virusesFiltHost['Host source'].unique()
```

## Now look at what we have

```
#filter
virusesGroup = virusesFiltHost.groupby(by=['family','alnlen','pident'], as_index=False).count()

#styling
sizeScatter = 10 * virusesGroup['count']
sns.set_style("darkgrid")
sns.set_palette("colorblind")
sns.set(rc={'figure.figsize':(20,10)})

g = sns.FacetGrid(virusesGroup, col="family", col_wrap=6)
g.map_dataframe(sns.scatterplot, "alnlen", "pident", alpha=.8, sizes=(100,500), size=sizeScatter)
for ax in g.axes.flat:
    ax.tick_params(axis='both', labelleft=True, labelbottom=True)
    ax.axhline(y=80, c='red', linestyle='dashed', label="_horizontal")
    ax.axvline(x=150, c='red', linestyle='dashed', label="_vertical")

g.fig.subplots_adjust(hspace=0.4)
g.set_axis_labels("Alignment length", "Percent Identity")
plt.legend(bbox_to_anchor=(5.0,1), loc=0, borderaxespad=2,ncol=6, shadow=True, labelspacing=1.5, borderpad=1.5)
plt.show()
```

## Now look at each family separately

```
virusesFiltHost[virusesFiltHost.family == 'Retroviridae']
```










