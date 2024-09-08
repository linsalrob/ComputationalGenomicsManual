# Installing software using CONDA

We are going to install all our software using conda and bioconda. First, we need to install conda:

## Download conda

Head to the [miniforge download page](https://github.com/conda-forge/miniforge) and download the appropriate installer.

You want :

 - **Linux Installer** 
 - **x86_64**
 - **Python** with the biggest number
 - **Miniconda3 Linux 64-bit**
 
Here is a quick way to download that script!

Right click on the appropriate link, and choose *Copy link address*. Go back to your terminal window (probably Putty) and type `wget` and a space, and then paste the URL that you just copied. Press return and it should download the file for you!

The file should be called `Miniconda3-latest-Linux-x86_64.sh` but in case it is not, just substitute the appropriate file name below. Remember that you can check with `ls -ltr` to see the newest file downloaded.

Run the miniconda installer:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

This will ask you some questions, and you can pretty much accept the default answer to all the questions.

Once the installer has finished the best way to continue is to log out, and then log back in. This will reset your account and you will have conda activated. At the bottom left of your screen you should see it say `(base)` which means that you are in the base conda installation.


## Install your first bioinformatics package

We are going to install `fastp` and test to see if it works. This will demonstrate how to install a conda package.

```bash
mamba install -c conda-forge -c bioconda fastp 
```

Again, this will work on resolution of the packages for you and ask if you are sure. Once it is complete, you should be able to issue the command:

```bash
fastp -v
```

to see the version of `fastp` that has been installed.

You can also run 

```bash
fastp -h
```

To get the full `fastp` help menu


## Install snakemake

For the next steps of this tutorial, we are going to use [snakemake](https://snakemake.readthedocs.io/en/stable/) to run some software. So we are going to use `mamba` to install that:

```bash
mamba install -c conda-forge -c bioconda snakemake
```

Once that has completed, `snakemake -v` should show you the current version.

## Install other bioinformatics packages

You can install pretty much any bioinformatics package using conda. The [anaconda website has a complete list](https://anaconda.org/bioconda/repo) and you can visit the [bioconda](https://bioconda.github.io/) page for more information about bioconda.

but to get started, you might want to install:

```bash
mamba install -c conda-forge -c bioconda fastp minimap2 samtools
```

## Adding channels

It is a pain to keep typing `-c conda-forge -c bioconda` so we can just add those two channels to our configuration

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
```

## Updating a package

With `conda` (or `mamba`) you can easily update a package if there is a newer version. For example, to update `snakemake` you would use:

```bash
mamba update snakemake
```

What this means is the _you_ are responsible for ensuring your software is up-to-date. Or not. If you are working on a set of data analyses you may want to keep all the software at the same version so that each time you do an analysis you get comparable answers. With `conda`, you have control over the update cycles, but don't forget from time-to-time you might want to update the software! 

## Environments

If you want to keep different versions of software or run different pipelines you can do that with conda, in what are called `environments`. Each one can have different software. `conda` is clever, because if you have the same software in two different environments you don't need an entire copy of the software. At this stage, you don't need to worry about that, and you can just install everything in the `base` environment. But if you start to run into installation issues, then remember you can separate things into different environments.

