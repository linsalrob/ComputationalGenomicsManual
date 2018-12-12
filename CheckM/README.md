# Checking Genome Contamination (Redundancy) and Completeness with CheckM

When you sequence a genome, and especially when you construct [metagenome assembled genomes](../CrossAssembly),you would like to know whether the genome is complete (i.e. it has all the genes you would expect to be there), and whether there is data from more than one organism in the sequence.

One approach you can take to explore your genome is to use [GenomePeek](../GenomePeek/), which will show you contamination, but will not show you completeness.

<a href="https://ecogenomics.github.io/CheckM/" target="_blank"><img src="images/checkm.png" width="600 px" alt="CheckM logo" title="The CheckM logo" /></a>

An alternative is to use [CheckM](https://ecogenomics.github.io/CheckM/) that does several computations. 

CheckM uses [prodigal](../ORFCalling/) to identify the genes in your sequences, 

First, it places your genome on a tree, and then it uses that phylogenetic placement to identify sets of genes that should be in your genome.

Next, it looks for those genes using [hmmer](../Databases) to search your genome. The completeness is an estimate of the fraction of genes that are expected to be there which were actually found. The contamination is based on identifying the number of single copy genes, that should only be there once.

There are lots of different workflows that you can complete using CheckM and  you should [check out their manual](https://github.com/Ecogenomics/CheckM/wiki) for more commands and workflows.

In the [previous steps](../CrossAssembly/) we made a bin based on contigs that are correlated with each other across multiple genomes. We can put those contigs into a directory, and use CheckM to test what is our coverage or completeness.

The simplest workflow that we use is:

```bash
checkm lineage_wf GenomeBins/ CheckMOut
```

*Note:* If you have more than one thread on your computer you can make this run a lot faster by using them. For example,

```bash
checkm lineage_wf -t 16 GenomeBins/ CheckMOut
```
will run `checkm` with 16 threads.

`GenomeBins` is the name of a directory that has the genome bins that you want to analyze.

You can find a detailed description of the steps that CheckM takes [in the lineage_wf](https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow) on the [CheckM wiki](https://github.com/Ecogenomics/CheckM/wiki). 

*Note;* In the  [previous steps](../CrossAssembly/) we filtered contigs based on their Pearson correlation coefficient. You can make repeated attempts at binning sequences with different correlation coefficients or other parameters, and then testing them with CheckM to see the result.


