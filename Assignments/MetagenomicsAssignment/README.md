# The metagenome assignment

For this assignment we are going to complete some metagenomics analyses. I have provided several [metagenomics datasets](../Datasets/). Please do **not** use the *Drinking Water* data set as it is a 16S sequencing data set and will not work for these analyses.  Also, as noted below, if you use the **Algae** data set, you will get the minimum marks possible as that is the example that we've worked through and you can just copy and paste the commands without thinking.

*Extra credit:* If you want to use your own data set, you are welcome to do so.

*Extra credit:* If you want to find another data set to use, you can [search the SRA](../Databases/SRA.html) for a metagenomic data set and use that instead. You should choose a metagenome that has at least three runs associated with it, as later in the assignment we will use those to create metagenome assembled genomes.

## Part 1. Describe your metagenomes

* Where did they come from?
* Who did the sequencing?
* What sequencing approach did they use?
* How much data did they generate?
* Why did they sequence it?
* Was there a question they were trying to answer?

* *Extra credit:* Using a different data set that is not one of those provided here will earn you extra credit! 
* *Miniumum credit!* If you use the Algae data set, you will definitely not get extra credit and will probably get the bare minimum credit since all you need to do is copy and paste the commands. Use a different data set!

## Part 2. Annotating the organisms present in the metagenome

First, we are going to identify the organisms present in the metagenome. There are several ways to do that, but I recommend [focus](https://edwards.sdsu.edu/FOCUS) as it is installed in the AWS instances.

* Generate a table that describes the genus/species of the organisms that are present in the samples.
* Are the same organisms the most abundant in each of the samples that you analyzed?
* Do you expect those organisms to be present based on what you know about the sample? For example, if your sample is from the marine environment, are those organisms typically found in the marine environment? If your sample is from the human gut, are those common gut bacteria?

There are several other ways to analyse the metagenomes, including [mg-rast](https://www.mg-rast.org/), [MGnify](https://www.ebi.ac.uk/metagenomics/), [CLARK](http://clark.cs.ucr.edu/), [MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan), [GenomePeek](https://edwards.sdsu.edu/GenomePeek), and plenty of others. In fact, you can [read about a host](https://www.nature.com/articles/nmeth.4458) of different software for analysing metagenomes in the CAMI paper: Sczyrba A, Hofmann P, Belmann P, Koslicki D. 2017. [Critical assessment of metagenome interpretation—a benchmark of metagenomics software](https://www.nature.com/articles/nmeth.4458). Nature.

* *Extra credit:* If you annotate the organisms present in the sample using another method: 
	a) what method did you use and how did you use it (e.g. was it on the AWS instance, a web site, somewhere else)
	b) did it give you similar results to FOCUS. What were the critical differences?
	c) Which method would you use again, and why?

## Part 3. Annotating the functions present in the metagenome

Just as with annotating the organisms present in the metagenome, there are several different methods to annotate the functions present in the metagenome.

One approach is to use [real time metagenomics](../RTMg/), either in the web version or the stand alone version (*pro tip*: the web version is limited in how many queries it can make at once. The standalone version is not limited. If there is a class running, you probably want to use the standalone version!)

You can also use [super-focus](https://edwards.sdsu.edu/SUPERFOCUS/). It is installed on the AWS instance, though before you start you will need to use this command to download the appropriate databases.

```bash
superfocus_downloadDB -a diamond
```

You can also use [mg-rast](https://www.mg-rast.org/), [MGnify](https://www.ebi.ac.uk/metagenomics/),  and of course other software described in the CAMI paper: Sczyrba A, Hofmann P, Belmann P, Koslicki D. 2017. [Critical assessment of metagenome interpretation—a benchmark of metagenomics software](https://www.nature.com/articles/nmeth.4458). Nature.

* How did you annotate the functions present in the metagenomic samples?
* How long did it take?
* What are the most abundant functions in each sample?
* Do those functions make sense given the environment that your sample comes from? e.g. if it is a marine sample are those functions you think would be important in the marine environment? If it is a gut sample, are those functions you think would be important in the gut?

* *Extra credit:* If you annotate the functions present in the sample using another method: 
a) what method did you use and how did you use it (e.g. was it on the AWS instance, a web site, somewhere else)
b) did it give you similar results to the first method you chose?. What were the critical differences?
c) Which method would you use again, and why?

*Extra credit:* How do you summarize the functions? We always see stacked bar charts (gross) or pie charts (grosser) for functions present in metagenomes. I always give additional credit for novel, unique, and different data visualization and presentation approaches (although this is not a data visualization class!)


## Part 4. Metagenome Assembled Genomes

As noted above the metagenome needs to have more than one run associated with it for this step to work. Ideally you will have a metagenome with 4 or more runs.

For this aim, we're going to use [cross-assembly](../CrossAssembly) to analyze the metagenomes and try to identify complete genomes present in the data. 

As you work through the steps associated with cross assembly, here are some questions to answer:

* What are the assembly statistics? For example:
	* How many contigs did you get from the assembly?
	* What is the N<sub>50</sub> and N<sub>75</sub>?
* How many reads mapped to the contigs?
* Why did the other reads not map to those contigs since they were all used in the assembly?
* Plot a graph of the abundance of reads mapping to each of the contigs 
* After you have run `crAss`, what are the most correlated contigs?
* *Extra credit:* If you pull out the reads associated with contigs that have Pearson correlation (r) > 0.95 and repeat the assembly of just those reads, do you get better assembly statistics?

## Part 5. Checking the metagenome annotated genomes with CheckM

Once we have assembled the genomes, we want to check the completeness and contamination of the bins. Take a set of highly correlated contigs and create a directory with them. They are a metagenome bin.

Next, we're going to run [CheckM](https://ecogenomics.github.io/CheckM/) on those contigs, using the [description here](../../CheckM).

* What is the completeness of your genome bin?
* What is the contamination of your genome bin?
* Do you think this represents a whole genome?

*Extra credit:* "Highly correlated contigs" is not a well defined term! You probably used Pearson correlation > 0.95 since that was in the previous question! What happens to completeness and contamination as you decrease the correlation coefficient of the contigs in your bin?

## Part 6. Visualization with anvi'o.

Finally, we're going to work through the anvi'o workflow to create  a visualization of our bins.

The workflow is [described here](../../ANVIO) and is also [thoroughly described on the anvi'o website](http://merenlab.org/2016/06/22/anvio-tutorial-v2/).

* What does the anvi'o visualization tell you about your metagenome bins?
* Paste a copy of the image into your report.

*Extra credit:* Can you use anvi'o to refine the quality of your bins? If you do that, what happens to the checkM scores?

## Part 7. Discussion
In the beginning you described the metagenome and why it was sequenced.

* What did you find?
* How has your analysis contributed to our understanding of this sample?
* Did you analysis recapitulate the results of the original authors?
* Did you find something they missed?
* Did they find something you missed?
