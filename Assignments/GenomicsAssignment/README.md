# Whole genome analysis assignments

This marks the beginning of several assignments where data will flow from one assignment to the other. We are going to study the genome of *Klebsiella pneumonia*. This is one of the so-called [ESKAPE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4871955/) pathogens that are a serious threat to public health because of the rise of multiple drug resistance. 

The ESKAPE pathogens are:

+ *Enterococcus faecium*
+ *Staphylococcus aureus*
+ *Klebsiella pneumoniae*
+ *Acinetobacter baumannii*
+ *Pseudomonas aeruginosa*
+ *Enterobacter* species


In the next series of assignments, we will recapitulate some of the amazing work performed by Professor Kat Holt and colleagues who did a global analysis of the emergence of *Klebsiella* as a pathogen. 

This is a particularly nasty pathogen that causes problems worldwide. It has virulence genes (that enable it to cause disease) and antibiotic ﻿﻿﻿﻿﻿﻿﻿﻿﻿resistance genes (that make it harder to treat!) and they hoped to learn many of the differences in these genomes

You can read more about their analysis in their open access paper published in [Proceedings of the National Academy of Sciences](http://www.pnas.org/content/early/2015/06/17/1501049112).

I particularly encourage you to read [this terrific blog post](https://holtlab.net/2015/06/23/population-genomics-of-klebsiella/) by Prof. Holt where she describes some of the key points of their work. During these assignments, we're going to try and recapitulate some of their findings!

Finally, the [metadata associated with the samples is available at microreact](https://microreact.org/project/VJdoJhfkx). You should take a look at that, as it will get you started with understanding the differences in the data.

## Part 1. Download and assemble the data

You can find a complete list of the sequence accessions for these samples [here](klebsiella.txt). These IDs are from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena) but you can also download them from the [NCBI sequence read archive](https://www.ncbi.nlm.nih.gov/sra/) which is what we will do here. (By the way, you will notice that we have trimmed out some of the sequences described in the [metadata](https://microreact.org/project/VJdoJhfkx) because they are already assembled! if you want to proceed with some of those, you can skip the assembly step.

### Assembly note
 
Sequence assembly is, by its very nature, extremely memory intensive. For complete assembly of the data using AWS, I recommend using a t2.medium machine that has two cores and 4 GB RAM or a t2.large machine that has four cores and 8 GB RAM. (When working with complex data like metagenomes that we will talk about later, even this is not enough and we step up to machines with hundreds of GB of RAM or even many TB of RAM!)

If you have a running AWS instance, you can change the state of the machine:

1. Log in to the AWS Console
2. Stop the machine by choosing the machine and clicking Stop from the Actions menu (under instance state). This will suspend but not delete it
3. Once it has stopped, you can change the type by choosing Actions --> Instance Settings --> Change Instance Type
4. Restart the instance as the new type.

### Getting started with the assignment

First, you need to download the data from the SRA using [fastq-dump](../Databases/SRA#fastq-dump). 

This should create one, two or three files for your data. Note that fastq-dump also leaves a copy of the data in the `~/ncbi/` directory that you can go ahead and delete to save space.

Next, we need to assemble that data. For this assignment we are going to use spades.

spades.py has a lot of options - run `spades.py` without any arguments to see what they are (pro tip: you may want to pipe the output of that to `less`). The main ones we will use are `-1` and `-2` for left and right paired end reads, and `-s` for unpaired reads. If `fastq-dump` gave you one file, use `-s` with that file name. If `fastq-dump` gave you two files, one will be used with the `-1` option and the other with the `-2` option, and finally if `fatq-dump` gave you three files, you will use the paired files with `-1` and `-2` and the unpaired file with `-s`. Note that with spades you can specify many multiples of paired end files, and many additional unpaired files.

Once you have assembled that data can you generate the data that describes:

* Unassembled Size (bp)
* Unassembled size (# reads)
* Number of Contigs
* Longest contig
* N<sub>50</sub>

## Part 2. Identify the open reading frames

Once you have assembled the genome into contigs, we are going to annotate the [open reading frames](../ORFCalling/) in the genome.

Describe:

* The number of predicted genes
* The length of the longest gene (bp)

## Part 3. Identify the RNA genes.

The next step is to identify the [RNA genes](../tRNA_rRNA) in the genome. We will identify both the rRNA and the tRNA genes.

Describe:

* The number of rRNA
* Number of tRNA


## Part 4. Identify the functions of the proteins.

Function of the longest protein	Number of hypothetical proteins	Number of proteins with known functions



