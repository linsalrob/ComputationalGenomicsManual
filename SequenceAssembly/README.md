# Sequence Assembly

The essential problem that sequence assembly is trying to overcome is that the average microbial genome is 2,000,000 bp, while the typical sequence read length is 150-300 bp for Illumina sequences and upto 50,000 bp for PacBio or Nanopore sequences.

With such (relatively) short sequences, how can we assemble a whole, or nearly whole, genome?

The answer is by repetitively sequencing the same thing over and over again! If we start each sequence at a random location, and we have enough sequences, eventually we can join those sequences together to form what we call **contigs**.

There are four types of sequence assembly algorithms:

1. Naive assemblers which just try and find all matching pairs of reads
2. Greedy assemblers which start with one read and keep adding reads until you can not find any more matches, and then start with the next read.
3. Overlap-layout-consensus assemblers which layout the reads looking for overlaps between them. The overlaps are usually refined by a Smith-Watermann search, and then a consensus constructed.
4. de Bruijn graph assemblers

This table describes some of the common sequence assemblers that you will run across. 

  **Name** | **Type** | **Sequencing Tech** | **Citation** | **Documentation** | **Homepage**
  --- | --- | --- | --- | --- | ---
SPAdes | genomes, single-cell, metagenomes, ESTs | Illumina, Solexa, Sanger, 454, Ion Torrent, PacBio, Oxford Nanopore |  [Nurk *et al.* 2013](https://link.springer.com/chapter/10.1007%2F978-3-642-37195-0_13) | [version 3.12 manual](http://cab.spbu.ru/files/release3.12.0/manual.html) |   [SPAdes](http://bioinf.spbau.ru/en/spades)
Velvet |  genomes  | Sanger, 454, Solexa, SOLiD | [Zerbino and Birney, 2008](https://genome.cshlp.org/content/18/5/821.long) | [version 1.12 manual](https://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf) |   [EBI](http://www.ebi.ac.uk/~zerbino/velvet/)
Canu |  genomes  | PacBio/Oxford Nanopore reads |  [Koren *et al.* 2017](https://genome.cshlp.org/content/27/5/722) | [manual for all versions](https://canu.readthedocs.io/en/latest/quick-start.html) | [Git repo](https://github.com/marbl/canu)
MaSuRCA | Any size, haploid/diploid genomes | Illumina and PacBio/Oxford Nanopore data, legacy 454 and Sanger data | [Zimin A, *et al.* 2017](https://www.ncbi.nlm.nih.gov/pubmed/28130360) |   [Git Repo](https://github.com/alekseyzimin/masurca) | [Git Repo](https://github.com/alekseyzimin/masurca)
Hinge | Small microbial genomes | PacBio/Oxford Nanopore reads | [Kamath *et al.* 2017](https://genome.cshlp.org/content/27/5/747.full) | [jupyter notebook](https://github.com/HingeAssembler/HINGE-analyses) |   [Git repo](https://github.com/HingeAssembler/HINGE)
Unicycler | Illumina-only data, and can optimize SPAdes | [Wick *et al.*](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595)  | [Git Repo](https://github.com/rrwick/Unicycler)  
Flye | De novo assembly of long reads (PacBio/Oxford Nanopore) | [Kolmogorov *et al.*](https://doi.org/10.1038/s41592-020-00971-x) | [Git Repo](https://github.com/fenderglass/Flye)
miniasm + minipolish | Long read assembler and polishing together  | [Wick and Holt](https://f1000research.com/articles/8-2138) | [Git repo](https://github.com/rrwick/Minipolish)
raven | Assembler for long, uncorrected reads | TBD | [Git repo](https://github.com/lbcb-sci/raven)
Trycycler | Not really an assembler, *per se*, but more an approach to merging assemblies. | [DOI:10.5281/zenodo.3965017](https://doi.org/10.5281/zenodo.3965017) | [git Repo](https://github.com/rrwick/Trycycler)

For the most comprehensive comparison of sequence assemblers, we encourage you to review [Wick RR, Holt KE. Benchmarking of long-read assemblers for prokaryote whole genome sequencing. F1000Research. 2019;8(2138).](https://f1000research.com/articles/8-2138)

We use the [St. Petersburg genome assembler, SPAdes](http://cab.spbu.ru/software/spades/) and the version installed on the AWS instances is 3.12.0 for which the [manual is here](http://cab.spbu.ru/files/release3.12.0/manual.html)

For Nanopore reads we typically use the CANU assembler.

## Running SPAdes

SPAdes is easy to run! The basic command is 

```
spades.py
```

The program takes a couple of inputs - your `fastq` files, for example that you download from [../Databases/SRA](../Databases/SRA).

If you have paired end reads, you need to add `-1` for the left pairs (the file called xxx\_1.fastq) and `-2` for the right pairs (the file called xxx\_2.fastq). Note that spades handles `gzip` compressed files, and you do not need to decompress them!

If you unpaired reads, you can specify that with the `-s` flag.

You also need to provide an output directory name where the results will be written using the `-o` flag.

Your final command might look something like:

```
spades.py -1 fastq/ERS011900_pass_1.fastq.gz -2 fastq/ERS011900_pass_2.fastq.gz -o assembly
```

## SPAdes output files

SPAdes makes a lot of files and directories in the output, and this summarizes what those files are. Of course, more details can be found in the [SPAdes manual](http://cab.spbu.ru/files/release3.12.0/manual.html)

* **`scaffolds.fasta`** contains the scaffolds generated by SPAdes and is the **output file you want to use**.
* the directory `/corrected/` contains reads corrected by BayesHammer in compressed fastq format
* `contigs.fasta` contains the contigs before they are scaffolded into scaffolds. Often this is similar to the scaffolds.fasta depending on how much scaffolding information there is
* `assembly_graph.gfa` contains the assembly graph and scaffolds paths in GFA 1.0 format
* `assembly_graph.fastg` contains the assembly graph in FASTG format
* `contigs.paths` contains paths in the assembly graph corresponding to contigs.fasta. This is how the graph is resolved into contigs.
* `scaffolds.paths` contains paths in the assembly graph corresponding to scaffolds.fasta.
* `K21`, `K33`, `K55`, etc are directories containing the de Bruijn graph assemblies for different lengths of *k*
* `before_rr.fasta` are the assembled contigs before repeat resolution has been applied.
* `dataset.info` and `input_dataset.yaml` contain information about the sequence read files that were supplied.
* `params.txt` is a summary of all the spades parameters
* `spades.log` is the log that was printed to the screen while SPAdes was running. This contains lots of information about the assembly process.


