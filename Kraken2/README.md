# Kraken2

Kraken2 uses _k_-mers to identify the taxonomy of the microbes in your sample. In essence, they have taken all complete genomes, and then identified all _k_-mers that are unique to each taxonomic level. Through some nifty computing, and special [data structures](https://www.youtube.com/watch?v=zgCnMvvw6Oo&list=PLpPXw4zFa0uKKhaSz87IowJnOTzh9tiBk), they have figured out how to search this very efficiently.

There are a wide range of pre-built [kraken databases](https://benlangmead.github.io/aws-indexes/k2) that you can download, so you do not need to go to the effort of building them yourself.


When [installing Kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual#installation), I recommend setting the `KRAKEN2_DB_PATH` and `KRAKEN2_DEFAULT_DB` variables, and then you do not need to specify them on the command line.


To run Kraken2, use this incantation:

```bash
kraken2 --paired --threads 4 --report kraken_taxonomy.txt --output kraken_output.txt \
	fastq/reads_1.fastq fastq/reads_2.fastq
```

This will output two files:

* `$SRR.kraken_output.txt` contains the [standard kraken output](https://github.com/DerrickWood/kraken2/wiki/Manual#output-formats):
    - A code (_C_ or _U_) indicating whether the read was classified or not
    - The read ID from the fastq file
    - The taxonomy ID assigned to the read if it is classified, or 0 if it is not classified
    - The length of the sequence in base pairs. Because we are using paired end reads, there are two lengths (R1\|R2)
    - A space-separated list of the lowest common ancestor for each sequence that indicates how many kmers map to which taxonomic IDs. Because we have paired end information, there is a `|:|` separator between the R1 and R2 information
* `$SRR.kraken_taxonomy.txt` contains the [standard kraken report](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format):
    - Percent of fragments at that taxonomic level
    - Number of fragments at that taxonomic level (the sum of fragments at this level and all those below this level)
    - Number of fragments exactly at that taxonomic level
    - A taxonomic level code:  `U`nclassified, `R`oot, `D`omain, `K`ingdom, `P`hylum, `C`lass, `O`rder, `F`amily, `G`enus, or `S`pecies. If the taxonomy is not one of these the number indicates the levels between this node and the appropriate node. See [the docs](https://github.com/DerrickWood/kraken2/wiki/Manual#sample-report-output-format) for more information.
    - NCBI Taxonomic name
    - Scientific name


For more information about Kraken2, [see the wiki page](https://github.com/DerrickWood/kraken2/wiki/Manual)


If you are using the HPC at Flinders University, the details on [this page](https://fame.flinders.edu.au) will show you how to install and use Kraken2


