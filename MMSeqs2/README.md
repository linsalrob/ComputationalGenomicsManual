# MMSeqs2

[MMSeqs2](https://github.com/soedinglab/MMseqs2) is a fast sequence searching algorithm that we use in replace of [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) or [diamond](https://github.com/bbuchfink/diamond)

MMSeqs2 has a lot of different options, and we have not yet included them all here, but you can find full details about MMSeqs2 in their [detailed manual](https://mmseqs.com/latest/userguide.pdf).

# MMSeqs2 databases

Like many tools, MMSeqs2 has precomputed databases that you can download. 

There is a [complete list on their website](https://github.com/soedinglab/MMseqs2/wiki#downloading-databases)

You can download a database with the `databases` command. For example, to download the [UniRef50](https://www.uniprot.org/help/uniref) database:

```bash
mkdir -p UniRef50
mmseqs databases --threads 8 UniRef50 UniRef50/UniRef50 /tmp
```

Some of the databases have taxonomy included with them, and that enables you to use `mmseqs easy-taxonomy` to explore the metagenome.

# Easy Taxonomy

We use the MMSeqs2 easy taxonomy a _lot_ for analysing metagenomes, especially by comparing to the [UniRef50](https://www.uniprot.org/help/uniref) database.

First, `mmseqs easy-taxonomy` _requires_ `fasta` files and does not work with `fastq` files. We have a [fast way to convert fastq to fasta](https://edwards.flinders.edu.au/fastq-to-fasta/) or you can find some tools online.

We also take advantage of `mmseqs` [sensitivity sweep](https://github.com/soedinglab/MMseqs2/wiki#set-sensitivity--s-parameter) but you should consider comparing [sensitivity and resources](https://github.com/soedinglab/MMseqs2/wiki#optimizing-sensitivity-and-consumption-of-resources). There is a lot of discussion on the [MMSeqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki) about setting sensitivity.


We typically use this command to run the easy taxonomy:


```bash
mkdir easy-taxonomy
mmseqs easy-taxonomy sequence.fasta UniRef50/UniRef50 easy-taxonomy/sequence_taxonomy /tmp --start-sens 1 --sens-steps 3 -s 7 --threads 32
```

The results will be in a series of files in the `easy-taxonomy` directory, whose names start with `sequence_taxonomy`:

SAGCFN_22_00809_S34_lca.tsv.gz  SAGCFN_22_00809_S34_report.gz  SAGCFN_22_00809_S34_tophit_aln.gz  SAGCFN_22_00809_S34_tophit_report.gz

- `sequence_taxonomy_lca.tsv.gz`: The lowest common ancestor of the sequences in tab separated text.

Example output:

```
R100400180029:20220829140225:V350082744:2:1145432:5:58/2/2      310915  species Pangasianodon hypophthalmus     2       2       1       0.540
```

Columns are:
1. the sequencing read
2. the taxonomy ID from [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree). For example, this is [310915](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/310915/)
3. the taxonomic clade. `Species` in this example
4. The organism name. `Pangasianodon hypophthalmus`
5. 


- `sequence_taxonomy_report.gz` a Kraken2 style output report

Example output:

```
0.8561  9653    9653    species 310915                                                          Pangasianodon hypophthalmus
```

- `sequence_taxonomy_tophit_aln.gz` the `blast m8` format 

Example output:

```
R100400180029:20220829140225:V350082744:2:1145432:5:58/2/2      UniRef50_UPI00147C5152  0.382   163     30      0       0       50      0       163     1.796E-26       108
```

The columns are:

1. Sequence Read
2. Match database ID. In this case from the [UniRef50](https://www.uniprot.org/) we have sequence [UniRef50_UPI00147C5152](https://www.uniprot.org/uniref/UniRef50_UPI00147C5152)
3. Similarity (38.2% identity)
4. Alignment length (163 bases) 
5. Gaps (30 bases)
6. Mismatches (0 bases)
7. Start on the sequence read (0)
8. End on the sequence read (50)
9. Start on the database sequence (0)
10. End on the database sequence (163)
11. E-value (1.796E-26)
12. Bit score


- `sequence_taxonomy_tophit_report.gz` the taxonomy and matches to all of the proteins

Example output

```
UniRef50_UPI00147C5152  6970    0.312   1849.374        0.367   310915  species Pangasianodon hypophthalmus
```

The columns are:

1. Database ID. In this case from the [UniRef50](https://www.uniprot.org/) we have sequence [UniRef50_UPI00147C5152](https://www.uniprot.org/uniref/UniRef50_UPI00147C5152)
2. Number of sequences aligning to target
3. Unique coverage of target uniqueAlignedResidues / targetLength
4. Target coverage alignedResidues / targetLength
5. Average sequence identity
6. [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/tree) identifier
7. Taxonomic level `species`
8. Taxonomic name, in this case `Pangasianodon hypophthalmus`






