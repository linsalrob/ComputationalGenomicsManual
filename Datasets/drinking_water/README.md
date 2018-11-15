# Drinking Water

A drinking water study from U. Adelaide that you can find on the [SRA with accession SRP059994](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP059994)

This is 16S amplicon dataset.

Publication:
Shaw JLA, Monis P, Weyrich LS, Sawade E, Drikas M, Cooper AJ. 2015. Using Amplicon Sequencing To Characterize and Monitor Bacterial Diversity in Drinking Water Distribution Systems. [Appl Environ Microbiol 81:6463–6473](http://aem.asm.org/content/81/18/6463.long)


## Abstract

Drinking water assessments use a variety of microbial, physical, and chemical indicators to evaluate water treatment efficiency and product water quality. However, these indicators do not allow the complex biological communities, which can adversely impact the performance of drinking water distribution systems (DWDSs), to be characterized. Entire bacterial communities can be studied quickly and inexpensively using targeted metagenomic amplicon sequencing. Here, amplicon sequencing of the 16S rRNA gene region was performed alongside traditional water quality measures to assess the health, quality, and efficiency of two distinct, full-scale DWDSs: (i) a linear DWDS supplied with unfiltered water subjected to basic disinfection before distribution and (ii) a complex, branching DWDS treated by a four-stage water treatment plant (WTP) prior to disinfection and distribution. In both DWDSs bacterial communities differed significantly after disinfection, demonstrating the effectiveness of both treatment regimes. However, bacterial repopulation occurred further along in the DWDSs, and some end-user samples were more similar to the source water than to the postdisinfection water. Three sample locations appeared to be nitrified, displaying elevated nitrate levels and decreased ammonia levels, and nitrifying bacterial species, such as Nitrospira, were detected. Burkholderiales were abundant in samples containing large amounts of monochloramine, indicating resistance to disinfection. Genera known to contain pathogenic and fecal-associated species were also identified in several locations. From this study, we conclude that metagenomic amplicon sequencing is an informative method to support current compliance-based methods and can be used to reveal bacterial community interactions with the chemical and physical properties of DWDSs

Project | Description
--- | ---
SRR2080423 | SA 3 tank
SRR2080425 | 1.5kmPostDis
SRR2080427 | SASourceWater2
SRR2080434 | SA 2 CT
SRR2080436 | WA1 outlet


We extracted this data with

```
for SRR_ID in SRR2080423 SRR2080425 SRR2080427 SRR2080434 SRR2080436
	fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip $SRR_ID
done
```

For some reason this data seems mixed up, as there are sequences with each tag in each file. 

We resplit the data with these tags.

The original number os counts per tag are:

Counts | Tag
--- | ---
303,446 | TCGCAGG
276,080 | GCTCGAA
312,009 | CTCGATG
293,741 | ACCAACT
269,702 | GGATCAA

## Sequences

* [barcodes](fastq/barcodes.fastq.gz)
* [sequences](fastq/sequences.fastq.gz)

## Metadata

* [Tab separated format](metadata.tsv)
* [Calypso Format](metadata_calypso.tsv)
