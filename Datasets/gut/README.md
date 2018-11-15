# Gut data

This data is from [SRA project SRP074153](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP074153)

This is a random community data set

This data comes from

Brooks B, Olm MR, Firek BA, Baker R, Thomas BC, Morowitz MJ, Banfield JF. 2017. Strain-resolved analysis of hospital rooms and infants reveals overlap between the human and room microbiome. [Nat Commun 8:1814](https://www.nature.com/articles/s41467-017-02018-w)

## Title

Metagenomes from 11 human infant fecal samples hospitalized in the same intensive care unit

## Abstract

Bacteria that persist in hospitals can contribute to the establishment of the microbiome in newborns and the spread of hospital-acquired diseases. Yet we know little about microbial communities in hospitals, or about the extent to which persistent vs. recently immigrated bacterial strains establish in the gastrointestinal tracts of hospitalized individuals.In combination with BioProject PRJNA273761 (10 infants / 55 samples) we analyzed strain-resolved genomes obtained from a total of 202 samples collected over a three-year period from 21 infants hospitalized in the same intensive care unit.Strains were rarely shared, consistent with prior analysis of a subset of these data. Enterococcus faecalis and Staphylococcus epidermidis, common gut colonists, exhibit diversity comparable to that of NCBI reference strains, suggesting no recent common ancestor for all populations in this hospital setting. Thus, we infer multiple introduction events for these species. Despite the rarity of shared strains, strains of five species exhibiting a degree of sequence variation consistent with in situ diversification were identified in different infants hospitalized three years apart. Three were also detected in multiple infants in the same year, suggesting that these strains are unusually widely dispersed and persistent in the hospital environment. Persistent strains were not significantly different from non-persistent strains with regards to pathogenicity potential including antibiotic resistance. Notably, non-identical siblings had multiple abundant strains in common, even 30 days after birth and antibiotic administration, suggesting overlapping strain sources and/or genetic selection. Our approach can be used in order to study microbial dynamics in hospitals and provides an important step towards directing health-promoting colonization in hospitalized individuals.

## Design

DNA was extracted using the MO BIO PowerSoil DNA Isolation kit; libraries were made Illuminas Nextera kit with average insert sizes of 500/900 bp

## Runs

SRR3466404 SRR3506419 SRR3506420 SRR3546776 SRR3546778 SRR3546779 SRR3546780 SRR3546781 SRR3546782



We downloaded the data from SRA like this:

```
for SRR_ID in SRR3466404 SRR3506419 SRR3506420 SRR3546776 SRR3546778 SRR3546779 SRR3546780 SRR3546781 SRR3546782; do
	fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip -N 5000 -X 255000 $SRR_ID
done
```

## Data

* [SRR3466404 pass 1](fastq/SRR3466404_pass_1.fastq.gz)
* [SRR3466404 pass 2](fastq/SRR3466404_pass_2.fastq.gz)

* [SRR3506419 pass 1](fastq/SRR3506419_pass_1.fastq.gz)
* [SRR3506419 pass 2](fastq/SRR3506419_pass_2.fastq.gz)

* [SRR3506420 pass 1](fastq/SRR3506420_pass_1.fastq.gz)
* [SRR3506420 pass 2](fastq/SRR3506420_pass_2.fastq.gz)

* [SRR3546776 pass 1](fastq/SRR3546776_pass_1.fastq.gz)
* [SRR3546776 pass 2](fastq/SRR3546776_pass_2.fastq.gz)

* [SRR3546778 pass 1](fastq/SRR3546778_pass_1.fastq.gz)
* [SRR3546778 pass 2](fastq/SRR3546778_pass_2.fastq.gz)

* [SRR3546779 pass 1](fastq/SRR3546779_pass_1.fastq.gz)
* [SRR3546779 pass 2](fastq/SRR3546779_pass_2.fastq.gz)

* [SRR3546780 pass 1](fastq/SRR3546780_pass_1.fastq.gz)
* [SRR3546780 pass 2](fastq/SRR3546780_pass_2.fastq.gz)

* [SRR3546781 pass 1](fastq/SRR3546781_pass_1.fastq.gz)
* [SRR3546781 pass 2](fastq/SRR3546781_pass_2.fastq.gz)

* [SRR3546782 pass 1](fastq/SRR3546782_pass_1.fastq.gz)
* [SRR3546782 pass 2](fastq/SRR3546782_pass_2.fastq.gz)

* [SRR3466404 pass 1](SRR3466404_pass_1.fastq.gz)
* [SRR3466404 pass 2](SRR3466404_pass_2.fastq.gz)

* [SRR3506419 pass 1](SRR3506419_pass_1.fastq.gz)
* [SRR3506419 pass 2](SRR3506419_pass_2.fastq.gz)

* [SRR3506420 pass 1](SRR3506420_pass_1.fastq.gz)
* [SRR3506420 pass 2](SRR3506420_pass_2.fastq.gz)

* [SRR3546776 pass 1](SRR3546776_pass_1.fastq.gz)
* [SRR3546776 pass 2](SRR3546776_pass_2.fastq.gz)

* [SRR3546778 pass 1](SRR3546778_pass_1.fastq.gz)
* [SRR3546778 pass 2](SRR3546778_pass_2.fastq.gz)

* [SRR3546779 pass 1](SRR3546779_pass_1.fastq.gz)
* [SRR3546779 pass 2](SRR3546779_pass_2.fastq.gz)

* [SRR3546780 pass 1](SRR3546780_pass_1.fastq.gz)
* [SRR3546780 pass 2](SRR3546780_pass_2.fastq.gz)

* [SRR3546781 pass 1](SRR3546781_pass_1.fastq.gz)
* [SRR3546781 pass 2](SRR3546781_pass_2.fastq.gz)

* [SRR3546782 pass 1](SRR3546782_pass_1.fastq.gz)
* [SRR3546782 pass 2](SRR3546782_pass_2.fastq.gz)
