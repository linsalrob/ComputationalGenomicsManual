# Ground Water

This data is from [SRA project SRP075429](https://www.ncbi.nlm.nih.gov/sra/?term=SRP075429)

This data comes from: Hernsdorf AW, Amano Y, Miyakawa K, Ise K, Suzuki Y, Anantharaman K, Probst A, Burstein D, Thomas BC, Banfield JF. 2017. Potential for microbial H2 and metal transformations associated with novel bacteria and archaea in deep terrestrial subsurface sediments. [ISME J 11:1915–1929](https://www.nature.com/articles/ismej201739)

### Title

Potential for microbial H2 and metal transformations associated with novel bacteria and archaea in deep terrestrial subsurface sediments

### Abstract

Geological sequestration in deep underground repositories is the prevailing proposed route for radioactive waste disposal. After the disposal of radioactive waste in the subsurface, H2 may be produced by corrosion of steel and, ultimately, radionuclides will be exposed to the surrounding environment. To evaluate the potential for microbial activities to impact disposal systems, we explored the microbial community structure and metabolic functions of a sediment-hosted ecosystem at the Horonobe Underground Research Laboratory, Hokkaido, Japan. Overall, we found that the ecosystem hosted organisms from diverse lineages, including many from the phyla that lack isolated representatives. The majority of organisms can metabolize H2, often via oxidative [NiFe] hydrogenases or electron-bifurcating [FeFe] hydrogenases that enable ferredoxin-based pathways, including the ion motive Rnf complex. Many organisms implicated in H2 metabolism are also predicted to catalyze carbon, nitrogen, iron and sulfur transformations. Notably, iron-based metabolism is predicted in a novel lineage of Actinobacteria and in a putative methane-oxidizing ANME-2d archaeon. We infer an ecological model that links microorganisms to sediment-derived resources and predict potential impacts of microbial activity on H2 consumption and retardation of radionuclide migration

This is a random community data set

### Design

DNA was extracted from the biomass retained on 0.2 um filters using the Extrap Soil DNA Kit Plus ver. 2 (Nippon Steel & Sumikin Eco-Tech Corporation) and sent for 150 bp paired-end sequencing with a 550 bp insert size by Hokkaido System Science Co., Ltd. using an Illumina HiSeq2000

### Runs

SRR3546457 SRR3546455 SRR3546454 SRR3546453 SRR3546452 SRR3546451 SRR3546450 SRR3546449


We downloaded from SRA like this

```
for SRR_ID in SRR3546457 SRR3546455 SRR3546454 SRR3546453 SRR3546452 SRR3546451 SRR3546450 SRR3546449; do
	fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip -N 5000 -X 255000 $SRR_ID
done
```

