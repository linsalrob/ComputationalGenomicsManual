
## Hands-on with Bin Chicken

### Bin Chicken installation and setup

See https://aroneys.github.io/binchicken for more details.

```bash
mamba create -n binchicken -c bioconda -c conda-forge 'binchicken>=0.12.5'
conda activate binchicken
binchicken build \
  --conda-prefix /storage/data/.conda \
  --singlem-metapackage /storage/data/metapackage \
  --checkm2-db /storage/data/checkm2
```

### Single-sample assembly with multi-sample binning

Start by running `binchicken single` to prepare the data to assemble each sample individually.

```bash
binchicken single \
    --forward 788707_20171213_S_R1.fastq.gz 788707_20180129_S_R1.fastq.gz 788707_20180313_S_R1.fastq.gz 788707_20181126_S_R1.fastq.gz \
    --reverse 788707_20171213_S_R2.fastq.gz 788707_20180129_S_R2.fastq.gz 788707_20180313_S_R2.fastq.gz 788707_20181126_S_R2.fastq.gz \
    --output single_assembly
```

The suggested assemblies with their respective binning samples can be found at `single_assembly/coassemble/target/elusive_clusters.tsv`.
In this case, only two of the samples are considered likely to recover genomes. These samples are 788707_20180313_S and 788707_20180129_S.
The other samples are probably too small (they are heavily subsampled) to recover genomes.

The actual assembly and binning can be run by adding `--run-aviary`.
Note that with 1 core, the assemblies will take ~30 minutes each.

```bash
binchicken single \
    --forward 788707_20171213_S_R1.fastq.gz 788707_20180129_S_R1.fastq.gz 788707_20180313_S_R1.fastq.gz 788707_20181126_S_R1.fastq.gz \
    --reverse 788707_20171213_S_R2.fastq.gz 788707_20180129_S_R2.fastq.gz 788707_20180313_S_R2.fastq.gz 788707_20181126_S_R2.fastq.gz \
    --output single_assembly --run-aviary --cores 5
```

The assembly and binning for each sample is found at `single_assembly/coassemble/coassemble/`.
Each sample should have a folder containing `assemble` for the assembly and `recover` for the binning.
The bins for each sample are found in `recover/bins`, with genome info at `recover/bins/bin_info.tsv`.

The recovered bins are likely only ~40-60% complete, with fairly high contamination.
This is probably due the small sample size, but the genome could still be analysed further with e.g. GTDBtk to find out their taxonomy.

### Coassembly with multi-sample binning

Now that we have run single-sample assembly for the decent samples, we can run coassembly across the dataset.
Because 

```bash
binchicken coassemble \
    --forward 788707_20171213_S_R1.fastq.gz 788707_20180129_S_R1.fastq.gz 788707_20180313_S_R1.fastq.gz 788707_20181126_S_R1.fastq.gz \
    --reverse 788707_20171213_S_R2.fastq.gz 788707_20180129_S_R2.fastq.gz 788707_20180313_S_R2.fastq.gz 788707_20181126_S_R2.fastq.gz \
    --output coassembly --max-coassembly-samples 5
```

The suggested coassemblies with their respective binning samples can be found at `coassembly/coassemble/target/elusive_clusters.tsv`.
Since they share single-copy marker genes, the samples 788707_20180129_S and 788707_20180313_S are suggested for coassembly.

We can run the actual coassembly and binning as before by adding `--run-aviary`.
Note that with 1 core, the coassembly will take ~1 hour.

```bash
binchicken coassemble \
    --forward 788707_20171213_S_R1.fastq.gz 788707_20180129_S_R1.fastq.gz 788707_20180313_S_R1.fastq.gz 788707_20181126_S_R1.fastq.gz \
    --reverse 788707_20171213_S_R2.fastq.gz 788707_20180129_S_R2.fastq.gz 788707_20180313_S_R2.fastq.gz 788707_20181126_S_R2.fastq.gz \
    --output coassembly --max-coassembly-samples 5 --run-aviary --cores 5
```

## References

- Alneberg, J., Bjarnason, B.S., de Bruijn, I., Schirmer, M., Quick, J., Ijaz, U.Z., Lahti, L., Loman, N.J., Andersson, A.F., Quince, C., 2014. Binning metagenomic contigs by coverage and composition. Nat. Methods 11, 1144–1146. https://doi.org/10.1038/nmeth.3103
- Ayling, M., Clark, M.D., Leggett, R.M., 2020. New approaches for metagenome assembly with short reads. Brief. Bioinform. 21, 584–594. https://doi.org/10.1093/bib/bbz020
- Kang, D.D., Froula, J., Egan, R., Wang, Z., 2015. MetaBAT, an efficient tool for accurately reconstructing single genomes from complex microbial communities. PeerJ. https://doi.org/10.7717/peerj.1165
- Kang, D.D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., Wang, Z., 2019. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ 7, e7359. https://doi.org/10.7717/peerj.7359
- Li, D., Liu, C.-M., Luo, R., Sadakane, K., Lam, T.-W., 2015. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics 31, 1674–1676. https://doi.org/10.1093/bioinformatics/btv033
- Mallawaarachchi, V., Wickramarachchi, A., Xue, H., Papudeshi, B., Grigson, S.R., Bouras, G., Prahl, R.E., Kaphle, A., Verich, A., Talamantes-Becerra, B., Dinsdale, E.A., Edwards, R.A., 2024. Solving genomic puzzles: computational methods for metagenomic binning. Brief. Bioinform. 25, bbae372. https://doi.org/10.1093/bib/bbae372
- Nissen, J.N., Johansen, J., Allesøe, R.L., Sønderby, C.K., Armenteros, J.J.A., Grønbech, C.H., Jensen, L.J., Nielsen, H.B., Petersen, T.N., Winther, O., Rasmussen, S., 2021. Improved metagenome binning and assembly using deep variational autoencoders. Nat. Biotechnol. 39, 555–560. https://doi.org/10.1038/s41587-020-00777-4
- Nurk, S., Meleshko, D., Korobeynikov, A., Pevzner, P.A., 2017. metaSPAdes: a new versatile metagenomic assembler. Genome Res. 27, 824–834. https://doi.org/10.1101/gr.213959.116
- Pan, S., Zhao, X.-M., Coelho, L.P., 2023. SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. Bioinformatics 39, i21–i29. https://doi.org/10.1093/bioinformatics/btad209
- Pan, S., Zhu, C., Zhao, X.-M., Coelho, L.P., 2022. A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. Nat. Commun. 13, 2326. https://doi.org/10.1038/s41467-022-29843-y
- Sieber, C.M.K., Probst, A.J., Sharrar, A., Thomas, B.C., Hess, M., Tringe, S.G., Banfield, J.F., 2018. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nat. Microbiol. 3, 836–843. https://doi.org/10/gfwwfg
- Wang, Z., You, R., Han, H., Liu, W., Sun, F., Zhu, S., 2024. Effective binning of metagenomic contigs using contrastive multi-view representation learning. Nat. Commun. 15, 585. https://doi.org/10.1038/s41467-023-44290-z
- Wu, Y.-W., Simmons, B.A., Singer, S.W., 2016. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinforma. Oxf. Engl. 32, 605–607. https://doi.org/10.1093/bioinformatics/btv638
