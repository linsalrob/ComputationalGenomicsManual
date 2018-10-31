
# Metagenomics introduction

<img src="images/Staley Konopka.png" align="right" alt="Paper describing the great plate count anomaly"  title="Paper describing the great plate count anomaly" />

In this course we are going to focus on microbial metagenomics, using DNA sequencing to understand microbes in different environments. Microbes don’t usually live alone – even the simplest environments we’ve studied have multiple microbes present (Tyson et al., 2004; Edwards et al., 2006; Rodriguez-Brito et al., 2010)

<img src="images/plate count anomaly.gif" align="left" alt="The great plate count anomaly" title="The great plate count anomaly shows many more bacteria are detected by staining (Fig 1) than by culturing (Fig 2)"  />

The *Great Plate Count Anomaly* is that not every microbe that we can detect &mdash; e.g. via direct staining &mdash; can be grown on petri dishes. This anomaly was first described in 1985, although observations describing this phenomena date back to 1932 from freshwater habitats and 1959 from marine environments (as reported in Staley and Konopka).

For example, Staley and Konopka present this data on the growth of bacteria in Lake Washington. The two panels represent the number of bacteria that were detected by staining with acridine orange, and the number of bacteria that were detected with plate counts. The *x*-axis is the month of sampling, and the *y*-axis is the depth in the lake.

More recently, the role of the microbiome or microbiota has become dominant in everything from human health and disease to agriculture to lifestyles. 

The concept of the microbiome or microbiota is not new – it has been around for at least 60 years. However, the application of cheap DNA sequencing to microbial environments has completely altered our understanding of the three fundamental questions that metagenomics addresses: What is in the environment, what are they doing, and how are they doing it.

The first flurry of papers about metagenomics emerged in the first years of the 21<sup>st</sup> Century (Breitbart et al., 2002, 2003; Tyson et al., 2004; Venter et al., 2004)

Those early studies were mostly descriptive – providing an overview of the organisms in the environment being studied.  Comparative studies quickly followed, with studies comparing different environments, and meta-analyses of metagenomics data sets (Tringe et al., 2005; Rodriguez-Brito, Rohwer & Edwards, 2006; Dinsdale et al., 2008).

<img align="right" src="images/who what where.png" title='the fundamental questions for metagenomics - who what where how?' />

As sequencing became cheaper, data sets grew, and more samples could be sequenced for each environment. This lead to the ability to cross-compare similar samples from different environments, and through using savvy statistics we are able to capture complete genomes from
metagenomics samples.

Throughout the course, we will work through several of these analyses, from simple samples where we identify the organisms that are present, through combining complex samples to identify complete genomes.


## Metagenomics, BACs, and metabarcoding

<img src="images/Handelsman1.png" title="The paper descibing metagenomics" />
<img src="images/Handelsman2.png" align="left" title="The paper descibing metagenomics" />

The term metagenomics was original coined to describe cloning genes from the environment, particularly from soils (Handelsman et al., 1998)]

As shown in Fig. 3 of that paper, the concept is to clone fragments of DNA into bacterial artificial chromosomes, essentially low copy number plasmids that can hold large pieces of DNA and then use those fragments to select for a desired biological activity.

[Thought Experiment: Why use bacterial artificial chromosomes
that can hold long pieces of DNA rather than just clone short fragments
of DNA into, for example, high copy number vectors?]{.c2 .c10}

[]{.c6 .c2}

[Thought Experiment]{.c4}[: How could you find a piece of DNA encoding a
particular gene if that gene did not have a biological activity that you
could select for? (Clue: what is the difference between a genetic
selection and a genetic screen?)]{.c2 .c10}

[]{.c6 .c2}

[]{.c6 .c2}

[]{.c6 .c2}

------------------------------------------------------------------------

[]{.c6 .c2}

[16S Sequencing]{.c1} {#h.53qm4vu06m5o .c18}
=====================

[The small subunit of RNA polymerase contains both RNA and protein. This
complex is essential for life in DNA-based organisms as it is involved
in transcribing DNA into RNA, the first step in the central dogma. RNA
polymerase was probably one of the earliest complexes to evolve, soon
after nucleic acid chemistry evolved in the primordial soup. RNA
polymerase was probably predated by ribonucleotide reductase, which
converts RNA monomers to DNA monomers, and these were probably both
preceded by RNaseP, a ribozyme that catalyzes the cleavage of RNA. Even
the primordial soup was competitive: RNaseP probably evolved to
eliminate competitive RNA molecules, which then escaped by being
converted to DNA, but the DNA needed to get back to RNA to be
active!]{.c6
.c2}[![](images/image11.png)]{style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 349.33px; height: 461.12px;"}

[Carl Woese and colleagues at the University of Illinois in Urbana
Champaign were among the first to propose that the subunits of RNA
polymerase could be used for phylogenetic purposes ]{.c2}[[(Sogin, Sogin
&
Woese)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/7yAL&sa=D&ust=1540407111837000){.c7}]{.c2}[,
and using the 16S gene lead Woese and Fox to first recognize the three
domains of life, in what has been called the most important paper in
microbiology
(]{.c2}[[http://www.pnas.org/content/74/11/5088.full](https://www.google.com/url?q=http://www.pnas.org/content/74/11/5088.full&sa=D&ust=1540407111837000){.c7}]{.c0}[)
]{.c2}[[(Woese & Fox,
1977)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/YRty&sa=D&ust=1540407111837000){.c7}]{.c2}[.
In follow up papers, Woese and colleagues expanded their ideas and
demonstrated the relationship between the three domains e.g.
]{.c2}[[(Woese et al.
1990)](https://www.google.com/url?q=https://paperpile.com/c/p2c9eR/wzr4&sa=D&ust=1540407111837000){.c7}]{.c2
.c15}

[![](images/image4.png)]{style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 624.00px; height: 317.33px;"}

[Subsequently, Norman Pace (who earned his PhD from the University of
Illinois at Urbana Champaign in 1967) and colleagues used 16S sequences
to characterize microbes from the marine environment ]{.c2}[[(Schmidt,
DeLong & Pace,
1991)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/4ABJ&sa=D&ust=1540407111838000){.c7}]{.c2}[.
Soon after, both Jed Fuhrman and Ed DeLong identified the abundant
Archaea in the marine environment, suggesting that they were not
restricted to the extreme environments where they had originally been
isolated ]{.c2}[[(Fuhrman, McCallum & Davis, 1992; DeLong,
1992)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/3Ect%2BsUKy&sa=D&ust=1540407111838000){.c7}]{.c2}[.]{.c6
.c2}

[Sanger sequencing generated reads of approximately 750bp, and so two
primers were used to amplifiy and sequence the 16S gene, 27F – 1492R. As
described by Frank et al.: nearly full length 16S rRNA genes were
amplified using the 1492r primer (5′-TACCTTGTTACGACTT) and one of the
following three 27f primer formulations: twofold-degenerate primer
27f-CM (5′-AGAGTTTGATCMTGGCTCAG, where M is A or C), fourfold-degenerate
primer 27f-YM (5′-AGAGTTTGATYMTGGCTCAG, where Y is C or T), or
sevenfold-degenerate primer 27f-YM+3. The sevenfold-degenerate primer
27f-YM+3 is four parts 27f-YM, plus one part each of primers specific
for the amplification of Bifidobacteriaceae (27f-Bif,
5′-AGGGTTCGATTCTGGCTCAG), Borrelia (27f-Bor, 5′-AGAGTTTGATCCTGGCTTAG),
and Chlamydiales (27f-Chl, 5′-AGAATTTGATCTTGGTTCAG) sequences
]{.c2}[[(Frank et al.,
2008)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/xMsQ&sa=D&ust=1540407111838000){.c7}]{.c2}[.]{.c6
.c2}

[However, the use of shorter sequences, especially from the Illumina
platform, requires the use of primers that amplify a shorter region. The
current “standard” region that is amplified and sequenced is the
 515f-806r region. You can find the current recommended PCR protocol on
the Earth Microbiome Project website at
]{.c2}[[http://www.earthmicrobiome.org/emp-standard-protocols/16s/](https://www.google.com/url?q=http://www.earthmicrobiome.org/emp-standard-protocols/16s/&sa=D&ust=1540407111838000){.c7}]{.c0}[.
Their recommended sequencing primers are: GTGYCAGCMGCCGCGGTAA and
GGACTACNVGGGTWTCTAAT ]{.c2}[[(Caporaso et al., 2012; Apprill et al.,
2015)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/AVqE%2Bf1ln&sa=D&ust=1540407111839000){.c7}]{.c2}

[Thought experiment:]{.c4}[ ]{.c2}[Choose a 16S sequencing paper at
random, provided it has a principal components analysis (PCA) type of
analysis. How much of the variance does the 16S sequences explain?]{.c2
.c10}

[In this figure from Findley ]{.c2}[et al]{.c2 .c10}[ ]{.c2}[[(Findley
et al.,
2013)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/q6F1&sa=D&ust=1540407111839000){.c7}]{.c2}[,
the first axis explains 10% of the variation, the second axis explains
3.6% of the variation and the third axis explains 1.3% of the variation.
However, their data is comprised of 5,000 taxa measured from 14 skin
sites from 10 people with 3 different skin types. In other words, they
52 variables (the skin sites and types), 10 replicates (the people) and
are have 5,000 taxa that they are trying to use to explain the
variation. Most of the variation can be explained by just a few taxa.
]{.c6
.c2}[![](images/image8.png)]{style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 275.53px; height: 225.20px;"}

[]{.c6 .c2}

### [16S Databases]{.c20} {#h.3fk18cc5dpvs .c18}

[Several groups have generated databases of 16S sequences that you can
use to compare your fragments to]{.c6 .c2}

-   [Greengenes:
    ]{.c2}[[http://greengenes.lbl.gov/](https://www.google.com/url?q=http://greengenes.lbl.gov/&sa=D&ust=1540407111840000){.c7}]{.c0}[ ]{.c6
    .c2}
-   [SILVA – ARB:
    ]{.c2}[[http://www.arb-silva.de/](https://www.google.com/url?q=http://www.arb-silva.de/&sa=D&ust=1540407111840000){.c7}]{.c0}[ ]{.c6
    .c2}
-   [VAMPS:
    ]{.c2}[[http://vamps.mbl.edu/](https://www.google.com/url?q=http://vamps.mbl.edu/&sa=D&ust=1540407111841000){.c7}]{.c0}[ ]{.c6
    .c2}
-   [Ribosomal Database Project (RDP):
    ]{.c2}[[http://rdp.cme.msu.edu/](https://www.google.com/url?q=http://rdp.cme.msu.edu/&sa=D&ust=1540407111841000){.c7}]{.c0}[ ]{.c6
    .c2}

[Each database has pros and cons, and will give you slightly different
taxonomic resolution. However, for most next generation sequencing
studies, you get to phylum level, or thereabouts when you compare short
fragments of the 16S sequence with these databases.]{.c6 .c2}

[]{.c4 .c9}

[Thought Experiment]{.c4}[: What is the most different eukaryote from
humans you can envision that is in the same phylum as us (hint: we are
in the Chordata phylum)]{.c2 .c10}

[]{.c6 .c2}

[There have been several approaches to extrapolate from 16S sequences to
the functions that are present in the environment ]{.c2}[[(Langille et
al.,
2013)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/RFoN&sa=D&ust=1540407111842000){.c7}]{.c2}[,
however many independent studies have demonstrated that these
extrapolations do not capture the true scope of the functions in the
environment ]{.c2}[[(Steven et al., 2012; Haggerty & Dinsdale,
2017)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/BzS1%2BfKtR&sa=D&ust=1540407111842000){.c7}]{.c2}[.
The horizontal transfer of different metabolic genes between highly
related organisms compared to the slow rate of evolution via point
mutation obfuscate the metabolic differences between closely related
organisms and make it impossible to extrapolate based on phyla, class,
order, or family (see the discussion about mutation rates below).]{.c6
.c2}

[]{.c6 .c2}

[Metagenomics (Random Community Genomics)]{.c1} {#h.ampufu25gh5c .c18}
===============================================

[Isolating bacteria from the environment is hard, and scientists,
especially scientists that know how to program computers are inherently
lazy]{.c2}^[\[1\]](#ftnt1){#ftnt_ref1}^[. In the early 2000’s,
sequencing became cheap enough, and computer programs became good enough
that many scientists eschewed isolating bacteria for just sequencing
environmental DNA in bulk.]{.c6
.c2}[![](images/image9.png)]{style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 429.33px; height: 334.46px;"}

[It’s not entirely true that all scientists are lazy! There are many
cases where people were interested in sequences from the environment,
but those organisms were hard to isolate. For example, many of the early
metagenomics papers focused on viruses, because we had a very incomplete
notion of viruses in the environment, were not sure what their hosts
were, and did know how to culture them.  In those cases, especially,
sequencing environmental DNA in bulk and then using computers to solve
the problem became an obvious way to avoid isolation and culturing.]{.c6
.c2}

[However, sequencing is not free (yet) and the computational issues are
not solved. Back in the early days of metagenomics, sequencing was more
expensive and the computational issues had not begun to be addressed.
The SEED group at the Fellowship for the Interpretation of Genomes (FIG)
was a collective of researchers in Chicago, Wisconsin, San Diego, and
elsewhere. We were building the SEED platform and the Rapid Annotation
using Subsystems Technology (RAST) that were designed to annotated
complete microbial genomes (see the discussion of databases below). At
the same time, we were also sequencing environmental microbial samples
and trying to annotate them. By leveraging subsystems we were able to
reduce the complexity of the data: rather than trying to explain tens or
hundreds of thousands of protein functions, we could summarize those
into a few subsystems that we could understand.]{.c6 .c2}

[By the end of 2007 there were almost 100 random metagenomes that had
been sequenced using pyrosequencing. At that time, we were running
blastx against the NCBI non-redundant database (see databases below) and
summarizing the data using subsystems. In our discussions it became
obvious that this was the first time we had a consistent set of
metagenomes from diverse environments all treated the same way. We
reanalyzed all the data using the same version of the database,
annotated using the same approach, and compared all the subsystems
across metagenomes. This became the nine-biomes paper published in
Nature ]{.c2}[[(Dinsdale et al.,
2008)](https://www.google.com/url?q=https://paperpile.com/c/nxzdOK/TUkb&sa=D&ust=1540407111843000){.c7}]{.c2}[.
]{.c6
.c2}[![](images/image3.png)]{style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 365.33px; height: 298.45px;"}

[Notice that in this comparison the first two components of the
Canonical Discriminant Analysis (CDA; similar to PCA) explain \~80% of
the variance of the data. In this case, we have 30 variables (the level
1 subsystems) that we are using to explain the 9 biomes.]{.c6
.c2}[![](images/image1.png)]{style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 400.00px; height: 318.00px;"}

[At that time it was still expensive to sequence metagenomes, and so
most studies only sequenced a single point. Recently, however, it is so
cheap to sequence metagenomes that most studies sequence multiple
different samples from the environment. This has lead to the notion of
binning reads from related sequencing efforts to create population
genomes (see below).]{.c2}

------------------------------------------------------------------------

<div>

[\[1\]](#ftnt_ref1){#ftnt1}[This is not just our opinion: see, for
example, the three virtues of a great programmer:
http://threevirtues.com/ ]{.c6 .c11}

</div>
