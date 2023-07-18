# Features allowed in INSDC sequence files.

These are the features that are allowed in DDBJ/EMBL/GenBank format files. Any other features are not valid in those files.

Feature Key | Description
--- | ---
assembly_gap | gap between two components of a genome or transcriptome assembly
C_region | constant region of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains; includes one or more exons depending on the particular chain
CDS | coding sequence; sequence of nucleotides that corresponds with the sequence of amino acids in a protein (location includes stop codon); feature includes amino acid conceptual translation.
centromere | region of biological interest identified as a centromere and which has been experimentally characterized
D-loop | displacement loop; a region within mitochondrial DNA in which a short stretch of RNA is paired with one strand of DNA, displacing the original partner DNA strand in this region; also used to describe the displacement of a region of one strand of duplex DNA by a single stranded invader in the reaction catalyzed by RecA protein
D_segment | Diversity segment of immunoglobulin heavy chain, and T-cell receptor beta chain
exon | region of genome that codes for portion of spliced mRNA, rRNA and tRNA; may contain 5'UTR, all CDSs and 3' UTR
gap | gap in the sequence
gene | region of biological interest identified as a gene and for which a name has been assigned
iDNA | intervening DNA; DNA which is eliminated through any of several kinds of recombination
intron | a segment of DNA that is transcribed, but removed from within the transcript by splicing together the sequences (exons) on either side of it
J_segment | joining segment of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains
mat_peptide | mature peptide or protein coding sequence; coding sequence for the mature or final peptide or protein product following post-translational modification; the location does not include the stop codon (unlike the corresponding CDS)
misc_binding | site in nucleic acid which covalently or non-covalently binds another moiety that cannot be described by any other binding key (primer_bind or protein_bind)
misc_difference | feature sequence is different from that presented in the entry and cannot be described by any other difference key (old_sequence, variation, or modified_base)
misc_feature | region of biological interest which cannot be described by any other feature key; a new or rare feature
misc_recomb | site of any generalized, site-specific or replicative recombination event where there is a breakage and reunion of duplex DNA that cannot be described by other recombination keys or qualifiers of source key (/proviral)
misc_RNA | any transcript or RNA product that cannot be defined by other RNA keys (prim_transcript, precursor_RNA, mRNA, 5'UTR, 3'UTR, exon, CDS, sig_peptide, transit_peptide, mat_peptide, intron, polyA_site, ncRNA, rRNA and tRNA)
misc_structure | any secondary or tertiary nucleotide structure or conformation that cannot be described by other Structure keys (stem_loop and D-loop)
mobile_element | region of genome containing mobile elements
modified_base | the indicated nucleotide is a modified nucleotide and should be substituted for by the indicated molecule (given in the mod_base qualifier value)
mRNA | messenger RNA; includes 5'untranslated region (5'UTR), coding sequences (CDS, exon) and 3'untranslated region (3'UTR)
ncRNA | a non-protein-coding gene, other than ribosomal RNA and transfer RNA, the functional molecule of which is the RNA transcript
N_region | extra nucleotides inserted between rearranged immunoglobulin segments.
old_sequence | the presented sequence revises a previous version of the sequence at this location
operon | region containing polycistronic transcript including a cluster of genes that are under the control of the same regulatory sequences/promoter and in the same biological pathway
oriT | origin of transfer; region of a DNA molecule where transfer is initiated during the process of conjugation or mobilization
polyA_site | site on an RNA transcript to which will be added adenine residues by post-transcriptional polyadenylation
precursor_RNA | any RNA species that is not yet the mature RNA product; may include ncRNA, rRNA, tRNA, 5' untranslated region (5'UTR), coding sequences (CDS, exon), intervening sequences (intron) and 3' untranslated region (3'UTR)
prim_transcript | primary (initial, unprocessed) transcript; may include ncRNA, rRNA, tRNA, 5' untranslated region (5'UTR), coding sequences (CDS, exon), intervening sequences (intron) and 3' untranslated region (3'UTR)
primer_bind | non-covalent primer binding site for initiation of replication, transcription, or reverse transcription; includes site(s) for synthetic e.g., PCR primer elements
propeptide | propeptide coding sequence; coding sequence for the domain of a proprotein that is cleaved to form the mature protein product.
protein_bind | non-covalent protein binding site on nucleic acid
regulatory | any region of sequence that functions in the regulation of transcription, translation, replication or chromatin structure
repeat_region | region of genome containing repeating units
rep_origin | origin of replication; starting site for duplication of nucleic acid to give two identical copies
rRNA | mature ribosomal RNA; RNA component of the ribonucleoprotein particle (ribosome) which assembles amino acids into proteins.
S_region | switch region of immunoglobulin heavy chains; involved in the rearrangement of heavy chain DNA leading to the expression of a different immunoglobulin class from the same B-cell
sig_peptide | signal peptide coding sequence; coding sequence for an N-terminal domain of a secreted protein; this domain is involved in attaching nascent polypeptide to the membrane leader sequence
source | identifies the biological source of the specified span of the sequence; this key is mandatory; more than one source key per sequence is allowed; every entry/record will have, as a minimum, either a single source key spanning the entire sequence or multiple source keys, which together, span the entire sequence.
stem_loop | hairpin; a double-helical region formed by base-pairing between adjacent (inverted) complementary sequences in a single strand of RNA or DNA.
STS | sequence tagged site; short, single-copy DNA sequence that characterizes a mapping landmark on the genome and can be detected by PCR; a region of the genome can be mapped by determining the order of a series of STSs
telomere | region of biological interest identified as a telomere and which has been experimentally characterized
tmRNA | transfer messenger RNA; tmRNA acts as a tRNA first, and then as an mRNA that encodes a peptide tag; the ribosome translates this mRNA region of tmRNA and attaches the encoded peptide tag to the C-terminus of the unfinished protein; this attached tag targets the protein for destruction or proteolysis
transit_peptide | transit peptide coding sequence; coding sequence for an N-terminal domain of a nuclear-encoded organellar protein; this domain is involved in post-translational import of the protein into the organelle
tRNA | mature transfer RNA, a small RNA molecule (75-85 bases long) that mediates the translation of a nucleic acid sequence into an amino acid sequence
unsure | a small region of sequenced bases, generally 10 or fewer in its length, which could not be confidently identified. Such a region might contain called bases (A, T, G, or C), or a mixture of called-bases and uncalled-bases ('N'). The unsure feature should not be used when annotating gaps in genome assemblies. Please refer to assembly_gap feature for gaps within the sequence of an assembled genome. For annotation of gaps in other sequences than assembled genomes use the gap feature.
V_region | variable region of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains; codes for the variable amino terminal portion; can be composed of V_segments, D_segments, N_regions, and J_segments
V_segment | variable segment of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains; codes for most of the variable region (V_region) and the last few amino acids of the leader peptide
variation | a related strain contains stable mutations from the same gene (e.g., RFLPs, polymorphisms, etc.) which differ from the presented sequence at this location (and possibly others)
3'UTR | 1) region at the 3' end of a mature transcript (following the stop codon) that is not translated into a protein; 2) region at the 3' end of an RNA virus (following the last stop codon) that is not translated into a protein
5'UTR | 1) region at the 5' end of a mature transcript (preceding the initiation codon) that is not translated into a protein; 2) region at the 5' end of an RNA virus genome (preceding the first initiation codon) that is not translated into a protein
