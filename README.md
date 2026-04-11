# Targeted Amplicon Sequencing Analyzer

#### About

This software was developed to facilitate analysis of mutations within targeted amplicons sequences by Illumina deep sequencing. This was originally developed to facilitate analysis of somatic hypermutation (SHM) of specific antibody sequences in B cells introduced by site-specific genome editing, as other tools focused on this sort of analysis look at the population level rather than a specific starting sequence.

tasAnalyzer can also be used to analyze the DNA mutations introduced at a specific site by a targeted nuclease. During my experiments, I found that some loci produce gene editing outcomes beyond simple insertions and deletions (indels). For loci that have close homeologs, gene conversion-like repair may result in the transfer of sequences from these genes to the targeted gene, likely via homology-directed repair. Other analysis softwares, both for Sanger sequencing and Illumina sequencing, miss these events and can drastically underestimate the rate of DNA change at these loci compared to when tasAnalyzer accounts for these gene coversion-like events.

Created by Geoffrey L. Rogers, PhD.
Department of Immunology and Immune Therapeutics, Keck School of Medicine, University of Southern California, Los Angeles, CA, USA

## **System Requirements**

#### To run locally:

R (tested in v4.5.2)

CRAN packages: BiocManager, R.utils, zip, data.table, openxlsx, stringr, ggplot2, ggseqlogo, ggrepel, seqinr, ape

Bioconductor packages: Biostrings, ShortRead, pwalign, msa, ggmsa

CRAN and Bioconductor packages will be automatically installed and loaded.

To merge paired end reads from Illumina sequencing, this software also requires installation of **PANDAseq** ([GitHub - neufeld/pandaseq: PAired-eND Assembler for DNA sequences](https://github.com/neufeld/pandaseq)).

## **Quick Start Guide**

* Specifications and options are provided using a spreadsheet template (Input.csv). A template is included in the inst folder. DO NOT delete any columns; for any unused values, leave those cells blank.
<br></br>
* Input the required details for each sample to be analyzed in the Input.csv file, then run tas_analyze() on that file. Settings optionally can be enabled with TRUE/FALSE arguments in the function.
<br></br>
* Outputs will be compressed into a .zip archive that can be downloaded from the "Download" tab.

## **Options for tas_analyze**

**Merge paired-end reads? -** Tells tasAnalyzer whether to call PANDAseq to merge paired-end reads.

**Enable somatic hypermutation module? -** Tells tasAnalyzer whether to identify AID hotspot motifs and measure mutations at those sites.

**Measure amino acid mutations? -** Tells tasAnalyzer whether to identify and classify mutations of the amino acid sequence.

**Number of sequences to align? -** For multiple sequence alignments, tells how many sequences to align. A default of 10 is recommended for readability and processing time, as this step can be quite intensive.

**Minimum read frequency for analysis? - Only for samples without UMIs.** The minimum frequency (%) of reads to be included in the analysis. A range of 1% (0.1) to 0.01% (0.0001) is recommended.

**Predict DNA repair pathway use? -** Tells tasAnalyzer whether to measure the frequences of different types of mutations and predict the rate of usage for the likely underlying DNA repair pathways. Useful for amplicons spanning gene editing nuclease target sites to estimate the activity of the nuclease.

**Plot phylogenetic trees of multiple sequence alignments? -** Tells tasAnalyzer whether to plot phylogenetic trees of the sequences aligned in the MSA module.

## Description of the Input.csv file

The Input.csv file specifies the analysis parameters for each sample used by the software. Each column is a parameter, and each row represents a different sample.

**SampleName -** Name identifier of each sample. Please avoid spaces in chosen names. Hyphens (or underscores) are appropriate in lieu of spaces.
<br></br>
**ForwardFASTQFileName -** File name of the R1 Illumina sequencing file for that sample. Files in .fastq.gz format are much smaller than uncompressed .fastq files and access speeds are comparable. 
<br></br>
**ReverseFASTQFileName -** File name of the R2 Illumina sequencing file for that sample. Files in .fastq.gz format are much smaller than uncompressed .fastq files and access speeds are comparable.
<br></br>
**MergedFASTQFileName -** File name of the merged .fastq.gz file if paired end merging is performed outside of tasAnalyzer. If using the scripts for PANDAseq merging, this column should be left blank and will be automatically filled by the software.
<br></br>
**ForwardExtensionType -** The type of extension on the forward primer. Can be either "Barcode" or "UMI".
<br></br>
**ForwardExtension -** DNA sequence (5' - 3') of the forward extension. Supports standard ambiguity codes.
<br></br>
**ForwardPrimer -** DNA sequence (5' - 3') of the forward primer used for amplification. In this box, only include the <u>sequence that binds to the target sequence</u>.
<br></br>
**ReverseExtensionType -** The type of extension on the reverse primer. Can be either "Barcode" or "UMI".
<br></br>
**ReverseExtension -** DNA sequence (5' - 3') of the reverse extension. Supports standard ambiguity codes.
<br></br>
**ReversePrimer -** DNA sequence (5' - 3') of the reverse primer used for amplification. In this box, only include the <u>sequence that binds to the target sequence</u>.
<br></br>
**ReferenceSequence -** The DNA sequence of your antibody (or other target sequence) to be analyzed by the software. If analyzing an antibody, it is recommended to omit the signal peptide portion of the sequence per convention.
<br></br>
**AmpliconLength -** The total length of your amplicon between the Illumina adapters but including the 5' barcode and 3' UMI sequences.
<br></br>
**MaxDeletion -** The maximum deletion size allowed in merged paired end reads by pandaseq. A default of 25 is recommended.
<br></br>
**MaxInsertion -** The maximum insertion size allowed in merged paired end reads by pandaseq. A default of 25 is recommended.
<br></br>
**InsertStart -** Identifies the start site of the desired sequence to be analyzed. For antibodies, this should be after the leader sequence at the start of the variable domain. To calculate: measure the number of nt from the beginning of the binding site of your forward primer on the top strand <u>through (including) the first nt of the desired sequence</u>. Add 7 (length of the barcode) to this measurement to get the final value of InsertStart.
<br></br>
**InsertEnd -** Identifies the end site of the desired sequence to be analyzed, as measured from the 3' end.  To calculate: measure the number of nt from the beginning of the binding site of your reverse primer on the bottom strand <u>through (including) the first nt of the desired sequence.</u> Add 12 (length of the UMI) to this measurement and make it negative to get the final value of InsertEnd.
<br></br>
Note: The **InsertStart** and **InsertEnd** coordinates are not affected by indels, since they are measured from each end of the sequence and reads are filtered and validated based on intact primer ends with appropriate barcodes and UMIs.
<br></br>
**Antibody -** If your sequence is an antibody, put "Yes" in this field. You will then input the coordinates of the different antibdoy regions for analysis. To omit these analyses, put "No" in this column.
<br></br>
**FR1Start -** The number of the first nt of the framework region 1 (FR1) within the sequence defined in **ReferenceSequence**. Should be 1.
<br></br>
**CDR1Start -** The number of the first nt of the complementarity determining region 1 (CDR1) within the sequence defined in **ReferenceSequence**.**FR2Start -**
<br></br>
**FR2Start -** The number of the first nt of the FR2 region within the sequence defined in **ReferenceSequence**.
<br></br>
**CDR2Start -** The number of the first nt of the CDR2 region within the sequence defined in **ReferenceSequence**.
<br></br>
**FR3Start -** The number of the first nt of the FR3 region within the sequence defined in **ReferenceSequence**.
<br></br>
**CDR3Start -** The number of the first nt of the CDR3 region within the sequence defined in **ReferenceSequence**.
<br></br>
**FR4Start -** The number of the first nt of the FR4 region within the sequence defined in **ReferenceSequence**.

## **Tool Output**

The tool outputs a variety of files, tables, and graphs that may be of use to the user:

#### If merging was performed by tasAnalyzer

* Merged .fastq.gz files for each sample (output of PANDAseq).
<br></br>
* PANDAseq logs and filtering statistics showing read counts after each step and UMI counts for each sample.

#### For all applications

* Filtered merged .fastq.gz files for each sample after applying barcode and UMI intergrity filters.

#### Results

* **Sequences.xlsx -** A table showing all sequences found in the sample. If appropriate, sequences are binned by UMI. Includes read/UMI counts, frequencies, DNA mutations, amino acid sequence, and amino acid mutations.
<br></br>
* **Mutations.xlsx -** Quantification of the amount of mutation observed at each nucleotide along the sequence (in %). MutAll includes all nts in the ReferenceSequence. Other measures require measure.shm = TRUE. MutCyt pulls out the cytosines of AID hotspot motifs (WR**C**H), whereas MutNonC is all other nts in the sequence. MotifSums shows calculated % of mutation across different target motifs and denominantors, as described by the row titles.
<br></br>
* **WRCH/WRCY tables.xlsx -** Created if the SHM module is enabled. Tables showing each of the AID hotspots identified by tasAnalyzer and the mutation frequency at that site. One tab for each sample.
<br></br>
* **Mutation Types.xlsx -** Created if the DNA Repair Pathway module is enabled. For each sample, classifies and counts the frequency of mutations based on the predicted underlying DNA repair pathway. Non-homologous end joining (NHEJ): insertions and -1 or -2 deletions. Microhomology-mediated end joining (MMEJ): deletions > -2. Base change: changes in the sequence without indels. Indel + Base change: Sequencing with both 1 or more indels and base change outcomes. Other: Sequences not falling into any other category.

#### Graphs

* Bar charts summarizing the percent of mutation at AID cytosines or at all nts for all samples analyzed.
<br></br>
* Box plots for each sample showing the mutation frequency of all AID cytosines vs. all other nts at the DNA or protein level. For proteins, if any nt in the codon is an AID cytosine, that residue is considered an AID cytosine. No distinction is made if the codon includes more than one AID cytosine.
<br></br>
* Mutagenesis bar charts for each sample showing the % mutation at every nt (or aa) in the sequence.
<br></br>
* Mutagenesis bar charts with the residues that are AID cytosines labeled with a red C above each bar.
<br></br>
* Histograms showing the distribution of the number of mutations per read after UMI binning.

#### MSA

* For each sample, a multiple sequence alignment (ClustalOmega) of the top 10 most common DNA or protein sequences after UMI normalization.
<br></br>
* Phylogenetic distance tree of the top 10 most common DNA sequences after UMI normalization, as above (if PhyloTree = TRUE).

#### Seqlogo

* Sequence logo plots. If the sequence was defined as an antibody in Input.csv, sequence logo plots are generated for each identified region of the antibody (CDRs and FR regions).

## **Change log**

**v0.1** - 4/10/2026 - snapshot version used in Huang, ..., Rogers, and Cannon, Nat Commun 2026

## **Appendix A: Experimental Procedures**

The preferred methods for sequencing analysis with this software are first described in Rogers et al. Nat. Biomed. Eng. 2024 Dec;8(12):1700-1714. doi: 10.1038/s41551-024-01240-4. Briefly, the region surrounding the sequence to be characterized is amplified by PCR, and UMIs were introduced in a single-step PCR reaction by dilution alongside a forward primer and a high-concentration amplificiation primer. Targeted amplicon sequencing is then performed by Illumina, with machine specs dependent on the length of the target sequence. 2x250 MiSeq sequencing can generate amplicons with good overlap that are able to span an antibody variable region. Barcodes and UMIs are optional and can be specified by the Input.csv file.

## **Appendix B: Primer Design**

###### Primers

The forward primer is laid out as follows:

PARTIAL ILLUMINA ADAPTER - [ForwardExtension] - ForwardPrimer

        ACACTCTTTCCCTACACGACGCTCTTCCGATCT[ForwardExtension]nnnnnnnnn

The reverse primer is laid out as follows:

PARTIAL ILLUMINA ADAPTER - [ReverseExtension] - ReversePrimer

        GACTGGAGTTCAGACGTGTGCTCTTCCGATCT[ReverseExtension]nnnnnn

The extensions on either primer are optional. These can be a barcode for sample indexing or a unique molecular identifier (UMI).

The software expects at most 1 UMI; dual UMIs are not currently supported.

A  7 bp barcode sequence is suggested for identification and enhanced demultiplexing. A list of 32 barcodes with good distance between them was subsetted from previously described methodology:

ATCGATT, CAGTCAA, GCTAGCC, TGACTGG, CAATTGC, GCCAATG, TGGCCAT, ATTGGCA, CGATGTA, GTCATAC, TAGCACG, ACTGCGT, GGTATCG, TTACAGT, AACGCTA, CCGTGAC, GTGATCT, TATCAGA, ACAGCTC, ACGGTCT, CGTTAGA, GTAACTC, TACCGAG, TACGTTC, ACGTAAG, CGTACCT, CTACTCG, GACGAGT, TCGTCTA, AGTAGAC, GGAGTAC, TTCTACG

A patterned 12 bp UMI [NNNYRNNNYRNN] is suggested was added to each sample to control for PCR bias/errors and enhance sequence precision and quantification. Using a defined pattern aids removal of reads with sequencing errors that could have their UMI mis-attributed, as previously described.

