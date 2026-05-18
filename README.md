# tasaR

## Toolbox for Amplicon Sequencing Analysis in R

# Note!

#### This is the early release of verion 0.2.0. Functionality has been re-implemented using S4 objects and an improved and more portable workflow for interactive use. PANDAseq merging has been vendored into this version, but portability and cross-platform functionality are still to be verified and implemented. A GUI wrapper in Shiny is also still to come for this version.


#### About

This software was developed to facilitate analysis of mutations within targeted amplicon sequences by Illumina deep sequencing. Most existing NGS analysis packages are more focused on aligning and mapping large populations of diverse sequences over an in-depth analysis of a single reference sequence. tasaR was built to close that gap -- facilitating close analysis of single amplicon sequences.

The initial use case for tasaR was to analyze somatic hypermutation of engineered antibody sequences. Mutations within a reference starting sequence are mapped, AID hotspot motifs are parsed, and the distribution of mutations across an antibody sequence are presented in a variety of outputs to facilitate downstream analyses. As a side effect, tasaR is also an effective software to analyze the activity of a targeted nuclease like Cas9. Specific modules build for paired amplicon analysis establish the reference sequence in a matched control sample and measure and classify rates and types of DNA mutations produced by the nuclease. Unlike other software such as CRISPResso2, tasaR identifies non-indel mutations such as base changes that can be introduced by gene conversion copying from closely-related homeologous sequences in the genome.

To facilitate paired amplicon analysis, tasaR also includes a wrapper for the paired sequence alignment software [PANDAseq](https://github.com/neufeld/pandaseq). This facilitates local sequence analysis in a cross-platform manner, and allows tasaR to serve as an all-in-one amplicon sequencing analysis platform.

Created by Geoffrey L. Rogers, PhD.
Department of Immunology and Immune Therapeutics, Keck School of Medicine, University of Southern California, Los Angeles, CA, USA

## **Dependencies**

R (tested in v4.5.2)

Depends: data.table

Imports: stringr, ggplot2, ggseqlogo, R.utils, zip, markdown, BiocGenerics, S4Vectors, ShortRead, ggh4x, ggmsa, methods, msa, Biostrings, pwalign, purrr

To build PANDAseq from source also requires the following system libraries: Standard math library (-lm),  Libtool Dynamic Loader (-lltdl),  Bzip2 (-lbz2), and Zlib (-lz).

## **Quick Start Tips**

* Run the command tasar::makeSettingsCSV(filedest) to write a template .csv file to the specified directory.

* Fill out the .csv file in your preferred spreadsheet editor. See later in this file for column descriptions.

* To measure mutations in a sequence, run the command tasaR::analyzeAmplicon(Settings.csv, row = x) for a single sample. Multiple samples can be analyzed at once using the function batchAnalyzeAmplicon().

* Nuclease activity can be measured using the function pairedAnalyzeAmplicon(). See help file for parameters with ?pairedAnalyzeAMplicon.

* Results can be exported with the functions exportTables, batchSummarize, and exportNucleaseAnalysis.

* Paired end reads can be merged with the function pandaseq_merge_files(). [Early implementation, portability has not been validated.]

## Description of the Input.csv file

The Input.csv file specifies the analysis parameters for each sample used by the software. Each column is a parameter, and each row represents a different sample.

**SampleName -** Name identifier of each sample. Please avoid spaces in chosen names. Hyphens (or underscores) are appropriate in lieu of spaces.
<br></br>
**IsAntibody -** Input TRUE or FALSE to identify whether the sequence is an antibody or not. If so, input coordinates of the regions to be analyzed at the end of the file in the FR1-FR4 columns.
<br></br>
**MergedFASTQFileName -** File name of the merged .fastq.gz file if paired end merging is performed outside of tasaR. If using the scripts for PANDAseq merging, this column should be left blank and will be automatically filled by the software.
<br></br>
**ReferenceSequence -** The DNA sequence of your antibody (or other target sequence) to be analyzed by the software. If analyzing an antibody, it is recommended to omit the signal peptide portion of the sequence per convention.
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
**AmpliconLength -** The total length of your amplicon between the Illumina adapters but including the 5' barcode and 3' UMI sequences.
<br></br>
**InsertStart -** Identifies the start site of the desired sequence to be analyzed. For antibodies, this should be after the leader sequence at the start of the variable domain. To calculate: measure the number of nt from the beginning of the binding site of your forward primer on the top strand <u>through (including) the first nt of the desired sequence</u>. Add 7 (length of the barcode) to this measurement to get the final value of InsertStart.
<br></br>
**InsertEnd -** Identifies the end site of the desired sequence to be analyzed, as measured from the 3' end.  To calculate: measure the number of nt from the beginning of the binding site of your reverse primer on the bottom strand <u>through (including) the first nt of the desired sequence.</u> Add 12 (length of the UMI) to this measurement and make it negative to get the final value of InsertEnd.
<br></br>
Note: The **InsertStart** and **InsertEnd** coordinates are not affected by indels, since they are measured from each end of the sequence and reads are filtered and validated based on intact primer ends with appropriate barcodes and UMIs.
<br></br>
**FR1Start -** The number of the first nt of the framework region 1 (FR1) within the sequence defined in **ReferenceSequence**. Should most likely be 1.
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

## **Results Output**

* Use exportTables to export .csv files of the analysis results.

* Use batchSummarzie to export summary .csv files and graphs comparing the overall results of samples in the batch.

* Use exportNucleaseAnalysis to export sequence alignment pdfs and graphs showing mutations around the nuclease cut site. Note: generating the alignment pdfs requires an installation of LaTeX on the system. If missing, load the package tinytex and run the command tinytex::install_tinytex().

## **Change log**

**v0.2.0** - 5/17/2026 - minimal release of refactored package with improved usability employing S4 objects and improved portable workflows

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
