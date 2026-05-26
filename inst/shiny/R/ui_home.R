home_page_ui <- function() {
               nav_panel(title = "",
                         icon = icon("home"),
                        value = "splash.app",
                        fluidRow(
                          div(
                            class = "hero-banner",
  
                            div(
                              class = "hero-content",
  
                              div(
                                class = "hero-title",
                                "tasaR"
                              ),
  
                              div(
                                class = "hero-subtitle",
                                "Toolbox for Amplicon Sequencing Analysis in R"
                              )
                            )
                          )
                        ),
                        fluidRow(
                          div(
                            class = "upload-card-shell",
                            card(
                              card_header(
                                "What is tasaR?",
                                class = "card-header-custom"
                              ),
                              card_body(
                                class = "card-text-body",
                                HTML("
      <div class='doc-block'>
        <h3>Overview</h3>
        <div style='height: 4px;'></div>
        <p>
          tasaR is an R package built to support the analysis of Illumina targeted amplicon sequencing experiments.
          This software performs in-depth analysis of single amplicons across multiple modules. tasaR provides 3 key functionalities:
        </p>
        <div style='height: 8px;'></div>
        <ul>
          <li>
            <strong>Merge Reads:</strong> An integrated workflow to merge paired end reads using the software PANDAseq.<sup>1</sup>
          </li>
          <li>
            <strong>Measure Mutations:</strong> Analyzes target sequences and quantifies mutagenesis across the reference sequence.
            Somatic hypermutation is measured by annotating AID hotspot cytosines and quantifying mutation specifically at those sites.
          </li>
          <li>
            <strong>Nuclease Analysis:</strong> Quantifies the DNA mutations created by a targeted nuclease at a specific site in the genome.
            Uses a matched control sequence to determine the reference sequence; supports heterozygous genomic DNA.
          </li>
        </ul>
  
        <div style='height: 12px;'></div>
  
        <h3>Features</h3>
        <div style='height: 4px;'></div>
        <ul>
          <li>Automated quality-aware merging of paired end reads.</li>
          <li>Quantify DNA and protein mutations across the sequence.</li>
          <li>Annotate AID hotspot motifs and quantify somatic hypermutation.</li>
          <li>Analyze indels and base change events around a nuclease cut site.</li>
          <li>Predict the usage of DNA repair pathways after nuclease-induced DNA scarring.<sup>2</sup></li>
          <li>Supports unique molecular identifier (UMI) tagging and barcodes to enhance sequence de-multiplexing.</li>
        </ul>
        <div style='height: 4px;'></div>
  
  
            <div class='footnote-divider'></div>
  
      <div class='footnotes'>
  
        <p>
          <sup>1</sup>
          Masella AP, et al. (2012). PANDAseq: paired-end assembler for Illumina sequences. BMC Bioinformatics 14:13:31.
        </p>
  
        <p>
          <sup>2</sup>
          Tatiossian KJ, et al. (2021). Rational Selection of CRISPR-Cas9 Guide RNAs for Homology-Directed Genome Editing. Mol Ther 3;29(3):1057-1069.
        </p>
  
      </div>
  
      </div>
    ")
  
                              )
                            )
                          )
                        )
                        )
}
