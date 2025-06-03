# This only plots transcripts with CDS now


# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggtranscript)
library(rtracklayer)
library(here)

# Arguments ---------------------------------------------------------------

APLP1_args <- list(
  gtf = here::here("raw_data", "IsoSeq", "APLP1_and_APLP2", "APLP1_corrected.gtf.cds.gff"),
  tx = here::here("raw_data", "IsoSeq", "APLP1_and_APLP2", "APLP1_classification_filtered.txt")
)
  
APLP2_args <- list(
  gtf = here::here("raw_data", "IsoSeq", "APLP1_and_APLP2", "APLP2_corrected.gtf.cds.gff"),
  tx = here::here("raw_data", "IsoSeq", "APLP1_and_APLP2", "APLP2_classification_filtered.txt")
)

# Load data ---------------------------------------------------------------

# GTF files
APLP1_gtf <- rtracklayer::import(APLP1_args$gtf)
APLP2_gtf <- rtracklayer::import(APLP2_args$gtf)

# transcript information
APLP1_tx <- read.table(APLP1_args$tx, header = TRUE, sep="\t")
APLP2_tx <- read.table(APLP2_args$tx, header = TRUE, sep="\t")

# Functions ---------------------------------------------------------------

# Function to plot transcripts from GTF
plot_transcripts <- function(gtf, gene) {
  
  #gtf_gene_of_interest <- gtf[which(gtf$gene_id == gene), ]
  gtf_gene_of_interest <- gtf[startsWith(gtf$gene_id, gene), ]
  
  # Ensure gff is a data frame
  gff <- as.data.frame(gtf_gene_of_interest)
  
  # Extract exons and CDS
  exons <- gff %>% filter(type == "exon")
  cds <- gff %>% filter(type == "CDS")
  
  # Add UTR regions
  utr_data <- add_utr(exons = exons, cds = cds, group_var = "transcript_id")
  
  # Rescale exons and introns with shorten_gaps
  rescaled_data <- shorten_gaps(
    exons = utr_data,
    introns = to_intron(utr_data, "transcript_id"),
    group_var = "transcript_id"
  )
  
  # Plot
  transcript_plot <- rescaled_data %>%
    filter(type == "CDS") %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range() +
    geom_range(
      data = rescaled_data %>% filter(type == "UTR"),
      height = 0.25,
      fill = "white"
    ) +
    geom_intron(
      data = to_intron(
        rescaled_data %>% filter(type != "intron"),
        "transcript_id"
      ),
      arrow.min.intron.length = 100,
      strand = unique(as.character(gff$strand))
    ) +
    labs(
      y = "Transcript ID",
      x = "Genomic Position (Rescaled)",
      fill = "Region Type"
    ) +
    ggtitle(paste("Unique", gene, "transcripts")) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "top",
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold")
    )
  
  return(transcript_plot)
}

# Main --------------------------------------------------------------------

# Filter transcripts
APLP1_filtered <- APLP1_gtf[APLP1_gtf$transcript_id %in% APLP1_tx$isoform,]
APLP2_filtered <- APLP2_gtf[APLP2_gtf$transcript_id %in% APLP2_tx$isoform,]

# plot
APLP1_plot <- plot_transcripts(
  gtf = APLP1_filtered, 
  gene = "ENSG00000105290"
)

APLP2_plot <- plot_transcripts(
  gtf = APLP2_filtered, 
  gene = "ENSG00000084234"
)

