# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggtranscript)
library(patchwork)
library(rtracklayer)
library(plyranges)
library(here)
library(R.utils)

# Arguments ---------------------------------------------------------------

Gene <- "APP"

# Load data ---------------------------------------------------------------

Filtered_tx <- read.table(here::here("raw_data", "IsoSeq", "APP_classification_filtered.txt"), header = TRUE, sep="\t")
LR_gff <- rtracklayer::import(here::here("raw_data", "IsoSeq", "APP_corrected.gtf.cds.gff"))

# Functions ---------------------------------------------------------------

get_TOI <- function(gff, pb_ids) {
  gff <- gff[gff$transcript_id %in% pb_ids]
  gff_exons <- gff[gff$type == "exon"]
  gff_cds <- gff[gff$type == "CDS"]
  
  list(
    exons = data.frame(gff_exons),
    cds = data.frame(gff_cds)
  )
}

annotate_TX <- function(tx, gene) {
  sapply(tx, function(x) {
    tmp <- Filtered_tx
    x$transcript_biotype <- tmp[match(x$transcript_id, tmp$isoform),]$Isoform_class
    x
  }, simplify = FALSE)
}

plot_for_class <- function(gene_name, isoform_class, TOI) {
  TOI <- TOI %>%
    lapply(function(x) dplyr::filter(x, transcript_biotype == isoform_class))
  
  exons <- TOI$exons
  introns <- exons %>% to_intron(group_var = "transcript_id")
  cds <- TOI$cds
  
  exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(height = 0.25) +
    geom_range(data = cds) +
    geom_intron(
      data = introns,
      arrow.min.intron.length = 1500,
      strand = unique(as.character(exons$strand))
    ) +
    labs(
      y = "Transcript ID",
      x = "Genomic position (hg38)"
    ) +
    ggtitle(paste(gene_name, "-", isoform_class, "transcripts")) +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"))
}

# Main --------------------------------------------------------------------

# Get list of unique Isoform_class values
isoform_classes <- unique(Filtered_tx$Isoform_class)

# Generate TOI data for all transcripts
TOI <- get_TOI(
  gff = LR_gff,
  pb_ids = Filtered_tx$isoform
) %>%
  annotate_TX(tx = ., gene = Gene)

# Create a list of plots, one for each Isoform_class
plots <- lapply(isoform_classes, function(class) {
  plot_for_class(Gene, class, TOI)
})

# Save plots to files with dynamic height
output_dir <- here::here("results", "plots")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save plots with height proportional to the number of unique transcript IDs
lapply(seq_along(plots), function(i) {
  # Get the number of unique transcript IDs for the current isoform class
  num_transcripts <- TOI$exons %>%
    filter(transcript_biotype == isoform_classes[i]) %>%
    pull(transcript_id) %>%
    unique() %>%
    length()
  
  # Calculate height based on the number of transcripts (e.g., 0.2 units per transcript)
  plot_height <- max(5, num_transcripts * 0.2) # Minimum height is 5
  
  # Clean the isoform class name to remove spaces and special characters
  safe_isoform_class <- gsub("[ /]", "_", isoform_classes[i])
  
  # Save the plot
  ggsave(
    filename = paste0(output_dir, "/", "02a", "_", safe_isoform_class, ".png"),
    plot = plots[[i]],
    width = 8,
    height = plot_height,
    dpi = 600,
    bg = "white"
  )
})
