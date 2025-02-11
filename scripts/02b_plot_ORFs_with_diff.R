# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggtranscript)
library(patchwork)
library(rtracklayer)
library(plyranges)
library(here)
library(R.utils)
library(ggforce)

# Arguments---------------------------------------------------------------

Gene <- "APP"
MANE_transcript <- "ENST00000346798"

# Load data ---------------------------------------------------------------

Filtered_tx <- read.table(here::here("raw_data", "IsoSeq", "APP_classification_filtered.txt"), header = TRUE, sep="\t")
LR_gff <- rtracklayer::import(here::here("raw_data", "IsoSeq", "APP_corrected.gtf.cds.gff"))

# also download and load reference annotation 
ref_path <- here::here(tempdir(), "gencode.v38.annotation.gtf.gz")

if(!file.exists(ref_path)) {
  
  download.file(
    url = paste0(
      "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/", 
      "gencode.v38.annotation.gtf.gz"
    ),
    destfile = ref_path
  )
  
}

ref <- rtracklayer::import(ref_path)

# Functions ---------------------------------------------------------------

get_lr_tx_of_interest <- function(lr, pb_ids) {
  
  lr <- lr[lr$transcript_id %in% pb_ids]
  lr_cds <- lr[lr$type == "CDS"]
  
  return(list(cds = lr_cds))
  
}

get_mane <- function(ref, mane_id) {
  
  mane <- ref[!is.na(ref$transcript_id)] 
  GenomicRanges::mcols(mane)[["transcript_id"]] <- 
    GenomicRanges::mcols(mane)[["transcript_id"]] %>% 
    stringr::str_remove("\\..*")
  
  mane <- mane[mane$transcript_id == mane_id, ]
  mane_cds <- mane[mane$type == "CDS"]
  
  return(list(cds = mane_cds))
  
}

plot_transcripts <- function(lr_cds, mane_cds, MANE_canonical = "MANE") {
  
  transcript_order <- c(
    ORFs$isoform,
    mane_cds$cds$transcript_id %>% unique()
  )
  
  lr_mane_cds_df <- c(lr_cds$cds, mane_cds$cds) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  
  # Rescale exons and introns with shorten_gaps
  rescaled_data <- shorten_gaps(
    exons = lr_mane_cds_df,
    introns = to_intron(lr_mane_cds_df, "transcript_id"),
    group_var = "transcript_id"
  )
  
  # Ensure introns are correctly defined for visualization
  rescaled_introns <- to_intron(rescaled_data %>% filter(type != "intron"), "transcript_id")
  
  # Add to_diff() to visualize differences compared to MANE
  rescaled_diffs <- to_diff(
    exons = rescaled_data %>% 
      dplyr::filter(
        transcript_id != MANE_transcript,
        type == "CDS"),
    ref_exons = rescaled_data %>% dplyr::filter(transcript_id == MANE_transcript,
                                                type == "CDS"),
    group_var = "transcript_id"
  )
  
  transcript_plot <- rescaled_data %>% 
    dplyr::filter(type == "CDS") %>%
    ggplot(aes(
      xstart = start, 
      xend = end, 
      y = transcript_id)
    ) + 
    geom_range(height = 0.75, show.legend = FALSE) + 
    geom_intron(
      data = rescaled_introns,
      aes(strand = strand), 
      arrow.min.intron.length = 100
    ) + 
    geom_range(
      data = rescaled_diffs,
      aes(
        fill = diff_type, 
        colour = diff_type
      ),
      height = 0.75,
      alpha = 0.2, 
      linetype = 2
    ) +
    scale_y_discrete(name = "Transcript ID") + 
    scale_x_continuous(name = "Genomic position (Rescaled)") + 
    scale_fill_discrete(
      name = paste0("Coding DNA sequence (CDS) in ", "MANE ", "(", MANE_transcript, ")", " transcript:"), 
      labels = c(
        paste0("In ", MANE_transcript), 
        paste0("Not in ", MANE_transcript)
      )
    ) + 
    scale_colour_discrete(
      name = paste0("Coding DNA sequence (CDS) in ", "MANE ", "(", MANE_transcript, ")", " transcript:"), 
      labels = c(
        paste0("In ", MANE_transcript), 
        paste0("Not in ", MANE_transcript)
      )
    ) + 
    ggtitle(label = paste0(Gene, " coding transcripts")) +
    theme_bw() + 
    theme(
      legend.position = "top", 
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  return(transcript_plot)
  
}

# Main --------------------------------------------------------------------

# Get all unique ORFs
ORFs <- Filtered_tx %>% 
  dplyr::select(isoform, ORF_seq, ORF_length) %>% 
  na.omit() %>% 
  dplyr::distinct(ORF_seq, .keep_all = TRUE) %>% 
  dplyr::arrange(ORF_length)

lr_cds <- get_lr_tx_of_interest(
  lr = LR_gff, 
  pb_ids = ORFs$isoform)

mane_cds <- get_mane(ref = ref, mane_id = MANE_transcript)

transcript_plot <- 
  plot_transcripts(
    lr_cds = lr_cds, 
    mane_cds = mane_cds
  )

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcript_plot, 
  filename = "02b_APP_transcript_CDS_plot.svg", 
  path = here::here("results", "plots"), 
  width = 18, 
  height = 12, 
  dpi = 600, 
  bg = "white"
)
