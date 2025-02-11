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
  lr_exons <- lr[lr$type == "exon"]
  lr_cds <- lr[lr$type == "CDS"]
  
  # return as list as we need both exons and cds
  lr_exons_cds <- list(
    exons = lr_exons, 
    cds = lr_cds
  )
  
  return(lr_exons_cds)
  
}

get_mane <- function(ref, mane_id) {
  
  # remove any NA transcript ids (i.e. type == "gene")
  mane <- ref[!is.na(ref$transcript_id)] 
  
  # remove to .XX after ENST
  GenomicRanges::mcols(mane)[["transcript_id"]] <- 
    GenomicRanges::mcols(mane)[["transcript_id"]] %>% 
    stringr::str_remove("\\..*")
  
  mane <- mane[mane$transcript_id == mane_id, ]
  mane_exons <- mane[mane$type == "exon"]
  mane_cds <- mane[mane$type == "CDS"]
  
  mane_exons_cds <- list(
    exons = mane_exons, 
    cds = mane_cds
  )
  
  return(mane_exons_cds)
  
}

plot_diff <- function(lr_exons_cds, 
                      mane_exons_cds, 
                      lr_mane_diffs, 
                      MANE_canonical = "MANE"
) {
  
  # merge mane and lr data and convert to data.frame() for plotting
  # convert transcript_id to factor to make sure mane is at top
  transcript_order <- c(
    lr_exons_cds$exons$transcript_id %>% unique(),
    mane_exons_cds$exons$transcript_id %>% unique()
  )
  
  lr_mane_exons_df <- c(lr_exons_cds$exons, mane_exons_cds$exons) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  lr_mane_cds_df <- c(lr_exons_cds$cds, mane_exons_cds$cds) %>% 
    as.data.frame() %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% 
        factor(levels = transcript_order)
    )
  
  # plot diff plot
  diff_plot <- lr_mane_exons_df %>% 
    ggplot(aes(
      xstart = start, 
      xend = end, 
      y = factor(transcript_id, levels = c("PB.4.188",
                                           "PB.4.510",
                                           "PB.4.998",
                                           "PB.4.251",
                                           "PB.4.557",
                                           "ENST00000346798"))
    )) + 
    geom_range(
      height = 0.25,
      fill = "white"
    ) +
    geom_range(
      data = lr_mane_cds_df
    ) + 
    geom_intron(
      data = to_intron(lr_mane_exons_df, "transcript_id"), 
      aes(strand = strand), 
      arrow.min.intron.length = 400
    ) + 
    geom_range(
      data = lr_mane_diffs, 
      aes(
        fill = diff_type, 
        colour = diff_type
      ), 
      alpha = 0.2, 
      linetype = 2
    ) + 
    scale_y_discrete(name = "Transcript ID") + 
    scale_x_continuous(name = "Genomic position (hg38)") + 
    scale_fill_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    scale_colour_discrete(
      name = paste0("Region in ", MANE_canonical, " transcript:"), 
      labels = c(
        paste0("In ", MANE_canonical), 
        paste0("Not in ", MANE_canonical)
      )
    ) + 
    facet_zoom(xlim = c(25895000, 25915000)) +
    ggtitle(label = paste0(Gene, " novel coding transcripts compared to MANE select")) +
    theme_bw() + 
    theme(legend.position = "top", plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  
  return(diff_plot)
  
}

# Main --------------------------------------------------------------------

lr_exons_cds <- get_lr_tx_of_interest(
  lr = LR_gff, 
  pb_ids = c("PB.4.557", "PB.4.251", "PB.4.998", "PB.4.510", "PB.4.188"))

mane_exons_cds <- get_mane(ref = ref, mane_id = MANE_transcript)

lr_mane_diffs <- 
  ggtranscript::to_diff(
    exons = lr_exons_cds$exons %>% as.data.frame(),
    ref_exons = mane_exons_cds$exons %>% as.data.frame(), 
    group_var = "transcript_id"
  )

lr_mane_diff_plot <- 
  plot_diff(
    lr_exons_cds = lr_exons_cds, 
    mane_exons_cds = mane_exons_cds, 
    lr_mane_diffs = lr_mane_diffs
  )

# Save data ---------------------------------------------------------------

ggsave(
  plot = lr_mane_diff_plot, 
  filename = "02c_APP_transcripts_diff_plot.svg", 
  path = here::here("results", "plots"), 
  width = 10, 
  height = 6, 
  dpi = 600, 
  bg = "white"
)
