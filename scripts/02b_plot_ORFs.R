# Load libraries ----------------------------------------------------------

library(tidyverse)
library(ggtranscript)
library(Biostrings)
library(ggforce)
library(rtracklayer)
library(plyranges)
library(here)
library(R.utils)
library(ggforce)

# Arguments---------------------------------------------------------------

Gene <- "APP"
MANE_transcript <- "ENST00000346798"

# Load data ---------------------------------------------------------------

tx_classification <- read.table(here::here("raw_data", "IsoSeq", "APP_classification_filtered.txt"), header = TRUE, sep="\t")
tx_gff <- rtracklayer::import(here::here("raw_data", "IsoSeq", "APP_corrected.gtf.cds.gff")) 

# also get the CCDS data from NCBI
ccds_path <- here::here(tempdir(), "CCDS_nucleotide.current.fna.gz")

if (!file.exists(ccds_path)) {
  
  download.file(
    url = "https://ftp.ncbi.nih.gov/pub/CCDS/current_human/CCDS_nucleotide.current.fna.gz",
    destfile = ccds_path
  )
  
}

R.utils::gunzip(ccds_path, remove = TRUE)

# Read CCDS FASTA file
CCDS <- Biostrings::readDNAStringSet(stringr::str_remove(ccds_path, "\\.gz"))

# Functions ---------------------------------------------------------------

# This function annotates Tx data with CCDS information and outputs the data.frame and GRanges as a list
annotate_tx_data <- function(tx, col_W_seqs, CCDS_seqs, gff) {
  
  # Ensure the column name exists in tx
  if (!col_W_seqs %in% colnames(tx)) {
    stop("Column name not found in tx")
  }
  
  # Match sequences
  matches <- match(tx[[col_W_seqs]], CCDS_seqs)
  
  # Bind results and get corresponding CCDS names
  tx <- tx %>%
    mutate(
      ccdsid = ifelse(
        !is.na(matches), 
        names(CCDS_seqs)[matches], 
        paste0("Novel CDS (", nchar(gsub("[^A-Za-z]", "", ORFik_aa)), "aa)") # this only allows for alpabetical as ORFik adds "*" to reperesent Ter
      )
    )
  
  # Get transcript and CCDS IDs
  tx_CCDS <- tx %>% dplyr::select(isoform, ccdsid)
  
  # Convert gff to a DataFrame to enable manipulation
  gff_df <- as.data.frame(mcols(gff))
  
  # Ensure transcript_id column exists before merging
  if (!"transcript_id" %in% colnames(gff_df)) {
    stop("Column 'transcript_id' not found in gff metadata")
  }
  
  # Merge with tx_CCDS to add ccdsid where transcript_id matches isoform
  gff_df <- gff_df %>% 
    dplyr::left_join(tx_CCDS, by = c("transcript_id" = "isoform"))
  
  # Assign the updated metadata back to the GRanges object
  mcols(gff) <- gff_df
  
  return(list(tx_annotated = tx, gff_annotated = gff))
}


# This function gets the Tx of interest from the GRanges object
# It only includes CDS exons
get_lr_tx_of_interest <- function(lr, pb_ids, exons) {
  
  lr <- lr[lr$transcript_id %in% pb_ids]
  lr_exons <- lr[lr$type %in% exons]
  
  return(lr_exons)
  
}


plot_transcripts <- function(lr_cds) {
  
  transcript_order <- ORFs$isoform
  
  lr_cds_df <- lr_cds %>% 
    as.data.frame() %>% 
    dplyr::filter(transcript_id %in% ORFs$isoform) %>% 
    dplyr::mutate(
      transcript_id = transcript_id %>% factor(levels = transcript_order)
    )
  
  # Rescale exons and introns with shorten_gaps
  rescaled_data <- shorten_gaps(
    exons = lr_cds_df,
    introns = to_intron(lr_cds_df, "transcript_id"),
    group_var = "transcript_id"
  ) %>% 
    dplyr::mutate(
      Novelty = ifelse(startsWith(ccdsid, "CCDS"), "In annotation", "Not in annotation")
    )
  
  
  # Plot trasnscripts
  transcript_plot <- rescaled_data %>% 
    dplyr::filter(type == "CDS") %>%
    ggplot(aes(
      xstart = start, 
      xend = end, 
      y = ccdsid)
    ) + 
    geom_range(height = 0.75, show.legend = FALSE) + 
    geom_intron(
      data = rescaled_data %>% dplyr::filter(type == "intron"),
      aes(strand = strand), 
      arrow.min.intron.length = 100
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
    facet_col(vars(Novelty), scale = "free_y", space = "free") +
    ggtitle(label = paste0(Gene, " open reading frames")) +
    theme_bw() + 
    theme(
      legend.position = "top", 
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(transcript_plot)
  
}

# Main --------------------------------------------------------------------

# Annotate with CCDS data
tx_data_CCDS <- annotate_tx_data(tx = tx_classification, 
                                 col_W_seqs = "ORFik_cds",
                                 CCDS_seqs = CCDS, 
                                 gff = tx_gff)
  


# Get all unique ORFs
ORFs <- tx_data_CCDS$tx_annotated %>%
  
  # This part is to only include ORFs that are concordant between the two callers
  filter(!is.na(ORF_seq) & !is.na(ORFik_aa)) %>% # Remove NAs
  mutate(ORFik_aa = sub("\\*$", "", ORFik_aa)) %>% # Remove trailing "*"
  filter(ORF_seq == ORFik_aa) %>%  # Keep only rows where values match
  
  dplyr::select(isoform, ORF_seq, ORF_length) %>% 
  na.omit() %>% 
  dplyr::distinct(ORF_seq, .keep_all = TRUE) %>% 
  dplyr::arrange(ORF_length)

# get tx with those ORFs
lr_cds <- get_lr_tx_of_interest(
  lr = tx_data_CCDS$gff_annotated, 
  pb_ids = ORFs$isoform, 
  exons = "CDS")

transcript_plot <- 
  plot_transcripts(
    lr_cds = lr_cds
  )

# Save data ---------------------------------------------------------------

ggsave(
  plot = transcript_plot, 
  filename = "02b_APP_transcript_CDS_plot.svg", 
  path = here::here("results", "plots"), 
  width = 14, 
  height = 8, 
  dpi = 600, 
  bg = "white"
)
