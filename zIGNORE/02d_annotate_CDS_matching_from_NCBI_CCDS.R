# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(Biostrings)
library(rtracklayer)

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
      CCDS_id = ifelse(
        !is.na(matches), 
        names(CCDS_seqs)[matches], 
        paste0("Novel CDS (", nchar(gsub("[^A-Za-z]", "", ORFik_aa)), "aa)") # this only allows for alpabetical as ORFik adds "*" to reperesent Ter
      )
    )
  
  # Get transcript and CCDS IDs
  tx_CCDS <- tx %>% dplyr::select(isoform, CCDS_id)
  
  # Convert gff to a DataFrame to enable manipulation
  gff_df <- as.data.frame(mcols(gff))
  
  # Ensure transcript_id column exists before merging
  if (!"transcript_id" %in% colnames(gff_df)) {
    stop("Column 'transcript_id' not found in gff metadata")
  }
  
  # Merge with tx_CCDS to add CCDS_id where transcript_id matches isoform
  gff_df <- gff_df %>% 
    dplyr::left_join(tx_CCDS, by = c("transcript_id" = "isoform"))
  
  # Assign the updated metadata back to the GRanges object
  mcols(gff) <- gff_df
  
  return(list(tx_annotated = tx, gff_annotated = gff))
}


# Main --------------------------------------------------------------------

tx_data_CCDS <- annotate_tx_data(tx = tx_classification, 
                                 col_W_seqs = "ORFik_cds", 
                                 CCDS_seqs = CCDS, 
                                 gff = tx_gff)

# Save data ---------------------------------------------------------------



