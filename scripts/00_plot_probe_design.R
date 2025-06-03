library(tidyverse)
library(here)
library(ggtranscript)
library(GenomicRanges)
library(rtracklayer)
library(cowplot)

# Arguments ---------------------------------------------------------------

# The probes are hg19 coordinates so need to liftover. This is the chain file to use
chainFile <- here::here("raw_data", "hg19ToHg38.over.chain")

# Load data ---------------------------------------------------------------

# Get GENCODE reference
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

R.utils::gunzip(ref_path, remove = TRUE)

ref <- rtracklayer::import(stringr::str_remove(ref_path, "\\.gz"))

# The probe design
probes <-
  read.table(here::here("raw_data", "IDT_probe_design_final_Charlie.csv"), header = T, sep = ",") %>% 
  dplyr::select(!X) %>% 
  GRanges()

# Functions ---------------------------------------------------------------

plot_by_gene <- function(gene_name, reference, probes, seqnames, start, end, strand) {
  
  # loci used to filter data
  locus_subset <- 
    data.frame(start = start,
               end = end)
  
  # filter ref for gene of interest
  gene <- 
    reference %>% 
    data.frame() %>% 
    dplyr::filter(gene_name == !!gene_name,
                  type != "gene")
  
  
  # prepare ggtranscript input
  exons <- gene[gene$type == "exon", ]
  introns <- exons %>% to_intron(group_var = "transcript_id")
  cds <- gene[gene$type == "CDS", ]
  
  # transcript plot
  transcript_plot <-
    exons %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(height = 0.25) +
    geom_range(
      data = cds
    ) +
    geom_intron(
      data = introns,
      arrow.min.intron.length = 3500,
      strand = "-"
    ) +
    labs(
      y = "Transcript ID",
      x = ""
    ) +
    scale_fill_brewer(palette = "Dark2") +
    xlim(locus_subset$start, locus_subset$end) +
    theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  # plot probes
  # Create a new variable 'row' to specify the row for each rectangle
  probes$y <- runif(nrow(probes))
  
  
  probe_plot <-
    probes %>%
    ggplot(aes(xmin = start, xmax = end, ymin = y - 0.015, ymax = y + 0.015)) +
    geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8, position = position_nudge(x = 0.05)) +
    xlim(locus_subset$start, locus_subset$end) +
    labs(
      y = "Probes",
      x = "Genomic position (hg38)"
    ) +
    theme_bw() +
    theme(axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 12, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  final_plot <- plot_grid(transcript_plot,
                          probe_plot,
                          ncol = 1,
                          align = "hv", 
                          rel_heights = c(1, 0.5), 
                          axis = "lr",
                          label_size = 18)
  
  return(final_plot)
  
}

# Main --------------------------------------------------------------------

# Liftover of probe coordinates
chain <- rtracklayer::import.chain(chainFile)

probes_hg38 <- 
  rtracklayer::liftOver(probes, chain) %>% 
  unlist() %>% 
  data.frame()

transcript_plot <- plot_by_gene(gene_name = "APP", 
                                reference = ref, 
                                probes = probes_hg38, 
                                start = 25879071, 
                                end = 26173261)

# Save data -------------------------------------------------------------------------------------------

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = transcript_plot,
    filename = paste0("00a_probe_design_plot.", ext),
    path = here::here("results", "plots"),
    width = 12,
    height = 6,
    dpi = 600
  )
}
