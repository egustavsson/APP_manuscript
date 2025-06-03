# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(plyranges)
library(ggtranscript)
library(cowplot)

# Arguments ---------------------------------------------------------------

args <-
  list(
    snRNA_files = list.files(here::here("raw_data", "snRNA_5prime"), 
                             full.names = TRUE, 
                             pattern = glob2rx("PFCTX*rev*")),
    CUTRUN_files = list.files(here::here("raw_data", "snRNA_5prime"), 
                              full.names = TRUE, 
                              pattern = "^CR116")
  )

# Load data ---------------------------------------------------------------

# Long-read transcripts
long_read <- rtracklayer::import(
  here::here("raw_data", "IsoSeq", "APP_corrected.gtf.cds.gff")) %>% 
  plyranges::filter(transcript_id %in% c(
    "PB.4.557",
    "PB.4.251"
  ))


# Functions ---------------------------------------------------------------

plot_suppotive_data_per_gene <-
  function(seqnames, start, end, strand, gene_id) {
    
    # Filter long read gtf/gff for gene of interest
    long_read <- long_read[long_read$gene_id == gene_id]
    
    
    # loci used to filter data
    locus_subset <- 
      GRanges(seqnames = seqnames,
              ranges = IRanges(start = start, 
                               end = end),
              strand = strand)
    
    # 5' single nuclear RNAseq data
    snRNA_data <- lapply(args$snRNA_files, function(x) {
      
      rtracklayer::import.bw(x) %>%
        subsetByOverlaps(., locus_subset) %>% 
        plyranges::mutate(
          sample = stringi::stri_replace_all_regex(x,
                                                   paste0(here::here("data", "5prime_CUTRUN"), "/"), 
                                                   "", 
                                                   vectorize_all = F),
          cell_type = case_when(grepl("Excitatoryneurons", sample) ~ "Excitatory neurons",
                                grepl("Inhibitoryneurons", sample) ~ "Inhibitory neurons",
                                grepl("Astrocytes", sample) ~ "Astrocytes",
                                grepl("Microglia", sample) ~ "Microglia",
                                grepl("OPC", sample) ~ "Oligodendrocyte progenitor cells",
                                grepl("Oligodendrocytes", sample) ~ "Oligodendrocytes")) %>% 
        data.frame()
      
    }
    ) %>% 
      bind_rows(!!!.)
    
    # FANTOM5 CAGE seq data
    CAGE_seq <- 
      rtracklayer::import.bed(here::here("raw_data", "hg38.cage_peak_phase1and2combined_coord.bed.gz")) %>% 
      subsetByOverlaps(., locus_subset)
    
    # Plot GBA and GBAP1 transcripts
    exons <- data.frame(long_read) %>% dplyr::filter(type == "exon")
    introns <- exons %>% to_intron(group_var = "transcript_id")
    CDS <- data.frame(long_read) %>% dplyr::filter(type == "CDS")
    
    APP_plot <-
      exons %>%
      ggplot(
        aes(xstart = start, xend = end, y = factor(transcript_id, levels = c("PB.4.557",
                                                                             "PB.4.251"
                                                                             )))) +
      geom_range(fill = "white",
                 height = 0.25) +
      geom_range(data = CDS) +
      geom_intron(data = introns, 
                  arrow.min.intron.length = 500, 
                  arrow = grid::arrow(ends = "first", length = grid::unit(0.1, "inches"))) +
      labs(y = "Transcript name",
           x = "") +
      xlim(start(locus_subset), end(locus_subset)) +
      theme_classic() +
      theme(axis.title = element_text(face = "bold", size = 14),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 10),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    
    # Plot CAGE seq coords (hg38)
    CAGE_plot <-
      CAGE_seq %>% 
      data.frame() %>% 
      ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = 1)) +
      geom_rect(colour = "black", fill = "black", show.legend = F, alpha=0.8) +
      xlim(start(locus_subset), end(locus_subset)) +
      labs(y = "CAGE peaks",
           x = "") +
      theme_classic() +
      theme(axis.title = element_text(size = 12),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 10,
                                       colour = "white"),
            axis.ticks.y = element_line(colour = "white"),
            strip.text.y = element_text(face = "bold",
                                        size = 10),
            strip.background = element_rect(fill ="lightgrey"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.spacing = unit(0, "cm"))
    
    # Plot 5' single nuclear RNAseq
    snRNA_plot <-
      snRNA_data %>% 
      ggplot(aes(xmin = start, xmax = end, ymin = 0, ymax = score)) +
      geom_rect(colour = "lightsalmon", fill = "lightsalmon", show.legend = F, alpha=0.8) +
      xlim(start(locus_subset), end(locus_subset)) +
      labs(y = "Expression levels (RPKM)",
           x = "Genomic coordinate (hg38)") +
      facet_wrap(vars(cell_type), 
                 ncol = 1, 
                 strip.position = "right",
                 labeller = label_wrap_gen(18)) +
      theme_classic() +
      theme(axis.title = element_text(face = "bold", size = 14),
            axis.text.x = element_text(size = 14),
            axis.ticks.x = element_line(),
            axis.text.y = element_text(size = 10),
            strip.text.y = element_text(face = "bold",
                                        size = 10),
            strip.background = element_rect(fill ="lightgrey"),
            panel.border = element_rect(colour = "black", fill=NA, size=1),
            panel.spacing = unit(0, "cm"))
    
    # Final plot
    APP_transcript_supportive_data_plot <-
      plot_grid(APP_plot, 
                CAGE_plot,
                snRNA_plot, 
                ncol = 1, 
                align = "hv", 
                rel_heights = c(1, 0.5, 3.2), 
                axis = "lr",
                label_size = 18)
    
    return(APP_transcript_supportive_data_plot)
  }


# Main --------------------------------------------------------------------

# Plot APP
APP_transcript_supportive_data_plot <-
  plot_suppotive_data_per_gene(seqnames = "chr21", start = 25879000, end = 25899500, strand = "-", gene_id = "ENSG00000142192.21")

# Save data ---------------------------------------------------------------

ggsave(plot = APP_transcript_supportive_data_plot, 
       filename = "03a_APP_transcript_supportive_data_plot.svg", 
       path = here::here("results", "plots"), 
       width = 10, 
       height = 14, 
       dpi = 600
)
