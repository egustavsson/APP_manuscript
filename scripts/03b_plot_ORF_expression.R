# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
library(rtracklayer)
library(readxl)

# Arguments---------------------------------------------------------------

Gene <- "APP"

args <-
  list(
    path_to_transcripts = here::here("raw_data", "IsoSeq", paste0(Gene, "_classification_filtered.txt")),
    path_to_samples = here::here("raw_data", "IsoSeq", "Charlies_samples.xlsx")
  )

# Load data ---------------------------------------------------------------

Transcripts <-
  read.table(args$path_to_transcripts, header = TRUE, sep="\t")

Samples <-
  read_excel(args$path_to_samples) %>% 
  dplyr::mutate(sample_id_matching = gsub(Barcode_ID, pattern = c("_3p"), replacement = "")) %>% 
  dplyr::mutate(Cell_category = case_when(Cell_type == "iPSC" ~ "iPSC",
                                          Cell_type == "Neuroepithelial" ~ "Neuroepithelial",
                                          Cell_type == "Neural progenitors" ~ "Neural progenitor",
                                          Cell_type == "Neuron" ~ "Neuron",
                                          Cell_type == "Astrocyte" & Treatment == "None" ~ "Astrocyte (untreated)",
                                          Cell_type == "Astrocyte" & Treatment != "None" ~ "Astrocyte (treated)",
                                          Cell_type == "Microglia" & Treatment == "None" ~ "Microglia (untreated)",
                                          Cell_type == "Microglia" & Treatment != "None" ~ "Microglia (treated)",
                                          Sample == "Total brain RNA" ~ "Total brain",
                                          Sample == "Human brain cerebral cortex polyA+RNA" ~ "Cerebral cortex")) %>% 
  drop_na(Cell_category) %>% 
  data.frame()


# Functions ---------------------------------------------------------------

plot_ORF_expression <- function(data, ORFs, Samples) {
  
  # Define all relevant sample categories
  Sample_categories <- c(
    "iPSC", "Neuroepithelial", "Neural progenitor", 
    "Neuron", "Astrocyte (untreated)", "Astrocyte (treated)", 
    "Microglia (untreated)", "Microglia (treated)",
    "Cerebral cortex", "Total brain"
  )
  
  # Process data for ORF expression
  Expression_per_ORF <- 
    data %>%
    dplyr::filter(ORF_seq %in% ORFs) %>% 
    dplyr::select(ORF_length, ORF_seq, starts_with("NFLR.")) %>% 
    dplyr::mutate(across(starts_with("NFLR."), ~(. / sum(.)) * 100)) %>%
    dplyr::rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)), dplyr::starts_with('NFLR.Clontech_5p..')) %>% 
    tidyr::pivot_longer(!c(ORF_seq, ORF_length), names_to = "Sample", values_to = "Relative_expression") %>%
    dplyr::left_join(., dplyr::select(Samples, 
                                      sample_id_matching, 
                                      Cell_category), 
                     by = c("Sample" = "sample_id_matching")) %>% 
    dplyr::mutate(Cell_category = factor(Cell_category, levels = Sample_categories)) %>% 
    dplyr::group_by(ORF_length, Cell_category, ORF_seq) %>% 
    dplyr::mutate(ORF_length = paste0(ORF_length, "aa"),
                  source = ifelse(Cell_category %in% c("Total brain", "Cerebral cortex"), 
                                  "CNS", "Cells")) %>% 
    dplyr::select(ORF_length, Cell_category, Relative_expression, source)
  
  # Plot
  Expression_per_ORF_plot <- 
    Expression_per_ORF %>%
    dplyr::mutate(color_code = ifelse(ORF_length %in% c("100aa", "191aa"), "Novel ORF", "Known ORF")) %>% 
    ggplot(aes(x = Cell_category, 
               y = Relative_expression,
               fill = color_code)) +
    geom_boxplot() +
    labs(title = "APP open reading frame expression in Cells and in CNS",
         x = "", 
         y = "Relative expression by ORF (%)",
         fill = "ORF annotation") +  # Set legend title
    scale_y_continuous(labels = function(x) paste0(x, '%')) +
    scale_x_discrete() +  # Ensure x-axis follows the factor levels
    scale_fill_manual(
      values = c("Novel ORF" = "green4", "Known ORF" = "lightgrey")
    ) +
    facet_grid(cols = vars(
      factor(source, levels = c("Cells", "CNS"))
    ), 
    rows = vars(ORF_length), 
    scales = "free", 
    space = "free_x") +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "top",
      legend.title = element_text(size = 16, face = "bold"),  # Bold legend title
      legend.text = element_text(size = 14, face = "bold"),   # Bold legend text
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 18), 
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face = "bold"),
      axis.text.y = element_text(size = 14),
      plot.margin = margin(t = 10, r = 10, b = 20, l = 22)
    )
  
  return(Expression_per_ORF_plot)
}


# Main --------------------------------------------------------------------

# Get all the ORFs of interest from iPSC data
ORF_of_interest <- Transcripts %>% 
  dplyr::filter(Isoform_class == "Coding known (alternate 3/5 end)") %>% 
  dplyr::select(ORF_length, ORF_seq) %>% 
  na.omit() %>% 
  unique() %>% 
  pull(ORF_seq)


# Plot expression by ORF
expression_by_ORF_plot <- plot_ORF_expression(data = Transcripts, 
                                              ORFs = ORF_of_interest,
                                              Samples = Samples)


# Save data ---------------------------------------------------------------

ggsave(
  plot = expression_by_ORF_plot, 
  filename = "03b_APP_expression_by_ORF_plot.svg", 
  path = here::here("results", "plots"), 
  width = 9, 
  height = 14, 
  dpi = 600, 
  bg = "white"
)
