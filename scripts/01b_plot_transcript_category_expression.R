# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggpubr)
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

plot_transcripts_per_gene <-
  function(data, gene, gene_name, samples, labelling) {
    
    # fill colour to use
    fill_colour <- c("Coding known (complete match)" = "#045a8d",
                     "Coding known (alternate 3/5 end)" = "#74add1",
                     "Coding novel" = "#4d9221",
                     "NMD novel" = "#d53e4f",
                     "Non-coding known" = "#b2abd2",
                     "Non-coding novel" = "#d8daeb")
    
    
    ## Plot number of transcripts per category ##
    Transcripts_per_category <-
      data %>%
      dplyr::select(isoform,
                    Isoform_class, 
                    starts_with("NFLR.")) %>%
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>%
      pivot_longer(!c(Isoform_class, isoform), 
                   names_to = "sample", 
                   values_to = "count") %>%
      dplyr::left_join(., 
                       dplyr::select(samples, 
                                     sample_id_matching, 
                                     Cell_category), 
                       by = c("sample" = "sample_id_matching")) %>% 
      subset(count != 0) %>% 
      dplyr::select(Cell_category, Isoform_class, isoform) %>% 
      unique() %>% 
      dplyr::count(Cell_category, Isoform_class) %>% 
      dplyr::mutate(source = ifelse(Cell_category %in% c("Total brain", "Cerebral cortex"), "Brain", "Cells"))
    
    
    Transcripts_per_category_plot <-
      Transcripts_per_category %>% 
      ggplot(aes(x = factor(Cell_category, 
                            levels = c("iPSC", 
                                       "Neuroepithelial", 
                                       "Neural progenitor", 
                                       "Neuron", 
                                       "Astrocyte (untreated)", 
                                       "Astrocyte (treated)",
                                       "Microglia (untreated)",
                                       "Microglia (treated)",
                                       "Total brain",
                                       "Cerebral cortex")), 
                 y = n, 
                 fill = Isoform_class)) +
      geom_col(show.legend = F, colour = "Black", width = 0.7) +
      scale_fill_manual(values = fill_colour) +
      scale_x_discrete(labels = c("Coding known (alternate 3/5 end)" = "Coding known\n(alternate 3/5 end)",
                                  "Coding known (complete match)" = "Coding known\n(complete match)")) +
      labs(x = "", y = "No. unique transcripts") +
      ggtitle(paste0("Number of unique ", gene_name, " transcripts")) +
      facet_grid(~factor(source, levels = c("Cells", "Brain")), scale = "free_x", space = "free") +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_text(size = 14),
            axis.title = element_text(size = 16), 
            strip.text = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(face = "bold", size = 12),
            axis.title.y = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
    
    
    ## Plot expression per category ##
    Expression_per_category <-
      data %>%
      dplyr::select(Isoform_class, 
                    associated_gene, 
                    starts_with("NFLR.")) %>%
      rename_with(~paste0("", gsub("NFLR.Clontech_5p..|_3p", "", .)),
                  starts_with('NFLR.Clontech_5p..')) %>%
      pivot_longer(!c(Isoform_class, associated_gene), 
                   names_to = "sample", 
                   values_to = "count") %>% 
      dplyr::left_join(., 
                       dplyr::select(samples, 
                                     sample_id_matching, 
                                     Cell_category), 
                       by = c("sample" = "sample_id_matching")) %>% 
      aggregate(count ~ Isoform_class + associated_gene + Cell_category + sample,
                data = .,
                FUN = "sum")
    
    Expression_per_category_plot <-
      Expression_per_category %>%
      group_by(Cell_category, Isoform_class) %>%
      summarise(mean_count = mean(count), sd_count = sd(count)) %>% 
      dplyr::mutate(source = ifelse(Cell_category %in% c("Total brain", "Cerebral cortex"), "Brain", "Cells")) %>% 
      ggplot(
        aes(x = factor(Cell_category, 
                       levels = c("iPSC", 
                                  "Neuroepithelial", 
                                  "Neural progenitor", 
                                  "Neuron", 
                                  "Astrocyte (untreated)", 
                                  "Astrocyte (treated)",
                                  "Microglia (untreated)",
                                  "Microglia (treated)",
                                  "Total brain",
                                  "Cerebral cortex")), 
            y = mean_count, 
            fill = Isoform_class)) +
      geom_col(position = "dodge", color = "black") +
      geom_errorbar(aes(ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                    position = position_dodge(0.9), width = 0.2) +
      labs(x = "", 
           y = "Expression per transcript category") +
      ggtitle(paste0(gene_name, " transcript expression by sample type")) +
      scale_fill_manual(name = "Transcript category",
                        values = fill_colour) +
      scale_y_continuous(labels = function(x) paste0(x, '%')) +
      facet_grid(~factor(source, levels = c("Cells", "Brain")), scale = "free_x", space = "free") +
      theme_bw() +
      theme(legend.position = "top",
            legend.title = element_text(size = 14),
            axis.title = element_text(size = 16),
            strip.text = element_text(face = "bold", size = 12),
            axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(face = "bold", size = 12),
            axis.title.y = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
    
    plot <- ggpubr::ggarrange(Transcripts_per_category_plot, Expression_per_category_plot, 
                              nrow = 1, 
                              common.legend = T, align = "v",
                              labels = labelling, 
                              font.label = list(size = 24),
                              widths = c(1, 1.5))
    
    return(plot)
    
    }

# Main --------------------------------------------------------------------

transcript_plot <-
  plot_transcripts_per_gene(data = Transcripts, 
                            gene = "ENSG00000142192.21", 
                            gene_name = Gene, 
                            samples = Samples,
                            labelling = c("B", "C"))


# Save data ---------------------------------------------------------------

file_extensions <- c("png", "svg")

for (ext in file_extensions) {
  ggsave(
    plot = transcript_plot,
    filename = paste0("01b_transcript_category_plot.", ext),
    path = here::here("results", "plots"),
    width = 16,
    height = 6,
    dpi = 600,
    bg = "white"
  )
}
