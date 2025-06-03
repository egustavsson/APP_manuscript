# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)

# Load data ---------------------------------------------------------------

qPCR_data <-
  read_tsv(
    here::here("results", "qPCR", "qPCR_reformatted.txt"), show_col_types = F)

# Functions ---------------------------------------------------------------

plot_qPCR <- function(data) {
    
    plot <- 
      data %>%
      ggplot(aes(x = factor(Genotype, levels = c("WT", "KO")), y = `Relative expression`, fill = Genotype)) +
      geom_boxplot(show.legend = FALSE, 
                   outlier.shape = NA,
                   na.rm = TRUE,
                   width = 0.7) +
      geom_point(show.legend = FALSE,
                 na.rm = TRUE,
                 size = 2,
                 alpha = 0.6,
                 shape = 21) +
      facet_wrap(vars(factor(Isoform, levels = c("APP - full length", "APP - C100"))), scales = "free_y") +
      labs(
        title = "APP expression", 
        x = "",
        y = "Relative expression") +
      scale_fill_brewer(palette = "Blues") +
      scale_colour_brewer(palette = "Blues") +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
      )
  
  return(plot)
}

# Main --------------------------------------------------------------------

# Create plots including multiple treatments
plot <- plot_qPCR(
  data = qPCR_data
)


# Save data ---------------------------------------------------------------

ggsave(plot = plot, 
       filename = "05_qPCR_plot.svg", 
       path = here::here("results", "plots"), 
       width = 5, 
       height = 4, 
       dpi = 600, 
       bg = "white"
)
