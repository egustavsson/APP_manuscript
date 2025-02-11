# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(rstatix)
library(ggpubr)

# Load data ---------------------------------------------------------------

MSD_data <-
  read_tsv(
    here::here("results", "MSD_abeta", "abeta_MSD_reformatted.txt"), show_col_types = F)

# Functions ---------------------------------------------------------------

plot_abeta <- function(data, split_by, Treatment) {
  
  plot_order <- c("WT", "HET", "KO")
  
  # Filter data based on the specified Treatment(s)
  filtered_data <- data %>% 
    filter(Treatment %in% !!Treatment)
  
  # Split the data by the specified column
  data_list <- split(x = filtered_data, f = filtered_data[[split_by]])
  
  # Create a list to store the plots
  plot_list <- list()
  
  # Loop through each element in data_list to create a plot
  for (name in names(data_list)) {
    subset_data <- data_list[[name]]
    
    plot <- 
      subset_data %>%
      ggplot(aes(x = factor(Genotype, levels = plot_order), y = Count, fill = Genotype)) +
      geom_boxplot(show.legend = FALSE, 
                   outlier.shape = NA,
                   na.rm = TRUE,
                   width = 0.7) +
      geom_point(show.legend = FALSE,
                 na.rm = TRUE,
                 size = 2, 
                 alpha = 0.6, 
                 shape = 21) +
      facet_grid(Treatment ~ Analyte) +  # Facet by Treatment and Analyte
      labs(
        title = unique(subset_data$Measurement), 
        x = "Genotype",
        y = unique(subset_data$unit)) +
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
    
    # Store the plot in the list with the name of the subset
    plot_list[[name]] <- plot
  }
  
  # Return the list of plots
  return(plot_list)
}

# Main --------------------------------------------------------------------

# Create plots including multiple treatments
plots <- plot_abeta(
  data = MSD_data, split_by = "Measurement", Treatment = c("Untreated", "BSI", "GSI")
)

arranged_plot <- ggpubr::ggarrange(
  plotlist = plots,
  align = "v",
  ncol = 1
)

# Save data ---------------------------------------------------------------

ggsave(plot = arranged_plot, 
       filename = "06_MSD_abeta_plot.svg", 
       path = here::here("results", "plots"), 
       width = 10, 
       height = 10, 
       dpi = 600, 
       bg = "white"
)
