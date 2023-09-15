# Load libraries ----------------------------------------------------------

library(tidyverse)
library(here)
library(ggforce)

# Load data ---------------------------------------------------------------

sqanti_class <-
  read.table(here::here("raw_data", "IsoSeq", "APP_classification_processed.txt"), header = TRUE, sep="\t")


# Functions ---------------------------------------------------------------

transcripts_depreciation_curve <-
  
  function(sqanti_output, gene_name) {
    
    # Generate data frame with number of transcripts by each given threshold
    n_samples <- sqanti_output %>%
      dplyr::select(starts_with("FL.")) %>% 
      ncol()
    
    counts.long <- sqanti_output %>%
      dplyr::select(isoform, 
                    starts_with("NFLR")) %>%
      pivot_longer(cols = c(-isoform), 
                   names_to = "sample", 
                   values_to = "read_count")
    
    counts.long$sample <- counts.long$sample %>% 
      str_replace("NFLR.*\\.", "") %>% 
      str_replace("_3p", "")
    
    threshs <- seq(from = 0, to = max(sqanti_output$NFLR_mean), by = 0.1)
    
    res <- data.frame()
    
    for(i in threshs){
      transcripts_per_sample <- counts.long %>%
        dplyr::filter(read_count >=i) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(n_transcripts = dplyr::n(), 
                         .groups = "keep") %>%
        dplyr::mutate(NFLR_threshold = i)
      
      res <- rbind(res, transcripts_per_sample)
      
      rm(transcripts_per_sample)
      
    }
    
    df <- res %>%
      dplyr::filter(sample %in% "NFLR_mean")
    
    sd <- res %>%
      dplyr::filter(!sample %in% "NFLR_mean") %>%
      dplyr::group_by(NFLR_threshold) %>%
      dplyr::summarise(sd = sd(n_transcripts), 
                       .groups = "keep")
    
    df <- left_join(df, 
                    sd, 
                    by = "NFLR_threshold")
    
    # Plot data
    plot <-
      df %>% 
      ggplot(aes(x=NFLR_threshold, 
                 y = n_transcripts)) +
      geom_line(group=1) +
      geom_ribbon(aes(NFLR_threshold, 
                      ymax = n_transcripts + sd, 
                      ymin = n_transcripts - sd), 
                  alpha = 0.2, 
                  fill = "blue", group=1) +
      ggtitle(paste0("Depreciation curve of ", gene_name, " transcripts")) +
      scale_x_continuous(name = expression(paste("Normalised expression (", NFLR[T], ")"))) +
      scale_y_continuous(name = "No. of transcripts detected") +
      guides(fill="none") +
      facet_zoom(xlim = c(0, 5),
                 ylim = c(0, 50),
                 horizontal = F) +
      
      theme(plot.title = element_text(face = "bold",
                                      size = 16,
                                      hjust = 0.5),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", size=0.4),
            axis.title = element_text(size = 14),
            axis.text.x = element_text(face = "bold",
                                       size = 10),
            axis.text.y = element_text(face = "bold",
                                       size = 10),
            axis.title.y = element_text(size = 14),
            strip.text.x = element_text(face = "bold",
                                        size = 12),
            legend.title = element_blank())
    
    return(plot)
  }

# Main --------------------------------------------------------------------

depreciation_curve_plot <-
  transcripts_depreciation_curve(sqanti_output = sqanti_class, 
                                 gene_name = "APP")

# Save data ---------------------------------------------------------------
  
ggsave(plot = depreciation_curve_plot,
       filename = "Depreciation_curve_plot.png", 
       path = here::here("results"), 
       width = 6, 
       height = 4, 
       dpi = 600
)

