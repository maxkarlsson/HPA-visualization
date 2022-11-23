
library(tidyverse)
library(pammtools)
library(multidplyr)
setwd("scripts/tau/")

calculate_tau_score <- 
  function(wide_data) {
    max_exp <- 
      apply(wide_data,
            MARGIN = 1,
            function(x) max(x, na.rm = T))
    
    N <- 
      apply(wide_data,
            MARGIN = 1,
            function(x) length(which(!is.na(x))))
    
    expression_sum <- 
      wide_data %>% 
      sweep(MARGIN = 1, 
            STATS = max_exp, 
            FUN = `/`) %>% 
      {1 - .} %>% 
      apply(MARGIN = 1,
            function(x) sum(x, na.rm = T))
    
    
    tau_score <- 
      (expression_sum / (N - 1)) %>% 
      enframe("gene", "tau_score")
    
    tau_score
  }


# Read consensus dataset (same data as for classification)
consensus_data <- 
  read_tsv("consensus_data.tsv")


# Scale data before tau calculation
tau_data <-
  consensus_data %>% 
  spread(cell_type_name, exp) %>%  
  column_to_rownames("ensg_id") %>%
  # Log scale data:
  {log10(. + 1)} %>%
  
  # Max-scale data (divide each gene by maximum valus):
  sweep(., 
        MARGIN=1, 
        apply(., 
              MARGIN = 1,
              function(x)
                max(x, na.rm = T))
        , `/`) 

# Tau is calculated for all genes, but tau scores should be removed for genes that are not expressed: 
tau <- 
  tau_data %>% 
  calculate_tau_score()

write_tsv(tau, "tau.tsv")

# Prepare data for tau silhouette plotting
tau_plot_data <- 
  tau_data %>% 
  as_tibble(rownames = "gene") %>% 
  gather(sample, value, -1) %>% 
  arrange(gene, -value) %>%
  group_by(gene) %>% 
  mutate(rank = row_number()) %>% 
  ungroup() %>% 
  left_join(tau,
            by = "gene") %>% 
  group_by(gene) %>% 
  slice(-1) %>% 
  mutate(rank_scaled = scales::rescale(rank, 0:1)) %>% 
  ungroup() %>% 
  select(gene, value, rank_scaled)


# Make an svg per each gene
tau_plot_data %>% 
  group_by(gene) %>%
  do({
    ggplot(., 
           aes(rank_scaled, ymin = value, y = value, ymax = Inf)) +
      
      geom_rect(xmin = 0, ymin = 0, 
                xmax = 1, ymax = 1,
                fill = "gray80",
                color = NA) +
      geom_stepribbon(fill = "gray20") +
      coord_fixed() +
      scale_x_continuous(limits = 0:1,
                         expand = expansion(0)) +
      scale_y_continuous(limits = 0:1, 
                         expand = expansion(0)) +
      theme_void() +
      theme(panel.border = element_rect(fill = NA, color = NA),
            panel.spacing = unit(1, "mm"), 
            strip.text = element_text(hjust = 0, color = NA))
    
    ggsave(paste0("svg/", unique(.$gene), ".svg"), 
           width = 1, height = 1, units = "cm")
    
    tibble()
  }) 