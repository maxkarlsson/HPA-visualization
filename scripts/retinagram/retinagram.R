

library(tidyverse)

# Source the HPA visualization library from github: maxkarlsson/HPA-visualization
source("../HPA_visualization_functions.R")

# read single cell consensus data. Three columns: ensg_id, cell_type_name, exp
sc_data <- read_tsv("single_cell_aggregated_data_v3_103.tsv")

# Read file with colors
consensus_colors <- 
  read_tsv("colors_consensus.tsv") %>% 
  select(sample, color) %>% 
  deframe()


single_cell_names <- 
  unique(sc_data$cell_type_name)

if(!all(single_cell_names %in% names(consensus_colors))) {
  stop(paste("Some cells not in palette:", paste(single_cell_names[!single_cell_names %in% names(consensus_colors)],
             collapse = ", ")))
}


sc_hclust <- 
  sc_data %>% 
  spread(cell_type_name, exp) %>% 
  column_to_rownames("ensg_id") %>% 
  cor(method = "spearman") %>% 
  {1 - .} %>% 
  as.dist() %>% 
  hclust(method = "ward.D2")

plot_angle <- 
  sc_hclust %>% 
  calculate_retina_cut_angle()


plot <- 
  plot_angle %>% 
  lapply(function(angle) {
    circular_dendrogram_retinastyle(sc_hclust,
                                    consensus_colors,
                                    preserve_height = F,
                                    scale_expansion = c(0.5, 0.5), 
                                    text_size = 2.5, 
                                    width_range = c(0.75, 4), 
                                    arc_strength = 0.3, 
                                    default_color = "gray80", 
                                    rotate_angle = angle, 
                                    shrink_angle = pi,
                                    flip_text = T, 
                                    text_vjust = 0.5,
                                    text_hnudge = 0.025,
                                    elbow = F)
  })
plot[[1]]
ggsave("sc circular dendro angle 1 alt 1.pdf",
       width = 8, height = 8)
plot[[2]]
ggsave("sc circular dendro angle 2 alt 1.pdf",
       width = 8, height = 8)


plot <- 
  plot_angle %>% 
  lapply(function(angle) {
    circular_dendrogram_retinastyle(sc_hclust,
                                    consensus_colors,
                                    preserve_height = T,
                                    scale_expansion = c(0.5, 0.5), 
                                    text_size = 2.5, 
                                    width_range = c(0.75, 4), 
                                    arc_strength = 0.3, 
                                    default_color = "gray80", 
                                    rotate_angle = angle, 
                                    shrink_angle = pi,
                                    flip_text = T, 
                                    text_vjust = 0.5,
                                    text_hnudge = 0.025,
                                    elbow = F)
  })
plot[[1]]
ggsave("sc circular dendro angle 1 alt 2.pdf",
       width = 8, height = 8)
plot[[2]]
ggsave("sc circular dendro angle 2 alt 2.pdf",
       width = 8, height = 8)



plot <- 
  plot_angle %>% 
  lapply(function(angle) {
    circular_dendrogram_retinastyle(sc_hclust,
                                    consensus_colors,
                                    preserve_height = T,
                                    scale_expansion = c(0.5, 0.5), 
                                    text_size = 2.5, 
                                    width_range = c(0.5, 0.5), 
                                    arc_strength = 0, 
                                    default_color = "gray80", 
                                    rotate_angle = angle, 
                                    shrink_angle = pi,
                                    flip_text = T, 
                                    text_vjust = 0.5,
                                    text_hnudge = 0.025,
                                    elbow = T)
  })
plot[[1]]
ggsave("sc circular dendro angle 1 alt 3.pdf",
       width = 8, height = 8)
plot[[2]]
ggsave("sc circular dendro angle 2 alt 3.pdf",
       width = 8, height = 8)
