
setwd("scripts/retinagram/")
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

# Calculate the two cut angles that doesn't interfere with edges
plot_angle <- 
  sc_hclust %>% 
  calculate_retina_cut_angle()

# Explanation to retinagram arguments: 
# clust             hclust object
# color_pal         named vector with hclust leaves as names and colors as values
# preserve_height   boolean whether to preserve hclust node heights or scale them equidistant across the plot radius
# scale_expansion   vector of length two with a factor to scale the plot axes with
# text_size         leaf text size
# width_range       vector of length two (or one) with the size range of edges
# arc_strength      the strength of the edge curvature (usually between 0 and 1)
# default_color     default color of edges
# rotate_angle      the angle (in radians) to rotate the tree. If shrink angle is not 0, run calculate_retina_cut_angle() to find angles that don't break the plot
# shrink_angle      the angle (in radians) to shrink the tree. 0 will generate a cirle, pi will generate a half circle, and so on
# flip_text         whether to flip the upside down. Useful if you intend to rotate the tree
# text_vjust        the vertical justification of the text
# text_hnudge       a numeric value to nudge leaf text to give a small gap between leaf texts and tree 
# elbow             boolean whether to plot elbows rather than curved edges


plot <- 
  plot_angle %>% 
  lapply(function(angle) {
    circular_dendrogram_retinastyle(sc_hclust,
                                    consensus_colors,
                                    preserve_height = F,
                                    scale_expansion = c(0.5, 0.5), 
                                    text_size = 2.5, 
                                    width_range = c(0.75, 4), 
                                    arc_strength = 0.5, 
                                    default_color = "gray80", 
                                    rotate_angle = angle, 
                                    shrink_angle = pi,
                                    flip_text = T, 
                                    text_vjust = 0.5,
                                    text_hnudge = 0.025,
                                    elbow = F)
  })
plot[[1]]
ggsave("sc circular dendro angle 1 alt 1.svg",
       width = 8, height = 8)
plot[[2]]
ggsave("sc circular dendro angle 2 alt 1.svg",
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
                                    arc_strength = 0.5, 
                                    default_color = "gray80", 
                                    rotate_angle = angle, 
                                    shrink_angle = pi,
                                    flip_text = T, 
                                    text_vjust = 0.5,
                                    text_hnudge = 0.025,
                                    elbow = F)
  })
plot[[1]]
ggsave("sc circular dendro angle 1 alt 2.svg",
       width = 8, height = 8)
plot[[2]]
ggsave("sc circular dendro angle 2 alt 2.svg",
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
ggsave("sc circular dendro angle 1 alt 3.svg",
       width = 8, height = 8)
plot[[2]]
ggsave("sc circular dendro angle 2 alt 3.svg",
       width = 8, height = 8)

