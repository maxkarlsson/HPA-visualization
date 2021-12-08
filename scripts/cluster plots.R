
source("scripts/cluster rasterizer.R")
load("data/louvain_clustering.Rdata")
load("data/pig_umap.Rdata")

plot_data <- 
  pig_umap %>% 
  rename(gene_name = features) %>%
  left_join(louvain_clustering,
            by = c("gene_name" = "gene")) %>% 
  mutate(cluster = factor(cluster))

plot_palette <- 
  set_names(colorRampPalette(ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13))(n_distinct(plot_data$cluster)),
            sort(unique(plot_data$cluster)))


pig_class <- 
  read_tsv("data/pig_compare_category_piggenome_92.tsv")

#### ------ gene -------
cluster_hulls <- 
  plot_data %$%
  generate_cluster_hulls(V1,
                         V2, 
                         element_id = gene_name,
                         cluster_membership = cluster, 
                         n = 500)

cluster_hulls$data %>% 
  ggplot(aes(V1, V2, color = sub_type)) +
  geom_point(size = 0.1) + 
  theme_bw()



cluster_hulls$density %>% 
  filter(cum_z < 0.95) %>% 
  ggplot(aes(z, cum_z, group = paste(cluster, sub_cluster)))+
  geom_line(alpha = 0.2) +
  facet_wrap(~density_type, ncol = 1) +
  theme_bw() +
  scale_x_log10()

cluster_hulls$density %>% 
  filter(cum_z < 0.95) %>%
  ggplot(aes(z)) +
  geom_density() +
  facet_wrap(~density_type, ncol = 1) +
  scale_x_log10() +
  theme_bw()


cluster_hulls$density_landmass %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster)) +
  geom_hex(data = plot_data,
           aes(V1, V2), 
           fill = "gray",
           bins = 100) +
  geom_tile() +
  facet_wrap(~density_type, ncol = 1) +
  theme_bw() + 
  coord_fixed() +
  scale_fill_manual(values = plot_palette)

cluster_hulls$density_landmass %>% 
  filter(cluster %in% c("30", "13")) %>%
  ggplot(aes(x_coord, y_coord, fill = cluster)) +
  geom_tile() +
  facet_wrap(~density_type, ncol = 1) +
  theme_bw() + 
  coord_fixed() +
  scale_fill_manual(values = plot_palette)


cluster_hulls$density_landmass %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster)) +
  geom_hex(data = plot_data,
           aes(V1, V2), 
           fill = "gray",
           bins = 100) +
  geom_tile() +
  facet_wrap(~density_type, ncol = 1) +
  theme_bw() + 
  coord_fixed() +
  scale_fill_manual(values = plot_palette)




cluster_hulls$hulls %>% 
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  geom_point(data = cluster_hulls$data %>% 
               filter(sub_type != "primary"),
             aes(V1, V2),
             inherit.aes = F,
             fill = "black",
             size = 0.1) +
  
  geom_polygon(color = "black", 
               size = 0.1) +
  stat_density_2d(data = plot_data,
                  geom = "polygon", 
                  aes(V1, V2, alpha = log10(..level..)),
                  h = 0.2 * 2,
                  n = 100,
                  inherit.aes = F,
                  fill = "white") +
  facet_wrap(~density_type, ncol = 1) +
  theme_bw() + 
  coord_fixed() +
  scale_fill_manual(values = plot_palette) +
  scale_alpha_continuous(range = c(0, 0.5))
ggsave("results/UMAP/Final hulls.pdf",
       width = 10, height = 10)


cluster_hulls$hulls %>% 
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  
  geom_polygon(color = "black", 
               size = 0.1,
               alpha = 0.5) +
  facet_wrap(~density_type, ncol = 1) +
  theme_bw() + 
  coord_fixed() +
  scale_fill_manual(values = plot_palette) +
  scale_alpha_continuous(range = c(0, 0.5))
ggsave("results/UMAP/Final hulls alpha.pdf",
       width = 10, height = 10)




p <- 
  cluster_hulls$hulls %>% 
  filter(density_type == "primary") %>%
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  geom_hex(data = cluster_hulls$data %>% 
             filter(sub_type != "outlier") %>%
             select(1:3),
           aes(V1, V2, fill = 1),
           bins = 100,
           inherit.aes = F,
           fill = "lightgray") +
  geom_polygon(color = "black", 
               show.legend = F,
               size = 0.1) +
  facet_wrap(~cluster, scales ="free") +
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  # coord_fixed() +
  scale_fill_manual(values = plot_palette)
ggsave("results/UMAP/final hulls separated.pdf",plot = p,
       width = 16, height = 12)





p1 <-
  cluster_hulls$data %>% 
  ggplot(aes(V1, V2, color = cluster)) +
  geom_point(size = 0.1, 
             show.legend = F) +
  theme_void() +
  coord_fixed() +
  scale_color_manual(values = plot_palette) +
  ggtitle("All points")

p2 <-
  cluster_hulls$data %>% 
  filter(sub_type != "outlier") %>%
  ggplot(aes(V1, V2, color = cluster)) +
  geom_point(size = 0.1, 
             show.legend = F) +
  theme_void() +
  coord_fixed() +
  scale_color_manual(values = plot_palette) +
  ggtitle("No outliers")

p3 <- 
  cluster_hulls$hulls %>% 
  filter(density_type == "primary secondary") %>%
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  geom_polygon(color = "black", 
               show.legend = F,
               size = 0.1) +
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  scale_fill_manual(values = plot_palette) +
  ggtitle("Complete hull")


p4 <- 
  cluster_hulls$hulls %>% 
  filter(density_type == "primary") %>%
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  geom_polygon(color = "black", 
               show.legend = F,
               size = 0.1) +
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  scale_fill_manual(values = plot_palette) +
  ggtitle("Mainland hull")


p1 | p2 | p3 | p4
ggsave("results/UMAP/UMAP graphical representation levels.pdf",
       width = 20, height = 6)




cluster_hulls$hulls %>% 
  filter(density_type == "primary secondary") %>%
  ggplot(aes(X, Y, group = paste(cluster, sub_cluster, landmass))) +
  geom_polygon(data = . %>% 
                 filter(cluster != 22),
               color = "gray",
               fill = "lightgray",
               show.legend = F,
               size = 0.1) +
  geom_point(data = cluster_hulls$data %>% 
               filter(cluster == 22),
             aes(V1, V2),
             inherit.aes = F,
             color = "red",
             size = 0.1, 
             show.legend = F) +
  geom_polygon(data = . %>% 
                 filter(cluster == 22),
               color = "black", 
               fill = NA,
               show.legend = F,
               size = 0.1) +
  

  
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  scale_fill_manual(values = plot_palette) +
  ggtitle("Highlight one group")
ggsave("results/UMAP/UMAP cluster highlight.pdf",
       width = 6, height = 6)


specificity_cluster_meta <- 
  louvain_clustering %>%
  group_by(cluster) %>%
  mutate(cluster_n = length(gene)) %>%
  ungroup() %>%
  left_join(pig_class,
            by = c("gene" = "ensg_id")) %>% 
  separate_rows(enhanced_tissues, sep = ",") %>% 
  mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues),
                                   "not enriched",
                                   enhanced_tissues),
         enhanced_tissues = factor(enhanced_tissues),
         cluster = factor(cluster)) %>%
  group_by(cluster, cluster_n, enhanced_tissues, .drop = F) %>% 
  count %>%
  ungroup() %>%
  mutate(fraq = n / cluster_n)
  
p1 <- 
  cluster_hulls$hulls %>% 
  filter(density_type == "primary secondary") %>%
  left_join(specificity_cluster_meta %>% 
              filter(enhanced_tissues == "testis")) %>% 
  ggplot(aes(X, Y, fill = fraq, group = paste(cluster, sub_cluster, landmass))) +
  
  geom_polygon(color = "black", 
               show.legend = F,
               size = 0.1) +
  
  
  
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  scale_fill_gradient(low = "white", high = "orangered") +
  # scale_fill_manual(values = plot_palette) +
  ggtitle("Highlight meta data per group")

p2 <- 
  cluster_hulls$data %>% 
  left_join(pig_class,
            by = c("element_id" = "ensg_id")) %>% 
  mutate(is_testis = grepl("testis", enhanced_tissues)) %>%
  ggplot(aes(V1, V2, color = is_testis)) +
  
  geom_point(show.legend = F,
             size = 0.1) +
  
  
  
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  # scale_fill_gradient(low = "white", high = "orangered") +
  scale_color_manual(values = c("TRUE" = "orangered", "FALSE" = "gray")) +
  ggtitle("Highlight meta data per element")

p2 | p1
ggsave("results/UMAP/UMAP cluster highlight metadata.pdf",
       width = 12, height = 6)


#### ------ cell images -------

cell_image_umap <- 
  read_csv("data/umap_results_fit_all_transform_all_sorted_20190422.csv") %>% 
  mutate(location_cluster = ifelse(grepl(" ", location_code),
                                   "multi", 
                                   location_code)) 



plot_palette_cell <- 
  set_names(colorRampPalette(ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13))(n_distinct(cell_image_umap$location_cluster)),
            sort(unique(cell_image_umap$location_cluster)))

#### ------ gene -------
cluster_hulls_cell <- 
  cell_image_umap %$%
  generate_cluster_hulls(x,
                         y, 
                         element_id = id,
                         cluster_membership = location_cluster, 
                         n = 500)

p1 <-
  cluster_hulls_cell$data %>% 
  arrange(rev(cluster)) %>%
  ggplot(aes(V1, V2, color = cluster)) +
  geom_point(size = 0.1, 
             show.legend = F) +
  theme_void() +
  coord_fixed() +
  scale_color_manual(values = plot_palette_cell) +
  ggtitle("All points")

p2 <-
  cluster_hulls_cell$data %>% 
  arrange(rev(cluster)) %>%
  filter(sub_type != "outlier") %>%
  ggplot(aes(V1, V2, color = cluster)) +
  geom_point(size = 0.1, 
             show.legend = F) +
  theme_void() +
  coord_fixed() +
  scale_color_manual(values = plot_palette_cell) +
  ggtitle("No outliers")

p3 <- 
  cluster_hulls_cell$hulls %>% 
  filter(density_type == "primary secondary") %>%
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  geom_polygon(color = "black", 
               show.legend = F,
               size = 0.1) +
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  scale_fill_manual(values = plot_palette_cell) +
  ggtitle("Complete hull")


p4 <- 
  cluster_hulls_cell$hulls %>% 
  filter(density_type == "primary") %>%
  ggplot(aes(X, Y, fill = cluster, group = paste(cluster, sub_cluster, landmass))) +
  geom_polygon(color = "black", 
               show.legend = F,
               size = 0.1) +
  theme_void() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(0, "mm")) +
  coord_fixed() +
  scale_fill_manual(values = plot_palette_cell) +
  ggtitle("Mainland hull")


p1 | p2 | p3 | p4
ggsave("results/UMAP/UMAP cellatlas graphical representation levels.pdf",
       width = 20, height = 6)
