
library(tidyverse)
library(magrittr)
library(patchwork)
library(pheatmap)
load("data/louvain_clustering.Rdata")
load("data/pig_umap.Rdata")
library(rayshader)


get_density <- 
  function(x, y, h = 0.5, n = 100, lims = c(range(x),
                                            range(y))) {
    g_density <- 
      MASS::kde2d(x, y, h = h, n = n, 
                  lims = lims) 
    
    
    g_density$z %>%
      as_tibble(rownames = "x") %>%
      gather(y, z, -x) %>% 
      mutate(y = gsub("V", "", y)) %>% 
      mutate_at(vars(x,y), as.integer) %>%
      left_join(enframe(g_density$x, "x", "x_coord"),
                by = "x") %>% 
      left_join(enframe(g_density$y, "y", "y_coord"),
                by = "y")
  }

plot_data <- 
  pig_umap %>% 
  rename(gene_name = features) %>%
  left_join(louvain_clustering,
            by = c("gene_name" = "gene")) %>% 
  mutate(cluster = factor(cluster))

plot_palette <- colorRampPalette(ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13))(n_distinct(plot_data$cluster))

plot_range <- 
  c(range(plot_data$V1),
    range(plot_data$V2)) * 1.05


plot_density <- 
  plot_data %>%
  group_by(cluster) %>%
  do({
    get_density(.$V1, 
                .$V2, 
                h = 0.5, 
                n = 300, 
                lims = plot_range) 
    
  }) %>% 
  ungroup()


plot_settings <- 
  crossing(n = 300,
           h = seq(0.15, 0.25, length.out = 3),
           z_lim = 10^seq(-2,-5))



plot_density_all <-
  plot_settings %>%
  select(h, n) %>%
  distinct() %>%
  group_by_all() %>%
  do({
    h <- .$h
    n <- .$n
    
    plot_data %>%
      group_by(cluster) %>%
      do({
        get_density(.$V1, 
                    .$V2, 
                    h = h, 
                    n = n, 
                    lims = plot_range) 
        
      })
  }) %>%
  ungroup() %>%
  left_join(plot_settings) %>%
  group_by(h, n, z_lim, cluster) %>%
  mutate(z_rel = z / max(z)) %>% 
  ungroup() %>% 
  filter(z_rel > z_lim)


   
plot_density_total <- 
  plot_settings %>%
  select(h, n) %>%
  distinct() %>%
  group_by_all() %>%
  do({
    h <- .$h
    n <- .$n
    plot_data %$%
      get_density(V1, 
                  V2, 
                  h = h, 
                  n = n, 
                  lims = plot_range) 
  }) %>%
  ungroup() %>%
  mutate(z_rel = z / max(z))
  
plot_density_total %>%
  left_join(plot_settings) %>%
  
  group_by(h, n, z_lim) %>% 
  mutate(z_rel = z / max(z)) %>% 
  ungroup() %>%
  filter(z_rel > z_lim) %>%
  ggplot(aes(x_coord, y_coord, fill = log10(z_rel + 1)))+
  geom_tile() + 
  scale_fill_viridis_c() +
  coord_fixed() +
  facet_grid(h ~ z_lim) +
  theme_bw()
ggsave("results/UMAP/UMAP density settings.pdf", width = 8, height = 6)


plot_density %>%
  group_by(cluster) %>% 
  mutate(z = z / max(z)) %>%
  filter(z > 0.05) %>%
  group_by(x, y) %>% 
  top_n(1, z) %>%
  slice(1) %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster)) +
  geom_tile() +
  theme_bw() + 
  coord_fixed() +
  scale_fill_manual(values = plot_palette)


p <- 
  plot_density_all %>% 
  group_by(h, n, z_lim, x, y) %>% 
  top_n(1, z) %>%
  slice(1) %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster))+
  geom_tile(show.legend = F) + 
  scale_fill_manual(values = plot_palette) +
  coord_fixed() +
  facet_grid(h ~ z_lim) +
  theme_bw()

ggsave("results/UMAP/UMAP density settings cluster.pdf", plot = p, width = 8, height = 6)


###############################3
plot_density_all %>% 
  filter(h == 0.2, z_lim == 0.001) %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster))+
  geom_tile(show.legend = F) + 
  scale_fill_manual(values = plot_palette) +
  coord_fixed() +
  facet_wrap(~ cluster) +
  theme_bw()

plot_density_all %>% 
  filter(h == 0.2, z_lim == 0.001) 


plot_data_sub <-
  plot_data %>%
  group_by(cluster) %>%
  mutate(sub_cluster = data.frame(V1, V2) %>%
           fpc::dbscan(eps = 0.2) %$%
           cluster) %>%
  ungroup()

subcluster_count <- 
  plot_data_sub %>% 
  group_by(cluster, sub_cluster) %>% 
  count() %>% 
  group_by(cluster) %>% 
  mutate(frac = n / sum(n))

keep_subcluster_fraction <- 
  subcluster_count %>% 
  filter(frac > 0.1) %>%
  filter(sub_cluster != 0) 

keep_subclusters <- 
  subcluster_count %>%
  top_n(1, n) %>% 
  slice(1) %>% 
  ungroup()
  
  
  
plot_data_sub %>% 
ggplot(aes(V1, V2, color = factor(cluster))) +
  geom_point(show.legend = F) +
  scale_color_manual(values = plot_palette) +
  coord_fixed() +
  theme_bw() +
  ggtitle("All points") +
  
  plot_data_sub %>% 
  filter(sub_cluster != 0) %>%
  ggplot(aes(V1, V2, color = factor(cluster))) +
  geom_point(show.legend = F) +
  scale_color_manual(values = plot_palette) +
  coord_fixed() +
  theme_bw() +
  ggtitle("No outliers") +
  
  plot_data_sub %>% 
  inner_join(keep_subcluster_fraction) %>%
  ggplot(aes(V1, V2, color = factor(cluster))) +
  geom_point(show.legend = F) +
  scale_color_manual(values = plot_palette) +
  coord_fixed() +
  theme_bw() +
  ggtitle("Only subcluster of 10% total") +
  
  plot_data_sub %>% 
  inner_join(keep_subclusters) %>%
  ggplot(aes(V1, V2, color = factor(cluster))) +
  geom_point(show.legend = F) +
  scale_color_manual(values = plot_palette) +
  coord_fixed() +
  theme_bw() +
  ggtitle("Only main cluster")
ggsave("results/UMAP/subcluster removal.pdf", width = 8, height = 8)  


plot_data_sub_frac <- 
  plot_data_sub %>%
  inner_join(keep_subcluster_fraction) 
  
plot_density_frac_all <-
  plot_settings %>%
  select(h, n) %>%
  distinct() %>%
  group_by_all() %>%
  do({
    h <- .$h
    n <- .$n
    
    plot_data_sub_frac %>%
      group_by(cluster) %>%
      do({
        get_density(.$V1, 
                    .$V2, 
                    h = h, 
                    n = n, 
                    lims = plot_range) 
        
      })
  }) %>%
  ungroup() %>%
  left_join(plot_settings) %>%
  group_by(h, n, z_lim, cluster) %>%
  mutate(z_rel = z / max(z)) %>% 
  ungroup() %>% 
  filter(z_rel > z_lim)

p <- 
  plot_density_frac_all %>% 
  group_by(h, n, z_lim, x, y) %>% 
  top_n(1, z) %>%
  slice(1) %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster))+
  geom_tile(show.legend = F) + 
  scale_fill_manual(values = plot_palette) +
  coord_fixed() +
  facet_grid(h ~ z_lim) +
  theme_bw()

ggsave("results/UMAP/UMAP density settings cluster frac.pdf", plot = p, width = 8, height = 6)



plot_density_frac_all %>% 
  filter(h == 0.2 & z_lim == 0.001) %>%
  group_by(h, n, z_lim, x, y) %>% 
  top_n(1, z) %>%
  slice(1) %>% 
  ggplot(aes(x_coord, y_coord, fill = cluster))+
  geom_tile(show.legend = F) + 
  scale_fill_manual(values = plot_palette) +
  coord_fixed() +
  facet_wrap(~ cluster) +
  theme_bw()
ggsave("results/UMAP/UMAP density cluster frac.pdf", width = 10, height = 10)
