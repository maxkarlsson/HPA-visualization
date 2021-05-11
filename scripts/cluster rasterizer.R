

library(tidyverse)
library(magrittr)
library(patchwork)
library(pheatmap)
library(sf)
library(sp)
library(concaveman)
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

generate_cluster_hulls <- 
  function(V1,
           V2, 
           element_id,
           cluster_membership, 
           n = 300,
           cum_z_lim = 0.95,
           frac_lim = 0.05,
           plot_range = c(range(V1),
                          range(V2)) * 1.05,
           bandwidth_fct = 1) {
    
    cluster_data <- 
      tibble(V1, V2, 
             element_id,
             cluster = cluster_membership)
    
    
    plot_range_tb <- 
      set_names(plot_range,
                c("xmin", 
                  "xmax",
                  "ymin", 
                  "ymax")) %>% 
      enframe() %>% 
      spread(name, value)
    
    plot_area <- 
      (plot_range_tb$xmax - plot_range_tb$xmin) * 
      (plot_range_tb$ymax - plot_range_tb$ymin)
    
    plot_bandwidth <- 
      bandwidth_fct * plot_area / 1000
    
    # Find subclusters
    plot_data_sub <-
      cluster_data %>%
      group_by(cluster) %>%
      mutate(sub_cluster = data.frame(V1, V2) %>%
               fpc::dbscan(eps = plot_bandwidth) %$%
               cluster) %>%
      ungroup() %>% 
      group_by(cluster, sub_cluster) %>%
      mutate(n_sub_genes = n_distinct(element_id)) %>% 
      ungroup() %>%
      left_join(select(., cluster, sub_cluster, n_sub_genes) %>% 
                  distinct() %>%
                  group_by(cluster) %>%
                  mutate(sub_type = case_when(sub_cluster == 0 ~ "outlier",
                                              n_sub_genes / sum(n_sub_genes) < frac_lim ~ "outlier",
                                              rank(-n_sub_genes, 
                                                   ties.method = "first") == 1 ~ "primary",
                                              T ~ "secondary")),
                by = c("cluster", "sub_cluster", "n_sub_genes")) %>% 
      ungroup() %>% 
      arrange(cluster)
    
    
    # Calculate plot density
    plot_density <- 
      plot_data_sub %>% 
      filter(sub_type != "outlier") %>%
      bind_rows("primary secondary" = ., 
                "primary" = filter(., sub_type == "primary"),
                .id = "density_type") %>%
      group_by(density_type, cluster, sub_cluster) %>%
      do({
        get_density(.$V1, 
                    .$V2, 
                    h = plot_bandwidth, 
                    n = n, 
                    lims = plot_range) 
        
      }) %>% 
      ungroup() %>% 
      filter(z > 1e-200) %>% 
      group_by(density_type, cluster, sub_cluster) %>% 
      mutate(z = z / sum(z)) %>% 
      arrange(cluster, sub_cluster, -z) %>% 
      mutate(cum_z = cumsum(z)) %>% 
      ungroup()
    
    
    # Filter pixels such that 95% of density is included
    # Each point is then assigned to the cluster with highest density
    plot_density_filtered <- 
      plot_density %>% 
      filter(cum_z < cum_z_lim) %>%
      group_by(density_type, x, y) %>% 
      top_n(1, z) %>%
      slice(1) %>% 
      ungroup()
    
    
    
    
    # Define a main landmass
    plot_density_landmass <-
      plot_density_filtered %>%
      group_by(density_type, cluster, sub_cluster) %>%
      mutate(landmass = data.frame(x_coord, y_coord) %>%
               fpc::dbscan(eps = plot_bandwidth / 2.5) %$%
               cluster) %>%
      group_by(density_type, cluster, sub_cluster, landmass) %>% 
      mutate(n_landmass_points = length(x)) %>%  
      ungroup() %>%
      left_join(select(., density_type, cluster, sub_cluster, landmass, n_landmass_points) %>% 
                  distinct() %>%
                  group_by(density_type, cluster, sub_cluster) %>%
                  
                  mutate(frac_landmass = n_landmass_points / sum(n_landmass_points),
                         landmass_type = case_when(rank(-n_landmass_points, 
                                                        ties.method = "first") == 1 ~ "primary",
                                                   T ~ "secondary"))) %>% 
      ungroup() %>% 
      arrange(cluster)
    
    
    
    
    plot_density_mainland_filtered <-
      plot_density %>% 
      left_join(plot_density_landmass %>% 
                  select(density_type, cluster, sub_cluster, x, y, landmass, frac_landmass, landmass_type)) %>%
      # filter(landmass_type != "secondary") %>%
      filter(frac_landmass > frac_lim) %>%
      filter(cum_z < cum_z_lim) %>%
      group_by(density_type, x, y) %>% 
      top_n(1, z) %>%
      slice(1) %>% 
      ungroup()
    
    ######
    plot_data_hulls <- 
      plot_density_mainland_filtered %>% 
      
      group_by(density_type, cluster, sub_cluster, landmass) %>% 
      do({
        st_as_sf(., coords=c('x_coord','y_coord')) %>%
          concaveman(concavity = 1, length_threshold = plot_bandwidth / 2) %$%
          st_coordinates(polygons) %>% 
          as_tibble()
      }) %>% 
      ungroup()
    
    
    list(data = plot_data_sub,
         hulls = plot_data_hulls,
         density = plot_density,
         density_landmass = plot_density_mainland_filtered)
  }


plot_data <- 
  pig_umap %>% 
  rename(gene_name = features) %>%
  left_join(louvain_clustering,
            by = c("gene_name" = "gene")) %>% 
  mutate(cluster = factor(cluster))

plot_palette <- 
  set_names(colorRampPalette(ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13))(n_distinct(plot_data$cluster)),
            sort(unique(plot_data$cluster)))

#### ------ -------
cluster_hulls <- 
  plot_data %$%
  generate_cluster_hulls(V1,
                         V2, 
                         element_id = gene_name,
                         cluster_membership = cluster, 
                         n = 1000)

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
                  h = plot_bandwidth * 2,
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
ggsave("results/UMAP/final hulls separated.pdf",
       width = 16, height = 12)

# plot_data_hulls_polygons <- 
#   plot_data_hulls %>% 
#   select(density_type, cluster, sub_cluster) %>% 
#   distinct() %>% 
#   unite(id, density_type, cluster, sub_cluster, remove = F)
# 
# plot_data_hulls_spts <-
#   plot_data_hulls %>%
#   distinct() %>%
#   left_join(plot_data_hulls_polygons) %>% 
#   split(.$id) %>% 
#   lapply(function(x) {
#     select(x, X, Y) %>%
#       Polygon() %>%
#       list()%>%
#       Polygons(1) %>%
#       list() %>%
#       SpatialPolygons()
#   }) 
# 
# plot_data_hulls_intersections <- 
#   plot_data_hulls_polygons %>% 
#   group_by(density_type) %>%
#   do({
#     g_ids <<- .
#     unique(g_ids$id) %>%
#     crossing(cluster1 = .,
#              cluster2 = .) %>% 
#       filter(cluster1 < cluster2) %>% 
#       group_by_all() %>% 
#       do({
#         g_data <<- .
#         
#         rgeos::gIntersection(spgeom1 = plot_data_hulls_spts[[g_data$cluster1]],
#                              spgeom2 = plot_data_hulls_spts[[g_data$cluster2]],
#                              byid = FALSE) %>%
#           list() %>%
#           enframe("temp", "spts") %>% 
#           select(-temp)
#         
#       })
#   })
#    
# inter1 <- raster::intersect(extent(r_list[[i]]), extent(p))

