

library(tidyverse)
library(magrittr)
library(patchwork)
library(pheatmap)
library(sf)
library(sp)
library(concaveman)
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



