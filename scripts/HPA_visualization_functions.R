
# ------ Data manipulation -----

convert_to_widedata <- 
  function(dat, xcol, ycol, valcol) {
    suppressMessages(require(tidyverse))
    
    dat %>% 
      select(all_of(c(xcol, ycol, valcol))) %>% 
      spread(xcol, valcol) %>% 
      column_to_rownames(ycol)
    
  }

remove_genes <- 
  function(wide_dat,
           rm_NA_genes = F,
           rm_not_expressed = F,
           not_expressed_lim = 1,
           rm_no_variance = F) {
    suppressMessages(require(tidyverse))
    
    
    if(rm_NA_genes) {
      
      wide_dat <- 
        wide_dat %>% 
        {.[complete.cases(.), ]}
    }
    
    if(rm_not_expressed) {
      
      wide_dat <- 
        wide_dat[apply(wide_dat,
                       MARGIN = 1, 
                       max) >= not_expressed_lim,] 
    }
    
    if(rm_not_expressed) {
      
      wide_dat <- 
        wide_dat[apply(wide_dat,
                       MARGIN = 1, 
                       max) >= not_expressed_lim,] 
    }
    
    if(rm_no_variance) {
      
      wide_dat <- 
        wide_dat[apply(wide_dat,
                       MARGIN = 1, 
                       sd) > 0,] 
    }
    
    wide_dat
  }



scale_data <- 
  function(wide_data, 
           logp1_scale = F,
           zscore_scale = F, 
           max_scale = F,
           sum_scale = F) {
    suppressMessages(require(tidyverse))
    
    if(zscore_scale & max_scale) stop("Scaling: Choose either max scaling or z-score scaling")
    
    if(logp1_scale) {
      wide_data <- 
        log10(wide_data + 1)
    }
    if(zscore_scale) {
      wide_data <- 
        wide_data %>% 
        t() %>% 
        scale() %>% 
        t()
    }
    
    if(max_scale) {
      gene_max <- 
        apply(wide_data, 
              MARGIN = 1,
              function(x)
                max(x, na.rm = T))
      
      wide_data <-
        sweep(wide_data, MARGIN=1, gene_max, `/`)
    }
    
    if(sum_scale) {
      gene_sum <- 
        apply(wide_data, 
              MARGIN = 1,
              function(x)
                sum(x, na.rm = T))
      
      wide_data <-
        sweep(wide_data, MARGIN=1, gene_sum, `/`)
    }
    
    wide_data
  }

extract_variable_features <- 
  function(wide_data, 
           log1p_data = F,
           maxvars = 5000,
           ...) {
    
    # Whether to log scale data for variable extraction or not. 
    # This does not affect output data
    if(log1p_data) {
      wide_data <- 
        log1p(wide_data)
    }
    
    # Calculate mean and standard deviation
    vardata <- 
      wide_data %>% 
      {tibble(var = rownames(wide_data), 
              mean = apply(., 
                           MARGIN = 1, 
                           mean),
              sd = apply(., 
                         MARGIN = 1, 
                         sd))}
    
    # Make loess model
    vardata_loess <- 
      loess(sqrt(sd) ~ mean, data = vardata,
            ...)
    
    vardata_res <- 
      vardata %>%
      mutate(pred_sd = predict(vardata_loess, newdata = .),
             high = sqrt(sd) > pred_sd,
             diff = sqrt(sd) - pred_sd,
             top = rank(-diff) < maxvars) 
    
    # vardata_res %T>%
    #   {group_by(., high) %>%
    #       count %>%
    #       print} %>%
    #   mutate(diff = sqrt(sd) - pred_sd,
    #          top = rank(-diff) < maxvars) %>% 
    #   ggplot(aes(mean, sqrt(sd), color  = top)) +
    #   geom_point() +
    #   stat_function(fun = function(x) predict(vardata_loess, newdata = tibble(mean = x)),
    #                 color = "black")
    
    highvar_vars <- 
      vardata_res %>% 
      filter(top) %>% 
      pull(var)
    
    wide_data[highvar_vars,]
    
  }

# ------ Network functions -----

calc_snn <- function(wide_data = NULL, long_data = NULL) {
  
  if((is.null(wide_data) & is.null(long_data)) |
     (!is.null(wide_data) & !is.null(long_data))) {
    stop("Supply either long or wide data")
  }
  
  if(!is.null(wide_data)) {
    
    snn <- 
      wide_data %>% 
      dist(method = "binary") %>% 
      as.matrix() %>% 
      {1 - .} %>% 
      as_tibble(rownames = "var1") %>% 
      gather(var2, snn, -1)
    
    
  } else if (!is.null(long_data)) {
    nn_searchspace <- 
      long_data %>% 
      select(temp1 = 1,
             temp2 = 2) %>% 
      mutate(order = temp1 < temp2,
             var1 = ifelse(order,
                           temp1,
                           temp2),
             var2 = ifelse(order,
                           temp2,
                           temp1)) %>% 
      select(var1, var2) %>% 
      distinct()
    
    stop("Long data is not supported yet.")
    
    }
  
  snn
  
}

# ------ Dimensionality reduction functions -----


do_pca <- 
  function(wide_data, npcs = NULL, ...) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(pcaMethods))
    
    if(is.null(npcs)) {
      npcs <- min(dim(wide_data))
    }
    
    wide_data %>% 
      t() %>% 
      pca(nPcs = npcs, ...) 
  }

get_pca_scores <- 
  function(pca_res, 
           use_R2cum_PCselection = F,
           R2cum_lim = 0.8,
           use_sDev_PCselection = F) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(pcaMethods))
    
    pc_lim <- ncol(pca_res@scores)
    
    if(use_R2cum_PCselection) {
      pc_lim <- 
        which(pca_res@R2cum > R2cum_lim)[1]
    } else if(use_sDev_PCselection) {
      pc_lim <- 
        rev(which(pca_res@sDev >= 1))[1]
    }
    
    pca_res %>% 
      scores() %>% 
      {.[,1:pc_lim]} 
  }

do_umap <- 
  function(wide_data, 
           seed = 42, 
           n_neighbors = 15,
           n_components = 2, 
           metric = "euclidean",
           min_dist = 0.01,
           clean_names = T,
           ret_nn = F,
           ...) {
    suppressMessages(require(tidyverse))
    suppressMessages(require(uwot))
    
    set.seed(seed)
    
    umap_res <- 
      umap(wide_data, 
           n_neighbors = n_neighbors,
           n_components = n_components,
           metric = metric,
           min_dist = min_dist,
           ret_nn = ret_nn,
           ...) %>% 
      as.data.frame()
    
    if(ret_nn) {
      umap_res <-
        list(embedding = umap_res[,1:n_components],
             nn_idx = umap_res[,which(str_detect(colnames(umap_res), "\\.idx"))],
             nn_dist = umap_res[,which(str_detect(colnames(umap_res), "\\.dist"))],
             dist_metric = metric)
    }
    
    if(clean_names) {
      if(ret_nn) {
        
        colnames(umap_res$embedding) <- paste0("UMAP", 1:ncol(umap_res$embedding))
        
        colnames(umap_res$nn_idx) <- paste0("nn", 1:ncol(umap_res$nn_idx))
        
        colnames(umap_res$nn_dist) <- paste0("nn", 1:ncol(umap_res$nn_dist))
        
      } else {
        colnames(umap_res) <- paste0("UMAP", 1:ncol(umap_res))
      } 
    }
    
    umap_res
  }

feature_umap_plot <- 
  function(wide_data,
           n_neighbors = 15,
           min_dist = 0.01) {
    
    wide_data_scaled <- 
      wide_data %>% 
      t() %>% 
      scale_data(zscore_scale = T) %>% 
      t()
    
    umap_res <-
      wide_data_scaled %>% 
      do_umap(n_neighbors = n_neighbors,
              min_dist = min_dist)
    
    
    
    umap_res %>% 
      as_tibble(rownames = "row_id") %>% 
      left_join(wide_data_scaled %>% 
                  as_tibble(rownames = "row_id") %>% 
                  gather(col, value, -row_id),
                by = "row_id") %>% 
      ggplot(aes(UMAP1, UMAP2, color = value)) +
      geom_point() +
      facet_wrap(~col) +
      theme_bw() +
      coord_fixed() +
      scale_color_viridis_c()
    
  }



HPA_sample_UMAP <- 
  function(wide_data) {
    wide_data %>% 
      # Remove genes that have max expression < 1 in dataset, 
      # and genes without variance in expression:
      remove_genes(rm_NA_genes = T,
                   rm_not_expressed = T, 
                   not_expressed_lim = 1,
                   rm_no_variance = T) %>% 
      
      # Extract highly variable variables: 
      # extract_variable_features(log1p_data = T) %>%
      
      # Scale data genewise, first log10(exp + 1), 
      # and then z-score:
      scale_data(logp1_scale = T, 
                 zscore_scale = T) %>% 
      
      # Perform PCA:
      do_pca() %>% 
      
      # Get PCA scores, selecting PCs that constitute 80% cumulative r2:
      get_pca_scores(use_R2cum_PCselection = T,
                     R2cum_lim = 0.8) %>% 
      
      # Perform UMAP with default settings:
      do_umap(n_neighbors = round(sqrt(dim(.)[1]))) %>% 
      
      # Add rownames as a column for saving:
      as_tibble(rownames = "sample")
  }

clusterbundle <- 
  function(x1, y1, x2, y2, 
           dist_method = "euclidean", 
           clustering_method = "average", 
           k = 15) {
    
    nodes <- 
      tibble(x = c(x1, x2),
             y = c(y1, y2)) %>% 
      distinct() %>% 
      mutate(node = as.character(row_number()))
    
    nodes_clust <- 
      nodes %>% 
      column_to_rownames("node") %>% 
      dist(method = dist_method) %>% 
      hclust(method = clustering_method) %>% 
      cutree(k = k) %>% 
      enframe("node", "cluster") %>% 
      mutate(cluster = as.character(cluster))
    
    nodes %>% 
      left_join(nodes_clust,
                by = "node") %>% 
      ggplot(aes(x, y, color = cluster, fill = cluster)) +
      geom_encircle(show.legend = F,
                    alpha = 0.2, 
                    expand = 0) +
      geom_point(show.legend = F)
    
    edges <- 
      tibble(x1, y1, 
             x2, y2) %>% 
      left_join(nodes, 
                by = c("x1" = "x",
                       "y1" = "y")) %>% 
      left_join(nodes, 
                by = c("x2" = "x",
                       "y2" = "y"),
                suffix = c("1", "2")) %>% 
      left_join(nodes_clust,
                by = c("node1" = "node")) %>% 
      left_join(nodes_clust,
                by = c("node2" = "node"),
                suffix = c("1", "2"))
    
    edges %>% 
      filter(cluster1 != cluster2) %>% 
      group_by(cluster1, cluster2) %>% 
      summarise(x1 = mean(x1),
                y1 = mean(y1),
                x2 = mean(x2),
                y2 = mean(y2),
                n = length(node1), 
                .groups = "drop") 
  }

# ------ Batch effect functions -----

do_limma <- 
  function(wide_data, 
           batch = NULL, 
           design = matrix(1,ncol(wide_data),1)) {
    
    require(limma)
    wide_data %>% 
      removeBatchEffect(batch = batch,
                        design = design,
                        method = 'ls') 
    
  }


# ------ Clustering functions -----

get_clustering_order <- 
  function(wide_data,
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T) {
    suppressMessages(require(tidyverse))
    order_row <- 
      rownames(wide_data)
    order_col <- 
      colnames(wide_data)
    
    
    if(cluster_rows) {
      order_row <- 
        wide_data %>% 
        do_cluster_rows(distance_method = distance_method,
                        clustering_method = clustering_method) %>% 
        get_clustering_labels()
    }
    
    if(cluster_cols) {
      order_col <- 
        wide_data %>% 
        t() %>% 
        do_cluster_rows(distance_method = distance_method,
                        clustering_method = clustering_method) %>% 
        get_clustering_labels()
    }
    
    list(row = order_row, 
         col = order_col)
  }

cluster_wide_data <- 
  function(wide_data,
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T) {
    suppressMessages(require(tidyverse))
    
    clustering_order <- 
      get_clustering_order(wide_data,
                           distance_method = distance_method,
                           clustering_method = clustering_method, 
                           cluster_rows = cluster_rows,
                           cluster_cols = cluster_cols)
    
    wide_data[clustering_order$row, clustering_order$col] 
    
  }

cluster_long_data <-  
  function(long_data,
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T, 
           fill = NA) {
    suppressMessages(require(tidyverse))
    
    wide_data <- 
      long_data %>% 
      select(1:3) %>% 
      spread(2, 3, fill = fill) %>% 
      column_to_rownames(names(long_data)[1])
    
    order_row <- 
      rownames(wide_data)
    order_col <- 
      colnames(wide_data)
    
    
    if(cluster_rows) {
      order1 <- 
        wide_data %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    if(cluster_cols) {
      order2 <- 
        wide_data %>% 
        t() %>% 
        dist(method = distance_method) %>% 
        hclust(method = clustering_method) %>% 
        with(labels[order])
    }
    
    long_data %>% 
      rename(v1 = 1, 
             v2 = 2,
             val = 3) %>% 
      mutate(v1 = factor(v1, order1),
             v2 = factor(v2, order2)) %>% 
      set_names(names(long_data))
    
  }

cluster_long_data_grouped <- 
  function(long_data,
           group_col = "ds2",
           distance_method = "euclidean",
           clustering_method = "ward.D2", 
           cluster_rows = T,
           cluster_cols = T) {
    
    clustering_order <- 
      lapply(unique(long_data[[group_col]]),
             function(group_) {
               
               long_data %>% 
                 rename(group = group_col) %>% 
                 filter(group == group_) %>% 
                 select(1:3) %>% 
                 spread(2, 3) %>% 
                 column_to_rownames(names(long_data)[1]) %>% 
                 get_clustering_order(distance_method = distance_method,
                                      clustering_method = clustering_method, 
                                      cluster_rows = cluster_rows,
                                      cluster_cols = cluster_cols)
             })
    
    clustering_order <- 
      clustering_order %>% 
      {list(row = map(., 
                      . %>% 
                        {.$row}),
            col = map(., 
                      . %>% 
                        {.$col}))} %>% 
      map(. %>% 
            unlist() %>% 
            unique()) 
    
    long_data %>% 
      rename(v1 = 1, 
             v2 = 2,
             val = 3) %>% 
      mutate(v1 = factor(v1, clustering_order$row),
             v2 = factor(v2, clustering_order$col)) %>% 
      set_names(names(long_data))
  }




do_cluster_rows <- 
  function(wide_data, 
           distance_method = "euclidean",
           clustering_method = "ward.D2") {
    
    wide_data %>% 
      dist(method = distance_method) %>% 
      hclust(method = clustering_method) 
    
  }

do_cor_clustering <- 
  function(wide_data, 
           cor_use = "everything",
           cor_method = "spearman",
           clustering_method = "average") {
    
    wide_data %>% 
      cor(method = cor_method,
          use = cor_use) %>% 
      {1 - .} %>% 
      as.dist() %>%
      hclust(method = clustering_method)
  }

get_clustering_labels <- 
  function(clust) {
    with(clust, labels[order])
  }

get_dendrogram_data <- 
  function(clust) {
    
    dendrodata <-
      ggdendro::dendro_data(clust)
    
    dendrodata$segments %>% 
      as_tibble() %>% 
      left_join(dendrodata$labels %>% 
                  as_tibble(),
                by = c("xend" = "x", "yend" = "y"))
  }



# ------ Retinagram functions -----

rotate_coords <- 
  function(x, y, rotate_angle, rotate_center = c(0, 0), radius_just = 0) {
    
    # Center data
    rotdata <- 
      tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             # Angle
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)) + rotate_angle,
             
             # Adjust hypotenuse
             hyp = hyp + radius_just,
             
             # New coordinates
             x = cos(angle) * hyp,
             y = sin(angle) * hyp,
             
             # Recenter coordinates
             x = x - rotate_center[1],
             y = y - rotate_center[2])
    
    
    rotdata
  }

shrink_rotation_coords <- 
  function(x, y, shrink_angle, rotate_center = c(0, 0)) {
    
    
    shrink_factor <- 
      1 - shrink_angle / (2 * pi)
    
    # Center data
    rotdata <- 
      tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)),
             
             # Shrink angle
             # angle = case_when(x == 0 & y == 0 ~ 0, 
             #                   quadrant %in% c(1,4) ~ angle + shrink_angle,
             #                   quadrant %in% c(2,3) ~ angle - shrink_angle),
             
             angle = angle * shrink_factor,
             
             # New coordinates
             x = cos(angle) * hyp,
             y = sin(angle) * hyp,
             
             # Recenter coordinates
             x = x - rotate_center[1],
             y = y - rotate_center[2])
    
    
    rotdata 
  }


calculate_retina_cut_angle <- 
  function(clust) {
    require(ggraph)
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    g <-
      ggraph(dendrogram, layout = 'dendrogram', circular = T)
    
    g_edgepoints <- 
      g$data %>% 
      as_tibble() %>% 
      filter(height == 0) %>% 
      left_join(cutree(clust, k = 2) %>% 
                  enframe("label", 
                          "cluster"),
                by = "label") %>% 
      mutate(angle = calculate_coord_angle(x, y))
    
    expand_grid(node1 = g_edgepoints$.ggraph.index,
                node2 = g_edgepoints$.ggraph.index) %>% 
      left_join(g_edgepoints %>% 
                  select(node1 = .ggraph.index,
                         angle1 = angle, 
                         cluster1 = cluster),
                by = "node1") %>% 
      left_join(g_edgepoints %>% 
                  select(node2 = .ggraph.index,
                         angle2 = angle, 
                         cluster2 = cluster),
                by = "node2") %>% 
      filter(cluster1 == 1,
             cluster2 == 2) %>% 
      group_by_all() %>% 
      mutate(dist = c(angle1 - angle2,
                      (angle1 - 2 * pi) - angle2,
                      angle1 - (angle2 - 2 * pi),
                      (angle1 - 2 * pi) - (angle2 - 2 * pi)) %>% 
               abs() %>% 
               min()) %>% 
      ungroup() %>% 
      arrange(dist) %>% 
      slice(1:2) %>% 
      mutate(cut_angle = (angle1 + angle2) / 2) %>% 
      pull(cut_angle) %>% 
      {2 * pi - .}
  }


calculate_coord_angle <- 
  function(x, y, rotate_center = c(0, 0)) {
    tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp))) %>% 
      pull(angle)
  }

circular_dendrogram_retinastyle <-
  function(clust, color_pal, 
           preserve_height = F,
           scale_expansion = c(0.25, 0.25), 
           text_size = 3, 
           width_range = c(1.5, 6), 
           arc_strength = 0.8, 
           default_color = "gray80", 
           rotate_angle = 0, 
           shrink_angle = 0,
           flip_text = F, 
           text_vjust = 0.5,
           text_hnudge = 0,
           elbow = F) {
    require(ggraph)
    require(igraph)
    require(viridis)
    require(tidyverse)
    require(magrittr)
    
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    
    if(preserve_height) {
      g <-
        ggraph(dendrogram, layout = 'dendrogram',
               height = height, circular = T)
    } else {
      g <-
        ggraph(dendrogram, layout = 'dendrogram',
               circular = T)
    }
    
    edge_data_temp <- 
      get_edges()(g$data) %>%
      as_tibble() %>% 
      mutate(row_number = row_number())
    
    
    edge_data_coords <- 
      edge_data_temp %>% 
      select(row_number, x, y, xend, yend) %>% 
      gather(dim, value, -row_number) %>% 
      mutate(dimtype = case_when(dim %in% c("xend", "yend") ~ "end",
                                 T ~ "start"),
             dim = gsub("end", "", dim)) %>% 
      spread(dim, value) 
    
    
    edge_data_coords_rotated <- 
      rotate_coords(edge_data_coords$x,
                    edge_data_coords$y, 
                    rotate_angle = rotate_angle, 
                    rotate_center = c(0, 0)) %>% 
      bind_cols(edge_data_coords %>% 
                  select(row_number, dimtype))
    
    edge_data_coords_shrunk <- 
      shrink_rotation_coords(edge_data_coords_rotated$x,
                             edge_data_coords_rotated$y, 
                             shrink_angle = shrink_angle, 
                             rotate_center = c(0, 0)) %>% 
      bind_cols(edge_data_coords %>% 
                  select(row_number, dimtype))
    
    edge_data_transformed <- 
      edge_data_coords_shrunk %>% 
      gather(dim, value, x, y) %>% 
      mutate(dim = ifelse(dimtype == "end", 
                          paste0(dim, dimtype), 
                          dim)) %>% 
      select(dim, row_number, value) %>% 
      spread(dim, value)
    
    
    
    edge_data <- 
      edge_data_temp %>% 
      select(-x, -y, -xend, -yend) %>% 
      left_join(edge_data_transformed,
                by = "row_number") %>% 
      left_join(color_pal %>%
                  enframe("label", "color"),
                by = c("node2.label" = "label")) %>%
      mutate(radius = xend^2 + yend^2) %>%
      arrange(-radius) %>% 
      mutate(edge_id = as.character(row_number()),
             rank_radius = unclass(factor(-radius))) 
    
    edge_id_colors <- 
      edge_data %>% 
      filter(!is.na(color)) %$%
      set_names(color, edge_id)
    
    for(rank_rad in 2:max(edge_data$rank_radius)) {
      edge_id_colors_new <- 
        left_join(edge_data %>%
                    select(edge_id, radius, xend, yend, rank_radius) %>%
                    filter(rank_radius == rank_rad),
                  edge_data %>%
                    select(edge_id, radius, x, y, rank_radius) %>%
                    filter(rank_radius < rank_rad),
                  by = c("xend" = "x", "yend" = "y")) %>%
        left_join(enframe(edge_id_colors),
                  by = c("edge_id.y" = "name")) %>%
        group_by(edge_id.x) %>% 
        summarise(color = ifelse(n_distinct(value) == 1 & any(value != default_color), 
                                 as.character(unique(value)),
                                 default_color)) %$%
        set_names(color, edge_id.x)
      edge_id_colors <- 
        c(edge_id_colors, edge_id_colors_new)
    }
    
    
    node_data_rotated <-
      g$data %>% 
      as_tibble() %$%
      rotate_coords(x,
                    y, 
                    rotate_angle = rotate_angle,
                    rotate_center = c(0, 0))  %>% 
      
      bind_cols(g$data %>% 
                  select(-x, -y))
    
    node_data_shrunk <- 
      node_data_rotated %$%
      shrink_rotation_coords(x,
                             y, 
                             shrink_angle = shrink_angle, 
                             rotate_center = c(0, 0)) %>% 
      bind_cols(g$data %>% 
                  select(-x, -y))
    
    node_data <- 
      node_data_shrunk
    
    
    text_degree <-
      ifelse(flip_text,
             180,
             0)
    
    
    label_data <-
      node_data %>%
      filter(label != "") %>%
      bind_cols(rotate_coords(.$x,
                              .$y,
                              rotate_angle = 0,
                              rotate_center = c(0, 0),
                              radius_just = text_hnudge) %>% 
                  set_colnames(paste(colnames(.), "adj", sep = "_"))) %>% 
      
      mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi + text_degree,
                                x < 0 ~ 360 - asin(y) * 180 / pi + text_degree)) %>%
      left_join(color_pal %>%
                  enframe("label", "color"),
                by = "label") 
    
    
    
    g <- 
      g +
      scale_edge_width(range = width_range, limits = c(0, 1)) +
      
      
      
      scale_edge_color_manual(values = edge_id_colors)  +
      label_data %>%
      {geom_node_text(data = .,
                      aes(x = x_adj,
                          y = y_adj,
                          label = label),
                      angle = .$degree,
                      hjust = case_when(.$x < 0 & !flip_text ~ 1,
                                        .$x >= 0 & !flip_text ~ 0,
                                        .$x < 0 & flip_text ~ 0,
                                        .$x >= 0 & flip_text ~ 1),
                      vjust = text_vjust,
                      size = text_size)}  +
      scale_x_continuous(expand = expansion(scale_expansion)) +
      scale_y_continuous(expand = expansion(scale_expansion)) +
      
      coord_fixed() +
      theme_void()
    
    if(elbow) {
      g + 
        geom_edge_elbow(data = edge_data,
                        aes(edge_color = edge_id,
                            edge_width = 1 - sqrt(xend^2 + yend^2)),
                        show.legend = F, 
                        lineend="round")
    } else {
      g +
        geom_edge_diagonal(data = edge_data,
                           aes(edge_color = edge_id,
                               edge_width = 1 - sqrt(xend^2 + yend^2)),
                           strength = arc_strength,
                           show.legend = F, 
                           lineend="round") 
    }
  }

# Deprecated:
circular_dendrogram_retinastyle_4 <-
  function(clust, color_mapping, label_col, color_col, 
           preserve_height = F,
           scale_expansion = c(0.25, 0.25), 
           text_size = 3, 
           width_range = c(1.5, 6), 
           arc_strength = 0.8, 
           default_color = "gray80", 
           rotate_angle = 0, 
           shrink_angle = 0,
           flip_text = F, 
           text_vjust = 0.5,
           elbow = F) {
    require(ggraph)
    require(igraph)
    require(viridis)
    require(tidyverse)
    require(magrittr)
    
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    # color_mapping <- 
    #   celltype_pal %>% 
    #   enframe()
    # label_col = "name" 
    # color_col = "value"
    # scale_expansion = c(0.25, 0.25) 
    # text_size = 3 
    # width_range = c(1.5, 6)
    # arc_strength = 0.5 
    # default_color = "gray80" 
    # rotate_circle = 3.7 
    # shrink_circle = 2
    # flip_text = F 
    # text_vjust = 0.3
    
    if(preserve_height) {
      g <-
        ggraph(dendrogram, layout = 'dendrogram',
               height = height, circular = T)
    } else {
      g <-
        ggraph(dendrogram, layout = 'dendrogram',
               circular = T)
    }
    
    edge_data_temp <- 
      get_edges()(g$data) %>%
      as_tibble() %>% 
      mutate(row_number = row_number())
    
    
    edge_data_coords <- 
      edge_data_temp %>% 
      select(row_number, x, y, xend, yend) %>% 
      gather(dim, value, -row_number) %>% 
      mutate(dimtype = case_when(dim %in% c("xend", "yend") ~ "end",
                                 T ~ "start"),
             dim = gsub("end", "", dim)) %>% 
      spread(dim, value) 
    
    
    edge_data_coords_rotated <- 
      rotate_coords(edge_data_coords$x,
                    edge_data_coords$y, 
                    rotate_angle = rotate_angle, 
                    rotate_center = c(0, 0)) %>% 
      bind_cols(edge_data_coords %>% 
                  select(row_number, dimtype))
    
    edge_data_coords_shrunk <- 
      shrink_rotation_coords(edge_data_coords_rotated$x,
                             edge_data_coords_rotated$y, 
                             shrink_angle = shrink_angle, 
                             rotate_center = c(0, 0)) %>% 
      bind_cols(edge_data_coords %>% 
                  select(row_number, dimtype))
    
    edge_data_transformed <- 
      edge_data_coords_shrunk %>% 
      gather(dim, value, x, y) %>% 
      mutate(dim = ifelse(dimtype == "end", 
                          paste0(dim, dimtype), 
                          dim)) %>% 
      select(dim, row_number, value) %>% 
      spread(dim, value)
    
    
    
    edge_data <- 
      edge_data_temp %>% 
      select(-x, -y, -xend, -yend) %>% 
      left_join(edge_data_transformed,
                by = "row_number") %>% 
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = c("node2.label" = "label")) %>%
      mutate(radius = xend^2 + yend^2) %>%
      arrange(-radius) %>% 
      mutate(edge_id = as.character(row_number()),
             rank_radius = unclass(factor(-radius))) 
    
    edge_id_colors <- 
      edge_data %>% 
      filter(!is.na(color)) %$%
      set_names(color, edge_id)
    
    for(rank_rad in 2:max(edge_data$rank_radius)) {
      edge_id_colors_new <- 
        left_join(edge_data %>%
                    select(edge_id, radius, xend, yend, rank_radius) %>%
                    filter(rank_radius == rank_rad),
                  edge_data %>%
                    select(edge_id, radius, x, y, rank_radius) %>%
                    filter(rank_radius < rank_rad),
                  by = c("xend" = "x", "yend" = "y")) %>%
        left_join(enframe(edge_id_colors),
                  by = c("edge_id.y" = "name")) %>%
        group_by(edge_id.x) %>% 
        summarise(color = ifelse(n_distinct(value) == 1 & any(value != default_color), 
                                 as.character(unique(value)),
                                 default_color)) %$%
        set_names(color, edge_id.x)
      edge_id_colors <- 
        c(edge_id_colors, edge_id_colors_new)
    }
    
    
    node_data_rotated <-
      g$data %>% 
      as_tibble() %$%
      rotate_coords(x,
                    y, 
                    rotate_angle = rotate_angle, 
                    rotate_center = c(0, 0))  %>% 
      
      bind_cols(g$data %>% 
                  select(-x, -y))
    
    node_data_shrunk <- 
      node_data_rotated %$%
      shrink_rotation_coords(x,
                             y, 
                             shrink_angle = shrink_angle, 
                             rotate_center = c(0, 0)) %>% 
      bind_cols(g$data %>% 
                  select(-x, -y))
    
    node_data <- 
      node_data_shrunk
    
    
    text_degree <-
      ifelse(flip_text,
             180,
             0)
    
    g <- 
      g +
      scale_edge_width(range = width_range, limits = c(0, 1)) +
      
      
      
      scale_edge_color_manual(values = edge_id_colors)  +
      node_data %>%
      filter(label != "") %>%
      mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi + text_degree,
                                x < 0 ~ 360 - asin(y) * 180 / pi + text_degree)) %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = "label") %>%
      {geom_node_text(data = .,
                      aes(label = label),
                      angle = .$degree,
                      hjust = case_when(.$x < 0 & !flip_text ~ 1,
                                        .$x >= 0 & !flip_text ~ 0,
                                        .$x < 0 & flip_text ~ 0,
                                        .$x >= 0 & flip_text ~ 1),
                      vjust = text_vjust,
                      size = text_size)}  +
      scale_x_continuous(expand = expansion(scale_expansion)) +
      scale_y_continuous(expand = expansion(scale_expansion)) +
      
      coord_fixed() +
      theme_void()
    
    if(elbow) {
      g + 
        geom_edge_elbow(data = edge_data,
                        aes(edge_color = edge_id,
                            edge_width = 1 - sqrt(xend^2 + yend^2)),
                        show.legend = F, 
                        lineend="round")
    } else {
      g +
        geom_edge_diagonal(data = edge_data,
                           aes(edge_color = edge_id,
                               edge_width = 1 - sqrt(xend^2 + yend^2)),
                           strength = arc_strength,
                           show.legend = F, 
                           lineend="round") 
    }
  }

# ------ Specificity functions -----

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

hpa_gene_classification <- 
  function(data, expression_col, tissue_col, gene_col, enr_fold, max_group_n, det_lim = 1) {
    data_ <- 
      data %>% 
      select(gene = gene_col,
             expression = expression_col,
             tissue = tissue_col) %>% 
      mutate(expression = round(expression, 4)) 
    
    if(any(is.na(data_$expression))) stop("NAs in expression column")
    if(any(is.na(data_$gene))) stop("NAs in gene column")
    if(any(is.na(data_$tissue))) stop("NAs in tissue column")
    
    n_groups <- length(unique(data_$tissue))
    
    gene_class_info <- 
      data_ %>%
      group_by(gene) %>%
      summarise(
        
        # Gene expression distribution metrics
        mean_exp = mean(expression, na.rm = T),
        min_exp = min(expression, na.rm = T),
        max_exp = max(expression, na.rm = T), 
        max_2nd = sort(expression)[length(expression)-1],
        
        # Expression frequency metrics
        n_exp = length(which(expression >= det_lim)),
        frac_exp = n_exp/length(expression[!is.na(expression)])*100,
        
        # Limit of enhancement metrics
        lim = max_exp/enr_fold, 
        
        exps_over_lim = list(expression[which(expression >= lim & expression >= det_lim)]),
        n_over = length(exps_over_lim[[1]]), 
        mean_over = mean(exps_over_lim[[1]]),
        min_over = ifelse(n_over == 0, NA,
                          min(exps_over_lim[[1]])),
        
        max_under_lim = max(expression[which(expression < min_over)], det_lim*0.1),
        
        
        exps_enhanced = list(which(expression/mean_exp >= enr_fold & expression >= det_lim)),
        
        
        
        
        # Expression patterns
        enrichment_group = paste(sort(tissue[which(expression >= lim & expression >= det_lim)]), collapse=";"),
        
        n_enriched = length(tissue[which(expression >= lim & expression >= det_lim)]),
        n_enhanced = length(exps_enhanced[[1]]), 
        enhanced_in = paste(sort(tissue[exps_enhanced[[1]]]), collapse=";"),
        n_na = n_groups - length(expression),
        max_2nd_or_lim = max(max_2nd, det_lim*0.1),
        tissues_not_detected = paste(sort(tissue[which(expression < det_lim)]), collapse=";"),
        tissues_detected = paste(sort(tissue[which(expression >= det_lim)]), collapse=";")) 
    
    
    gene_categories <- 
      gene_class_info %>%
      
      mutate(
        spec_category = case_when(n_exp == 0 ~ "not detected", 
                                  
                                  # Genes with expression fold times more than anything else are tissue enriched
                                  max_exp/max_2nd_or_lim >= enr_fold ~ "tissue enriched", 
                                  
                                  # Genes with expression fold times more than other tissues in groups of max group_n - 1 are group enriched
                                  max_exp >= lim &
                                    n_over <= max_group_n & n_over > 1 &
                                    mean_over/max_under_lim >= enr_fold ~ "group enriched", 
                                  
                                  # Genes with expression in tissues fold times more than the mean are tissue enhance
                                  n_enhanced > 0 ~ "tissue enhanced", 
                                  
                                  # Genes expressed with low tissue specificity
                                  T ~ "low tissue specificity"), 
        
        
        dist_category = case_when(frac_exp == 100 ~ "detected in all",
                                  frac_exp >= 31 ~ "detected in many",
                                  n_exp > 1 ~ "detected in some",
                                  n_exp == 1 ~ "detected in single",
                                  n_exp == 0 ~ "not detected"),
        
        spec_score = case_when(spec_category == "tissue enriched" ~ max_exp/max_2nd_or_lim,
                               spec_category == "group enriched" ~ mean_over/max_under_lim, 
                               spec_category == "tissue enhanced" ~ max_exp/mean_exp)) 
    
    
    
    
    ##### Rename and format
    gene_categories %>%
      mutate(enriched_tissues = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ enrichment_group,
                                          spec_category == "tissue enhanced" ~ enhanced_in),
             n_enriched = case_when(spec_category %in% c("tissue enriched", "group enriched") ~ n_enriched,
                                    spec_category == "tissue enhanced" ~ n_enhanced)) %>%
      select(gene, 
             spec_category, 
             dist_category, 
             spec_score,
             n_expressed = n_exp, 
             fraction_expressed = frac_exp,
             max_exp = max_exp,
             enriched_tissues,
             n_enriched,
             n_na = n_na,
             tissues_not_detected,
             tissues_detected) 
    
    
    
  }	

# ------ Enrichment functions -----
perform_ORA <-
  function(gene_associations,
           database,
           universe,
           n_clusters = 5,
           minGSSize = 10,
           maxGSSize = Inf,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.2) {
    require(clusterProfiler)
    require(multidplyr)
    
    if(n_clusters != 1) {
      worker_cluster <- new_cluster(n = n_clusters)
      cluster_library(worker_cluster, c("dplyr",
                                        "tidyverse"))
      cluster_copy(worker_cluster, c("enricher",
                                     "universe",
                                     "database",
                                     "minGSSize",
                                     "maxGSSize",
                                     "pvalueCutoff",
                                     "qvalueCutoff" ))
      
      pre_out <- 
        gene_associations %>%
        group_by(partition) %>%
        partition(worker_cluster) 
    } else {
      pre_out <- 
        gene_associations %>%
        group_by(partition)
    }
    
    outdata <-
      pre_out %>% 
      do({
        g_data <- .
        pull(g_data, gene) %>%
          enricher(universe = universe,
                   TERM2GENE = database, 
                   minGSSize = minGSSize,
                   maxGSSize = maxGSSize,
                   pvalueCutoff = pvalueCutoff,
                   qvalueCutoff = qvalueCutoff) %>%
          as_tibble()
      }) %>%
      ungroup() %>%
      collect()
    
    if(n_clusters != 1) rm(worker_cluster)
    outdata
  }

# ------ Alluvial functions -----

multi_alluvial_plot <- 
  function(data, vars, chunk_levels, pal, color_by = c(1, 3, 3)) {
    
    selvars = vars
    
    if(!is.null(names(vars))) {
      vars = names(vars)
    }
    
    alluv_1 <-
      data %>%
      ungroup() %>%
      select(selvars) %>% 
      ungroup() %>%
      mutate(row_n = row_number()) %>%
      gather(bar, chunk, -row_n) %>%
      left_join(tibble(bar = vars, 
                       color_vars = color_by), 
                by = "bar") %>% 
      group_by(row_n) %>%
      mutate(chunk_color = chunk[match(vars[color_vars], bar)]) %>% 
      ungroup() %>%
      
      mutate(chunk = factor(chunk, levels = chunk_levels),
             bar = factor(bar, levels = vars)) %>%
      
      
      ggplot(aes(x = bar, stratum = chunk, alluvium = row_n,
                 y = 1)) +
      
      geom_flow(aes(fill = chunk_color), 
                show.legend = F) +
      geom_stratum(aes(fill = chunk), 
                   show.legend = F, color = NA) +
      
      scale_x_discrete(expand = c(.1, .1), position = "top") +
      scale_fill_manual(values = pal) + 
      
      
      theme(axis.text.x = element_text(size = 18, face = "bold"),
            axis.text.y = element_blank(), 
            axis.ticks = element_blank(), 
            panel.background = element_blank(), 
            axis.title = element_blank())
    
    
    
    
    flow_data <-
      ggplot_build(alluv_1)$data[[1]] %>%
      as_tibble() %>%
      {
        if("side" %in% names(.)) {
          .
        } else{
          mutate(.,
                 side = case_when(flow == "from" ~ "start",
                                  flow == "to" ~ "end"))
        }}
    
    
    stratum_data <- 
      ggplot_build(alluv_1)$data[[2]]
    
    flow_data_labels <-
      flow_data %>% 
      as_tibble() %>% 
      
      select(x, stratum, group, side, ymin, ymax) %>% 
      pivot_wider(names_from = side, values_from = c(x, stratum, ymin, ymax)) %>%
      
      mutate_at(c("x_end", "ymax_end", "ymin_end", "x_start", "ymax_start", "ymin_start"), as.numeric) %>% 
      group_by(stratum_start, stratum_end, x_start, x_end) %>%
      summarise(y_end = (min(ymin_end) + max(ymax_end)) / 2, 
                y_start = (min(ymin_start) + max(ymax_start)) / 2, 
                size = max(ymax_start) - min(ymin_start))
    
    alluv_1 <- 
      alluv_1 +
      geom_text(data = flow_data_labels,
                aes(x = x_start + 1/6,
                    y = y_start, 
                    label = size), 
                inherit.aes = F, 
                size = 3, 
                hjust = 0) +
      geom_text(data = flow_data_labels,
                aes(x = x_end - 1/6,
                    y = y_end, 
                    label = size), 
                inherit.aes = F, 
                size = 3, 
                hjust = 1) +
      
      # Stratum label
      
      geom_text(data = stratum_data,
                aes(x = x, 
                    y = y,
                    label = paste(stratum, 
                                  ymax - ymin, sep = "\n")), 
                size = 4, 
                inherit.aes = F)
    
    
    alluv_1
  }



# ------ Color functions -----

ramp_color_bias <- 
  function(color1, color2, biases) {
    require(tidyverse)
    
    colors_rgb <- 
      colorRamp(c(color1, color2))(biases)
    
    
    rgb(colors_rgb[,1],
        colors_rgb[,2], 
        colors_rgb[,3], 
        maxColorValue = 255) 
  }

ramp_pal_bias <- 
  function(colors, biases) {
    require(tidyverse)
    
    colors_rgb <- 
      colorRamp(colors)(biases)
    
    
    rgb(colors_rgb[,1],
        colors_rgb[,2], 
        colors_rgb[,3], 
        maxColorValue = 255) 
  }


lab_to_hex <- 
  function(l, a, b) {
    
    rgb_code <- 
      tibble(l, a, b) %>% 
      convertColor(from='Lab', to='sRGB')
    
    rgb(red = rgb_code[,1],
        green = rgb_code[,2], 
        blue = rgb_code[,3],
        maxColorValue = 1)
  }

mix_2_colors <- 
  function(color1, color2, mix = 0) {
    sapply(1:length(color1),
           function(i) {
             colorRamp(c(color1[i], color2[i]))(mix[i]) %>% 
               rgb(maxColorValue = 255)
           }) 
  }

place_in_color_space <- 
  function(l = 70, a, b) {
    
    lab_to_hex(l = l,
               a = scales::rescale(a, to = c(-100, 100)),
               b = scales::rescale(b, to = c(-100, 100)))
    
    
  }

spread_colors <- 
  function(color, n, colorrange = 0.5) {
    
    if(n == 1) {
      return(color)
    } else {
      colorRamp(c("white", color, "black"))(seq(0.5 - colorrange / 2,
                                                0.5 + colorrange / 2, 
                                                length.out = n)) %>% 
        as_tibble() %$% 
        rgb(V1, V2, V3, maxColorValue = 255) 
      
    }
    
  }

# ------ Cluster hull functions -----

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
           plot_range = rep(c(min(c(V1, V2)),
                              max(c(V1, V2))),
                            2) * 1.05,
           poly_concavity = 1,
           poly_smoothing = 1,
           relative_bandwidth = 1/200) {
    require(tidyverse)
    require(magrittr)
    require(sf)
    require(sp)
    require(concaveman)
    
    # Compile indata
    cluster_data <- 
      tibble(V1, V2, 
             element_id,
             cluster = cluster_membership)
    
    # Calculate plot range
    plot_range_tb <- 
      set_names(plot_range,
                c("xmin", 
                  "xmax",
                  "ymin", 
                  "ymax")) %>% 
      enframe() %>% 
      spread(name, value)
    
    # Calculate diagonal length
    plot_diagonal <- 
      sqrt((plot_range_tb$xmax - plot_range_tb$xmin)^2 + 
             (plot_range_tb$ymax - plot_range_tb$ymin)^2)
    
    # Set bandwidth to fraction of diagonal
    plot_bandwidth <- 
      relative_bandwidth * plot_diagonal
    
    # Find subclusters
    subclusters <-
      cluster_data %>%
      group_by(cluster) %>%
      mutate(n_cluster_genes = n_distinct(element_id), 
             sub_cluster = data.frame(V1, V2) %>%
               fpc::dbscan(eps = plot_bandwidth) %$%
               cluster) %>%
      ungroup() %>% 
      group_by(cluster, sub_cluster) %>%
      mutate(n_sub_genes = n_distinct(element_id)) %>% 
      ungroup() 
    
    # Classify subclusters
    subclusters_classes <- 
      subclusters %>% 
      select(cluster, sub_cluster, n_cluster_genes, n_sub_genes) %>% 
      distinct() %>% 
      group_by(cluster) %>%
      mutate(n_sub_genes = ifelse(sub_cluster == 0, 
                                  0, 
                                  n_sub_genes),
             sub_type = case_when(sub_cluster == 0 ~ "outlier",
                                  n_sub_genes / n_cluster_genes < frac_lim ~ "outlier",
                                  rank(-n_sub_genes, 
                                       ties.method = "first") == 1 ~ "primary",
                                  T ~ "secondary")) %>% 
      select(cluster, sub_cluster, sub_type)
    
    subclusters_classed <- 
      subclusters %>% 
      left_join(subclusters_classes,
                by = c("cluster", "sub_cluster"))
    
    
    # Calculate plot density
    plot_density <- 
      subclusters_classed %>% 
      filter(sub_type != "outlier") %>%
      group_by(cluster, sub_cluster, sub_type) %>%
      do({
        get_density(.$V1, 
                    .$V2, 
                    h = plot_bandwidth, 
                    n = n, 
                    lims = plot_range) 
        
      }) %>% 
      ungroup() %>% 
      filter(z > 1e-200) %>% 
      group_by(cluster, sub_cluster) %>% 
      mutate(z = z / sum(z)) %>% 
      arrange(cluster, sub_cluster, -z) %>% 
      mutate(cum_z = cumsum(z)) %>% 
      ungroup()
    
    
    # Filter pixels such that 95% of density is included
    # Each point is then assigned to the cluster with highest density
    plot_density_filtered <- 
      plot_density %>% 
      filter(cum_z < cum_z_lim) %>%
      group_by(x, y) %>% 
      top_n(1, z) %>%
      slice(1) %>% 
      ungroup()
    
    
    # Calculate size of landmass
    plot_density_landmass <-
      plot_density_filtered %>%
      group_by(cluster, sub_cluster) %>%
      mutate(landmass = data.frame(x_coord, y_coord) %>%
               fpc::dbscan(eps = plot_bandwidth) %$%
               cluster) %>%
      group_by(cluster, sub_cluster, landmass) %>% 
      mutate(n_landmass_points = length(x)) %>%  
      ungroup() %>% 
      group_by(cluster) %>%
      mutate(n_total_points = length(x)) %>% 
      ungroup() 
    
    # Classify landmasses
    plot_density_landmass_classes <- 
      plot_density_landmass %>% 
      select(cluster, sub_cluster, landmass, n_landmass_points, n_total_points) %>% 
      distinct() %>%
      group_by(cluster, sub_cluster) %>%
      
      mutate(frac_landmass = n_landmass_points / n_total_points,
             landmass_type = case_when(rank(-n_landmass_points, 
                                            ties.method = "first") == 1 ~ "primary",
                                       T ~ "secondary")) %>% 
      ungroup() %>% 
      select(cluster, sub_cluster, landmass, landmass_type, frac_landmass)
    
    plot_density_landmass_classed <- 
      plot_density_landmass %>%
      left_join(plot_density_landmass_classes,
                by = c("cluster", 
                       "sub_cluster",
                       "landmass")) %>% 
      arrange(cluster)
    
    plot_density_mainland_filtered <-
      plot_density_landmass_classed %>% 
      filter(frac_landmass > frac_lim) %>%
      filter(cum_z < cum_z_lim) %>%
      group_by(x, y) %>% 
      top_n(1, z) %>%
      slice(1) %>% 
      ungroup()
    
    # Create polygons
    # Poly smoothing: How small distances should be further detailed - 
    # higher values --> less detailed
    # poly concavity: How convex polygons should be - 
    # higher values --> less detailed
    
    plot_data_hulls <- 
      plot_density_mainland_filtered %>% 
      
      group_by(cluster, sub_cluster, landmass, sub_type) %>% 
      do({
        st_as_sf(., coords=c('x_coord','y_coord')) %>%
          concaveman(concavity = poly_concavity, 
                     length_threshold = plot_bandwidth * poly_smoothing) %$%
          st_coordinates(polygons) %>% 
          as_tibble()
      }) %>% 
      ungroup() %>% 
      mutate(polygon_id = paste(cluster, sub_cluster, landmass, sep = "_"))
    
    
    
    # plot_data_hulls %>% 
    #   filter(sub_type == "primary") %>% 
    #   group_by(cluster) %>% 
    #   do({
    #     g_data <<- .
    #     
    #     g_data %>% 
    #       select(X, Y) %>% 
    #       # column_to_rownames("element_id") %>% 
    #       dist() %>% 
    #       as.matrix() %>% 
    #       as_tibble() %>% 
    #       colSums() %>% 
    #       enframe("element_id", "sumdist") %>% 
    #       arrange(sumdist) 
    #       
    #   })
    
    plot_density_center <- 
      plot_density %>% 
      filter(sub_type == "primary") %>% 
      group_by(cluster) %>% 
      top_n(1, z) %>% 
      slice(1) %>% 
      ungroup() %>% 
      select(cluster, x = x_coord, y = y_coord)
    
    
    
    list(hulls = plot_data_hulls,
         density = plot_density,
         landmass_pixels = plot_density_mainland_filtered,
         center_density = plot_density_center)
  }



# ------ Normalization -----

calc_tmm_normfactors <- 
  function (object, method = c("TMM", "quantile"), refColumn = NULL, 
            logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, 
            Acutoff = -1e+10, quantile = 0.75) {
    method <- match.arg(method)
    if (is.matrix(object)) {
      if (is.null(refColumn)) 
        refColumn <- 1
      data <- object
      libsize <- colSums(data)
    } else {
      stop("calcNormFactors() only operates on 'matrix' objects")
    }
    
    if(refColumn == "median") {
      ref <- 
        apply(data, MARGIN = 1, median)
    } else {
      ref <- data[, refColumn]
    }
    
    f <- switch(method, TMM = apply(data, 2, NOISeq:::.calcFactorWeighted, 
                                    ref = ref, logratioTrim = logratioTrim, 
                                    sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff), 
                quantile = NOISeq:::.calcFactorQuantile(data, libsize, q = quantile))
    f <- f/exp(mean(log(f)))
    return(f)
  }
