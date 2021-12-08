
library(tidyverse)
library(pcaMethods)
library(uwot)

#######################################################
n_points = 1000
n_samples = 100
n_groups = 2
centers <-
  rnorm(n_points * n_groups) %>% 
  enframe("i", "mean") %>% 
    mutate(group = rep(1:n_groups, each = n_points),
           gene = rep(1:n_points, n_groups))

sds <-
  rnorm(n_points * n_groups) %>% 
  abs() %>% 
  enframe("i", "sd") %>% 
  mutate(group = rep(1:n_groups, each = n_points),
         gene = rep(1:n_points, n_groups))

settings <- 
  left_join(centers, 
            sds)

indata <- 
  settings %>% 
  group_by(group, gene) %>% 
  do({
    g_data <- .
    
    tibble(sample = paste(g_data$group, 1:n_samples, sep = "_"),
           value = rnorm(n_samples, 
                         mean = g_data$mean,
                         sd = g_data$sd))
  }) %>% 
  ungroup() %>% 
  mutate(value = round(value, 2),
         sample = paste0("sample", sample),
         group = paste0("group", group))



indata_wide <- 
  indata %>% 
  select(sample, gene, value) %>% 
  spread(sample, value)

indata_meta <- 
  indata %>% 
  select(sample, group) %>% 
  distinct()

#######################################################

HPA_UMAP <- 
  function(wide_data, 
           seed = 1,
           filter_zero_sd = F,
           n_epochs = 500,
           n_neighbors = 15) {
    
    set.seed(seed)
    
    if(filter_zero_sd) {
      row_sd <- 
        wide_data %>% 
        column_to_rownames(colnames(wide_data)[1]) %>% 
        apply(MARGIN = 1, 
              sd)
      wide_data <- 
        wide_data %>% 
        slice(which(row_sd != 0))
    }
    
    pca_res <- 
      wide_data %>% 
      column_to_rownames(colnames(wide_data)[1]) %>% 
      t() %>% 
      scale() %>% 
      pca(nPcs = dim(.)[1])
    
    
    pc_lim <- 
      which(pca_res@R2cum > 0.8)[1]
    
    pc_lim_sd <- 
      rev(which(pca_res@sDev > 1))[1]
    
    pca_res@scores[,1:pc_lim] %>% 
      umap(n_neighbors = n_neighbors,
           n_epochs = n_epochs) %>% 
      as_tibble() %>% 
      set_names(paste0("UMAP", 1:ncol(.))) %>% 
      mutate(sample = rownames(pca_res@scores)) %>% 
      select(sample, everything())
  }

umap_res <- HPA_UMAP(indata_wide)

umap_res %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point() +
  coord_fixed()

celline_data <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/proteinatlas_celline_hpa_103.tsv")

celline_wide <- 
  celline_data %>% 
  select(ensg_id, celline, ntpm) %>% 
  spread(celline, ntpm)


umap_res_cl <- HPA_UMAP(celline_wide, 
                        filter_zero_sd = T,
                        n_neighbors = 30)




umap_res_cl %>% 
  ggplot(aes(UMAP1, UMAP2, label = sample)) +
  geom_point() + 
  geom_text()

