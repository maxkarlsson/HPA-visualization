
a <- 

  long_data %>% 
  # Convert data to wide format:
  convert_to_widedata("sample", "gene_id", "exp") %>% 
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

b <- 
  long_data %>% 
  # Convert data to wide format:
  convert_to_widedata("sample", "gene_id", "exp") %>% 
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
  t() %>% 
  # Perform UMAP with default settings:
  do_umap(n_neighbors = round(sqrt(dim(.)[1]))) %>% 
  
  # Add rownames as a column for saving:
  as_tibble(rownames = "sample")





a %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point() +
  coord_fixed() +
  theme_bw()

b %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point() +
  coord_fixed() +
  theme_bw()

celline_meta <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/hpa_celline_all_samples_gene_tpm_103.tsv") %>% 
  select(sample = 3, celline = 2) %>% 
  distinct()

p1 <-  
  a %>% 
  left_join(celline_meta) %>% 
  ggplot(aes(UMAP1, UMAP2, group = celline)) +
  geom_path() +
  geom_point() +
  geom_text(data = . %>% 
              group_by(celline) %>% 
              summarise(UMAP1 = median(UMAP1),
                        UMAP2 = median(UMAP2)),
            aes(UMAP1, UMAP2, label = celline),
            inherit.aes = F, 
            size = 2) +
  # geom_path() +
  coord_fixed() +
  theme_bw() +
  ggtitle("With PCA")

p2 <- 
  b %>% 
  left_join(celline_meta) %>% 
  ggplot(aes(UMAP1, UMAP2, group = celline)) +
  geom_path() +
  geom_point() +
  geom_text(data = . %>% 
              group_by(celline) %>% 
              summarise(UMAP1 = median(UMAP1),
                        UMAP2 = median(UMAP2)),
            aes(UMAP1, UMAP2, label = celline),
            inherit.aes = F, 
            size = 2) +
  # geom_path() +
  coord_fixed() +
  theme_bw() +
  ggtitle("No PCA")

library(patchwork)
p1 + p2



plot_data <- 
  a %>% 
  left_join(b,
            suffix = c("_a", "_b"),
            by = "sample") %>% 
  gather(variable, value, -1)


plot_data %>% 
  select(sample, var1 = variable) %>% 
  distinct() %>% 
  expand_grid(var2 = unique(.$var1)) %>% 
  # filter(grepl("UMAP1", var1),
  #        grepl("UMAP2", var2)) %>% 
  left_join(plot_data,
            by = c("sample",
                   "var1" = "variable")) %>% 
  left_join(plot_data,
            by = c("sample",
                   "var2" = "variable"),
            suffix = c("1", "2")) %>% 
  ggplot(aes(value1, value2)) +
  geom_point() +
  coord_fixed() +
  facet_grid(var2 ~ var1)

  
list("With PCA" = a, 
     "No PCA" = b) %>% 
  bind_rows(.id = "type") %>% 
  # mutate(UMAP1 = ifelse(type == "With PCA",
  #                       UMAP1 + 4,
  #                       UMAP1),
  #        UMAP2 = ifelse(type == "With PCA",
  #                       UMAP2 + 4,
  #                       UMAP2)) %>%
  ggplot(aes(UMAP1, UMAP2, group = sample, color = type)) +

  geom_path(alpha = 0.4, 
            color = "gray") +
  geom_point() +
  theme_bw()
  

list("With PCA" = a, 
     "No PCA" = b) %>% 
  bind_rows(.id = "type") %>% 
  left_join(celline_meta) %>% 
  group_by(type, celline) %>% 
  mutate(UMAP1 = (UMAP1 - mean(UMAP1)),
         UMAP2 = (UMAP2 - mean(UMAP2))) %>% 
  ggplot(aes(UMAP1, UMAP2, group = celline, color = type)) +
  
  geom_path(alpha = 0.4, 
            color = "gray") +
  ggalt::geom_encircle(aes(UMAP1, UMAP2, fill = type),
                       inherit.aes = F, alpha = 0.5, 
                       s_shape=1, expand=0) +
  facet_wrap(~type) +
  geom_point() +
  theme_bw()


  