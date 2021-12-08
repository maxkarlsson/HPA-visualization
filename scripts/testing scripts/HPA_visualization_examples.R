
library(tidyverse)
source("scripts/HPA_visualization_functions.R")
c("tidyverse", 
  "uwot", 
  "pcaMethods")



long_data <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/hpa_celline_all_samples_gene_tpm_103.tsv")



celline_data <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/For gene clustering/proteinatlas_celline_hpa_103.tsv")


pca_scores <- 
  celline_data %>% 
  convert_to_widedata("celline", "ensg_id", "ntpm") %>% 
  remove_genes(rm_not_expressed = T, 
               rm_no_variance = T) %>% 
  scale_data(logp1_scale = F, 
             zscore_scale = T) %>% 
  do_pca() %>% 
  get_pca_scores(use_R2cum_PCselection = T,
                 R2cum_lim = 0.8) 

umap_data <- 
  pca_scores %>% 
  do_umap()

heatmap_data <- 
  celline_data %>% 
  filter(ensg_id %in% unique(ensg_id)[1:1000]) %>% 
  convert_to_widedata("celline", "ensg_id", "ntpm") %>% 
  remove_genes(rm_not_expressed = T, 
               rm_no_variance = T) %>% 
  scale_data(logp1_scale = T, 
             zscore_scale = T) %>% 
  cluster_wide_data()


long_data

umap_data %>% 
  left_join(meta_data) %>% 
  ggplot(aes(UMAP1, UMAP2, group = anno)) +
  geom_point() +
  geom_path() +
  coord_fixed() +
  theme_bw()

single_cell_consensus_data <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/single cell/consensus.tsv")

single_cell_consensus_data %>% 
  convert_to_widedata("cell_type_name", "ensg_id", "exp") %>% 
  do_cor_clustering() %>% 
  get_clustering_labels() %>% 
  paste0(collapse = "\n") %>% 
  cat
  

single_cell_data <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/single cell/cluster_data.tsv")

single_cell_meta <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/single cell/annotation.tsv")



single_wide <- 
  single_cell_data %>% 
  unite(id, assay_id, cluster_id) %>% 
  select(ensg_id, id, ntpm) %>% 
  convert_to_widedata("id", "ensg_id", "ntpm") 

single_UMAP <- 
  single_wide %>% 
  HPA_sample_UMAP()

""
library(ggrepel)
library(plotly)

plot <- 
  single_UMAP %>% 
  left_join(single_cell_meta %>% 
              unite(sample, assay_id, cluster_id)) %>% 
  mutate(label = paste(tissue_name, cell_type_name)) %>% 
  ggplot(aes(UMAP1, UMAP2, color = group_name, group = label)) +
  geom_point() +
  geom_text(data = . %>% 
              group_by(cell_type_name) %>% 
              summarise(UMAP1 = median(UMAP1),
                        UMAP2 = median(UMAP2)),
            aes(UMAP1, UMAP2, label = cell_type_name),
            inherit.aes = F, 
            size = 2) +
  # geom_path() +
  coord_fixed() +
  theme_bw()

ggplotly(plot)
##############

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
  extract_variable_features(log1p_data = T)



b <- 
  a %>% 
  # Scale data genewise, first log10(exp + 1), 
  # and then z-score:
  scale_data(logp1_scale = T, 
             zscore_scale = T) 


cc <- 
  b%>% 
  
  # Perform PCA:
  do_pca() %>% 
  
  # Get PCA scores, selecting PCs that constitute 80% cumulative r2:
  get_pca_scores(use_R2cum_PCselection = T,
                 R2cum_lim = 0.8)

a_data <- 
  tibble(neighs = 5:15) %>% 
  group_by(neighs) %>% 
  do({
    neigh <- .$neighs
    cc %>% 
      # Perform UMAP with default settings:
      do_umap(n_neighbors = neigh) %>% 
      
      # Add rownames as a column for saving:
      as_tibble(rownames = "sample")
  })

a_data %>% 
  ggplot(aes(UMAP1, UMAP2)) +
  geom_point() +
  facet_wrap(~neighs) +
  coord_fixed()



vardata <- 
  a %>% 
  log1p() %>% 
  {tibble(mean = apply(., 
                       MARGIN = 1, 
                       mean),
          sd = apply(., 
                     MARGIN = 1, 
                     sd))}


vardata_loess <- 
  loess(sd ~ mean, data = vardata)

vardata %>% 
  mutate(pred_sd = predict(vardata_loess, newdata = .),
         high = sd > pred_sd) %T>% 
  {group_by(., high) %>% 
      count %>% 
      print} %>% 
  ggplot(aes(mean, sd, color  = high)) +
  geom_point() +
  geom_smooth(color = "black")
