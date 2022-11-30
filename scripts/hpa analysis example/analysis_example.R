

setwd("scripts/hpa analysis example/")

library(tidyverse)
library(tidygraph)
library(ggraph)
library(geomtextpath)

source("../HPA_visualization_functions.R")

# ------ read data --------


sc_cluster_annotation <- read_tsv("data/single_cell_annotation_v3_103.tsv")
sc_class <- read_tsv("data/single_cell_category_v3_103.tsv", na = "NULL") %>% 
  rename(specificity_category = category_name)
sc_cluster_data <- read_tsv("data/single_cell_cluster_data_v3_103.tsv")
sc_consensus_data <- read_tsv("data/single_cell_aggregated_data_v3_103.tsv")



sc_consensus_data_wide <- 
  sc_consensus_data %>% 
  spread(cell_type_name, exp) %>% 
  column_to_rownames("ensg_id")

sc_cluster_data_wide <- 
  sc_cluster_data %>% 
  unite(id, assay_id, cluster_id) %>% 
  select(ensg_id, id, ntpm) %>% 
  spread(id, ntpm) %>% 
  column_to_rownames("ensg_id")


# ------ establish theme -------


spec_category_levels <- 
  c('tissue enriched',
    'group enriched', 
    'tissue enhanced', 
    'low tissue specificity',  
    'not detected',
    
    'Tissue enriched',
    'Group enriched', 
    'Tissue enhanced', 
    'Low tissue specificity',  
    'Not detected')

gene_category_pal <- 
  c("tissue enriched" = "#e41a1c",
    "group enriched" = "#FF9D00",
    "tissue enhanced" = "#984ea3",
    "low tissue specificity" = "grey40",
    
    "detected in all" = "#253494",
    "detected in many" = "#2c7fb8",
    "detected in some" = "#41b6c4",
    "detected in single" = "#a1dab4",
    
    "not detected" = "grey", 
    "not detected " = "grey", 
    
    "Tissue enriched" = "#e41a1c",
    "Group enriched" = "#FF9D00",
    "Tissue enhanced" = "#984ea3",
    "Low tissue specificity" = "grey40",
    
    "Detected in all" = "#253494",
    "Detected in many" = "#2c7fb8",
    "Detected in some" = "#41b6c4",
    "Detected in single" = "#a1dab4",
    
    "Not detected" = "grey", 
    "Not detected " = "grey")

consensus_colors <- 
  read_tsv("data/colors_consensus.tsv") %>% 
  select(sample, color) %>% 
  deframe()

# ------ pca + umap --------

plot_data <- 
  sc_consensus_data_wide %>% 
  remove_genes(rm_NA_genes = T, 
               rm_not_expressed = T, 
               not_expressed_lim = 1, 
               rm_no_variance = T) %>% 
  scale_data(logp1_scale = T, zscore_scale = T) 

plot_pca <- 
  plot_data %>% 
  do_pca(npcs = 20) %>% 
  get_pca_scores() 

plot_pca %>% 
  as_tibble(rownames = "celltype") %>% 
  ggplot(aes(PC1, PC2, color = celltype)) +
  geom_point(show.legend = F) +
  coord_fixed() +
  scale_color_manual(values = consensus_colors) +
  theme_bw()


plot_umap <- 
  plot_pca %>% 
  do_umap(n_neighbors = 10) %>% 
  as_tibble(rownames = "celltype") 


plot_umap %>% 
  ggplot(aes(UMAP1, UMAP2, color = celltype)) +
  geom_point(show.legend = F) +
  coord_fixed() +
  scale_color_manual(values = consensus_colors) +
  theme_bw()

# ------ classification using IT's version--------

plot_data <- 
  sc_class %>% 
  group_by(specificity_category) %>% 
  count() %>% 
  mutate(specificity_category = factor(specificity_category, spec_category_levels)) 

plot_data %>% 
  ggplot(aes(1, n, label = n, fill = specificity_category)) + 
  geom_col(show.legend = F) +
  geom_textpath(position = position_stack(vjust = 0.5), 
                aes(x = 1.5),
                angle = 90) +
  # facet_wrap(~dataset_id) +
  scale_fill_manual(values = gene_category_pal) +
  coord_polar("y") +
  theme_void()


plot_data %>% 
  ggplot(aes(1, n, label = n, fill = specificity_category)) + 
  geom_col(show.legend = F) +
  geom_textpath(position = position_stack(vjust = 0.5), 
                hjust = 1,
                aes(x = 1.5)) +
  scale_fill_manual(values = gene_category_pal) +
  scale_x_discrete(expand = expansion(c(0,1)))+
  coord_polar("y") +
  theme_void()

# Network
plot_data <-
  sc_class %>%
  mutate(specificity_category = factor(specificity_category, spec_category_levels),
         enhanced_tissues = str_replace_all(enhanced_tissues, 
                                            ", ",
                                            "¤") %>% 
           str_replace_all(",", ";") %>% 
           str_replace_all("¤", ", "),
         group_node = enhanced_tissues) %>% 
  filter(specificity_category %in% c("Tissue enriched",
                                     "Group enriched")) %>% 
  
  select(sample = enhanced_tissues, group_node, 
         specificity_category) %>% 
  group_by_all() %>% 
  count() %>% 
  ungroup() %>%
  separate_rows(sample, sep = ";") %>% 
  arrange(sample) %>% 
  filter(n > 2 | specificity_category == "Tissue enriched") %>%
  group_by(sample, specificity_category) %>% 
  mutate(edge_rank = rank(-n, ties.method = "min")) %>% 
  # filter(sample == "brain") %>% 
  group_by(group_node) %>% 
  mutate(high_rank = any(edge_rank <= 2)) %>% 
  ungroup() %>% 
  filter(high_rank) %>% 
  unite(group_node, group_node, specificity_category, n, remove = F) 

plot_data %>% 
  select(sample, group_node) %>% 
  as_tbl_graph(directed = F) %>%
  left_join(as_tibble(.) %>% 
              separate(name, 
                       into = c("node_name", "node_type", "n"),
                       sep = "_", remove = F)) %>% 
  mutate(node_type = ifelse(is.na(node_type),
                            "sample", 
                            node_type),
         node_color = ifelse(node_type == "sample",
                             name, 
                             node_type),
         n = as.numeric(n)) %>% 
  ggraph(layout = "nicely") +
  geom_edge_link() +
  geom_node_point(data = . %>% 
                    filter(node_type == "sample"),
                  aes(color = node_color),
                  show.legend = F,
                  size = 6) +
  geom_node_point(data = . %>% 
                    filter(node_type != "sample"),
                  aes(color = node_color,
                      size = sqrt(n)),
                  show.legend = F) +
  geom_node_text(data = . %>% 
                   filter(node_type == "sample"),
                 aes(label = str_wrap(name, width = 10)),
                 lineheight = 0.8,
                 size = 2) +
  geom_node_text(data = . %>% 
                   filter(node_type != "sample"),
                 aes(label = n),
                 color = "white",
                 size = 2) +
  scale_color_manual(values = c(consensus_colors, gene_category_pal)) +
  scale_size_continuous(range = c(0, 10),
                        limit = c(0, max(sqrt(plot_data$n)))) +
  
  coord_fixed() +
  theme_void()


#  network cytoscape input
net_data <- 
  plot_data %>% 
  select(sample, group_node) %>% 
  as_tbl_graph(directed = F) %>%
  left_join(as_tibble(.) %>% 
              separate(name, 
                       into = c("node_name", "node_type", "n"),
                       sep = "_", remove = F)) %>% 
  mutate(node_type = ifelse(is.na(node_type),
                            "sample", 
                            node_type),
         node_color = ifelse(node_type == "sample",
                             name, 
                             node_type),
         n = as.numeric(n)) 


node_data <-   
  net_data %>% 
  as_tibble() %>% 
  mutate(node_id = row_number()) %>% 
  select(node_id, everything()) %>% 
  mutate(node_label = ifelse(node_type == "sample",
                             name, 
                             n)) %>% 
  mutate(node_color = c(consensus_colors, 
                        gene_category_pal)[match(node_color, 
                                                 names(c(consensus_colors, 
                                                         gene_category_pal)))],
         node_color2 = ifelse(node_type == "sample", 
                              "#BFBFBF",
                              node_color),
         node_size = sqrt(n))


edge_data <- 
  net_data %>% 
  activate(edges) %>% 
  as_tibble() 

# Save nodes and edges to use as input for cytoscape:
# write_csv(node_data, "sc class net nodes.csv")
# write_csv(edge_data, "sc class net edges.csv")




# ------ classification by manual calculation version--------

manual_class <- 
  sc_consensus_data %>%
  hpa_gene_classification(expression_col = "exp",
                          tissue_col = "cell_type_name", 
                          gene_col = "ensg_id",
                          enr_fold = 4, 
                          max_group_n = 10)
  
plot_data <- 
  manual_class %>% 
  group_by(spec_category) %>% 
  count() %>% 
  mutate(spec_category = factor(spec_category, spec_category_levels)) 

plot_data %>% 
  ggplot(aes(1, n, label = n, fill = spec_category)) + 
  geom_col(show.legend = F) +
  geom_textpath(position = position_stack(vjust = 0.5), 
                hjust = 1,
                aes(x = 1.5)) +
  scale_fill_manual(values = gene_category_pal) +
  scale_x_discrete(expand = expansion(c(0,1)))+
  coord_polar("y") +
  theme_void()

# ------ tau --------

detected_genes <- 
  sc_class %>% 
  filter(specificity_category != "Not detected")

sc_tau_score <- 
  sc_consensus_data_wide[detected_genes$ensg_id,] %>% 
  {log10(. + 1)} %>% 
  calculate_tau_score()

sc_tau_score %>% 
  ggplot(aes(tau_score)) +
  geom_density(color = NA,
               fill = "darkcyan",
               alpha = 0.5) +
  theme_bw()


sc_class %>% 
  left_join(sc_tau_score,
            by = c("ensg_id" = "gene")) %>% 
  mutate(specificity_category = factor(specificity_category, spec_category_levels)) %>% 
  ggplot(aes(tau_score, specificity_category, fill = specificity_category)) +
  geom_violin(draw_quantiles = 0.5,
              show.legend = F) +
  scale_fill_manual(values = gene_category_pal) +
  theme_bw()
  
# ------ compare classifications --------

# Here we simulate a classification that we can compare the real one with
class_probs <- 
  table(sc_class$specificity_category)/nrow(sc_class)

sc_class_2 <- 
  sc_class %>% 
  mutate(scramble = sample(c(T, F), nrow(.), prob = c(0.1, 0.9), replace = T)) %>% 
  group_by(scramble) %>% 
  mutate(specificity_category = ifelse(scramble,
                                       sample(names(class_probs), nrow(.),
                                              prob = class_probs, replace = T),
                                       specificity_category)) %>% 
  ungroup()


sc_class %>% 
  select(1,2) %>% 
  left_join(sc_class_2 %>% 
              select(1,2),
            by = "ensg_id", 
            suffix = c("_1", "_2")) %>% 
  multi_alluvial_plot(vars = c("v22" = "specificity_category_1", 
                               "'v23'" = "specificity_category_2"), 
                      chunk_levels = c('Tissue enriched', 'Group enriched', 
                                       'Tissue enhanced', 'Low tissue specificity', 
                                       'Not detected'), 
                      pal = c(gene_category_pal), 
                      color_by = c(1, 1)) 



# Here we compare IT's classification with the manually calculated one


sc_class %>% 
  select(1,2) %>% 
  left_join(manual_class %>% 
              select(1,2) %>% 
              mutate(spec_category = str_to_sentence(spec_category)),
            by = c("ensg_id" = "gene"), 
            suffix = c("_1", "_2")) %>% 
  multi_alluvial_plot(vars = c("IT" = "specificity_category", 
                               "Manual" = "spec_category"), 
                      chunk_levels = c('Tissue enriched', 'Group enriched', 
                                       'Tissue enhanced', 'Low tissue specificity', 
                                       'Not detected'), 
                      pal = c(gene_category_pal), 
                      color_by = c(1, 1)) +
  ggtitle("The categories are exactly the same",
          "If they weren't it would be very troubling")


# ------ tmm normalization --------

# Here we will perform tmm adjustment from ptpm data: 

# Data input is a wide format matrix so let's create that format: 
wide_data <- 
  sc_cluster_data %>% 
  select(ensg_id, assay_id, cluster_id, ptpm) %>% 
  left_join(select(., assay_id, cluster_id) %>% 
              distinct() %>% 
              unite(id, assay_id, cluster_id, sep = ";", remove = F)) %>% 
  select(ensg_id, id, ptpm) %>% 
  spread(id, ptpm) %>% 
  column_to_rownames("ensg_id") %>% 
  as.matrix()

# Calculate and apply sample wise tmm factors
tmm_factors <- 
  wide_data  %>% 
  calc_tmm_normfactors(method = "TMM", 
                       
                       # Use median column as reference distribution:
                       refColumn = "median",
                       
                       # Trim parameters:
                       logratioTrim = 0.3, 
                       sumTrim = 0.3, 
                       
                       # Weighting should be done only if count data:
                       doWeighting = F) %>% 
  enframe("sample", "tmm_factor")

wide_data_tmm <- t(t(wide_data) / tmm_factors$tmm_factor)


sc_cluster_data_2 <- 
  wide_data_tmm %>% 
  as_tibble(rownames = "ensg_id") %>% 
  gather(id, ntpm, -1) %>% 
  left_join(select(., id) %>% 
              distinct() %>% 
              separate(id, into = c("assay_id", "cluster_id"), sep = ";", remove = F))  %>% 
  select(ensg_id, assay_id, cluster_id, ntpm)


joined_data <- 
  sc_cluster_data %>% 
  mutate(assay_id = as.character(assay_id),
         cluster_id = as.character(cluster_id)) %>% 
  left_join(sc_cluster_data_2,
            by = c("ensg_id", "assay_id", "cluster_id"),
            suffix = c("_1", "_2"))
            

# We see the largest difference between IT's calculated nTPM and our manually calculated nTPM is virtually 0
# That's really good!
joined_data %>% 
  mutate(diff = ntpm_1 - ntpm_2) %>% 
  arrange(-abs(diff))
  






         