suppressMessages(library(RMariaDB))
suppressMessages(library(tidyverse))

# Read the script at its location
source("scripts/HPA_visualization_functions.R")

args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (SQL query).n", call. = FALSE)
}

## Connect to MySQL
rmysql.settingsfile <- "../.config.txt"

long_data <- 
  dbConnect(RMariaDB::MariaDB(),
            default.file = rmysql.settingsfile, 
            group = 'atlas') %>% 
  dbGetQuery(con, args[1]) %>% 
  as_tibble()

heatmap_data <- 
  long_data %>% 
  # Convert data to wide format:
  convert_to_widedata("sample", "gene_id", "exp") %>% 
  # Scale data to max value per gene
  scale_data(max_scale = T) %>% 
  # Cluster data
  cluster_wide_data()
  
## OUTPUT
cat(format_tsv(heatmap_data))
