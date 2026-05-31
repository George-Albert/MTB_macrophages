

#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(xlsx)
  library(ggpubr)
  library(shiny)
}
# cmdstanr::rebuild_cmdstan()
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### Set dir
main_dir <- getwd()
setwd(main_dir)
input_dir      <- "Analyses/Inputs"
output_dir     <- "Analyses/Outputs"
eqtl_input_dir <- file.path(input_dir,"006_EQTL_data")

#################
##  Functions  ##
#################

# Load the results of the Hill function fits
type_of_cluster_name <- c("Gradual_UP","Gradual_DOWN","Transient_UP","Transient_DOWN")
# names <- paste0("Hill_fits_wo_constrains_",type_of_cluster_name[1],"_responses.rds")
names <- paste0("Hill_fits_w_constrains_",type_of_cluster_name[1],"_responses.rds")

file_dir <- file.path(eqtl_input_dir,names)

results_list <- lapply(file_dir, readRDS)
names(results_list) <- type_of_cluster_name[1]

# Load Fit by cluster data frame
cluster_vec <- c(3,8,9,14)
cluster_name <- paste0("Cluster_", cluster_vec)


fit_by_cluster_list <- lapply(cluster_vec, function(i) {
  # fit_filename <- paste0("Hill_fits_Cluster_", i, ".txt")
  fit_filename <- paste0("Hill_fits_w_constrains_Cluster_", i, ".txt")
  read.table(file.path(eqtl_input_dir, fit_filename),header = TRUE, sep = "\t")
})

fit_by_cluster_list <- setNames(fit_by_cluster_list, cluster_name)

### Inspect the data
### These are the df with the fits parameters
df3_fit  <- fit_by_cluster_list[[1]]
df3_fit <- df3_fit %>% 
  drop_na()

df <- df3_fit

lower_Tao <- min(df$Tao, na.rm = TRUE)
upper_Tao <- max(df$Tao, na.rm = TRUE)
epsilon <- 0.01

df <- df %>%
  filter(R2 >= 0)
  # filter(Tao >= lower_Tao + epsilon &
  #        Tao <= upper_Tao - epsilon)

# Minimum number of individuals  to consider per gene
min_individuals_vec <- c(30,40, 50, 60, 70)

# Thresholds for R² values
r2_thresholds <- seq(0.5, 0.9, by = 0.02)

# Results data frame to store the number of genes for each combination of R² 
# threshold and minimum individuals
results <- expand.grid(
  R2_threshold = r2_thresholds,
  Min_Individuals = min_individuals_vec
)

# Function to calculate the number of genes for each combination of R² threshold
results$Num_Genes <- mapply(function(r2_thresh, min_indiv) {
  
  df %>%
    group_by(Gene) %>%
    summarise(n_indiv = sum(R2 >= r2_thresh, na.rm = TRUE), .groups = "drop") %>%
    filter(n_indiv >= min_indiv) %>%
    nrow()
}, results$R2_threshold, results$Min_Individuals)

results$Id <- paste0("R2_",results$R2_threshold, "_Ind_", results$Min_Individuals,"_num_genes_", results$Num_Genes)

list_min_ind <- list()
list_min_ind <- mapply(function(r2_thresh, min_indiv) {
  
  filtered_genes <- df %>%
    group_by(Gene) %>%
    summarise(n_indiv = sum(R2 >= r2_thresh, na.rm = TRUE), .groups = "drop") %>%
    filter(n_indiv >= min_indiv) %>%
    pull(Gene)
  
  # Filter the original data frame to include only the selected genes and individuals with R² >= threshold
  
  df %>%
    filter(Gene %in% filtered_genes, R2 >= r2_thresh)
  
}, results$R2_threshold, results$Min_Individuals, SIMPLIFY = FALSE)

list_min_ind <- setNames(list_min_ind, results$Id)


total_genes <- length(unique(df$Gene))
results$Porcentage_Genes <- round((results$Num_Genes / total_genes) * 100, 2)

# Plot the results
plt <- ggplot(results, aes(x = R2_threshold, y = Num_Genes, color = as.factor(Min_Individuals))) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Genes vs. R² threshold and minimum individuals",
    x = "Minimum R² threshold",
    y = "Number of genes",
    color = "Min Individuals"
  ) +
  scale_y_continuous(breaks = seq(0, max(results$Num_Genes), by = 20))+
  theme_minimal() +
  scale_color_brewer(palette = "Set1")+
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold"),
  )

plt

# Save the results to a file

write.table(results, file = file.path(input_dir,"006_EQTL_data",
                                      "genes_r2_threshold_and_min_individuals.txt"), 
            sep = "\t", quote = FALSE)
saveRDS(list_min_ind, file = file.path(input_dir,"006_EQTL_data",
                                       "genes_r2_threshold_and_min_individuals.rds"))


# Save the plot
r2_th_dir <- file.path(output_dir,"001_Figures", "012_R2_thresholds_vs_genes")
dir.create(r2_th_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file = file.path(r2_th_dir, "genes_vs_r2_threshold_and_min_individuals.pdf"), width = 8, height = 6)
print(plt)
dev.off()



























