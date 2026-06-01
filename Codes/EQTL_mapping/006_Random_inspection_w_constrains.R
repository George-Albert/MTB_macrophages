#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(openxlsx)
  library(ggpubr)
}

### Set dir
main_dir <- getwd()
setwd(main_dir)
input_dir      <- "Analyses/Inputs"
output_dir     <- "Analyses/Outputs"
eqtl_input_dir <- file.path(input_dir,"006_EQTL_data")

#################
##  Functions  ##
#################

random_inspection_plot <- function(summary_list, stats_df) {
  
  for (summary_name in names(summary_list)) {
    
    # Extract the cluster ID from the summary name
    cluster_id <- gsub(".*Cluster_", "Cluster_", summary_name)
    cluster_id
    
    total_genes <- length(unique(stats_df$Gene))
    gene_idx <- unique(stats_df$Gene)
    if (total_genes < 3) next
    
    max_genes <- 25
    n_to_select <- min(max_genes, total_genes)
    indices <- round(seq(1, total_genes, length.out = n_to_select))
    selected_genes <- gene_idx[indices]
    
    for (i in seq_along(parameters)) {
      
      param     <- parameters[i]
      # Filter the fit_df
      filtered_df <- stats_df %>%
        filter(Gene %in% selected_genes, !is.na(.data[[param]]))
      
      if (nrow(filtered_df) == 0) next
      
      # Plot
      p_before <- ggviolin(filtered_df, x = "Gene", y = param, fill = "Gene") +
        geom_jitter(width = 0.1, alpha = 0.5) +
        # ggtitle(paste0(param, " before"))
        # scale_y_log10("Log_scale") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10)) +
        labs(title = paste(summary_name,"-",param, " before"))
      
      # Save the plot
      filename <- file.path(random_pvalues_path_before, paste0(summary_name, "_", param, ".pdf"))
      ggsave(filename, plot = p_before, width = 16, height = 8)
      
      
      #########################
      ##  Filtering process  ##
      #########################
      
      filtered_df <- stats_df %>%
        filter(Gene %in% selected_genes, !is.na(.data[[param]]))
      
      
      # --- Filtering step: remove extreme 1% per gene
      df_filtered <- filtered_df %>%
        group_by(Gene) %>%
        mutate(
          lower = quantile(.data[[param]], 0.01, na.rm = TRUE),
          upper = quantile(.data[[param]], 0.99, na.rm = TRUE)
        ) %>%
        filter(.data[[param]] > lower, .data[[param]] < upper) %>%
        ungroup()
      
      if (nrow(df_filtered) == 0) next
      
      # Plot after filtering
      p_after <- ggviolin(df_filtered, x = "Gene", y = param, fill = "Gene") +
        geom_jitter(width = 0.1, alpha = 0.5) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10)) +
        labs(title = paste(summary_name, "-", param, "- after filtering"))
      
      ggsave(file.path(random_pvalues_path_after, paste0(summary_name, "_", param, "_after.pdf")),
             plot = p_after, width = 16, height = 8)
      
      
      
      }
    }
  }
 
# Create a function to plot the distribution of parameters
plot_distribution <- function(stats_df, param, output_path, before_after = "before",width = 8, height = 6,color=c("red", "blue")) {         
  
  gene_colors <- stats_df %>%
    distinct(Gene) %>%
    arrange(Gene) %>%
    mutate(ColorGroup = rep(c(color[1], color[2]), length.out = n()))
  
  stats_df1 <- stats_df %>%
    left_join(gene_colors, by = "Gene")
  
  param <- rlang::sym(param) # Convert the parameter name to a symbol for ggplot
  p <- ggplot(stats_df1, aes(x = Gene, y = !!param, color = ColorGroup)) +
    geom_point(size = 3, alpha = 0.8) +            # puntos más grandes y semitransparentes
    labs(title = paste(param, "per Gene"),
         x = "Gene", y = param) +
    scale_color_manual(values = unique(stats_df1$ColorGroup)) +
    theme_minimal(base_size = 14) +                # fuente más grande para poster
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10), 
      axis.title = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14)
    )
  
  pdf(file = file.path(output_path, paste0("Distribution_", param,"_",before_after, ".pdf")), width = width, height = height)
  print(p)
  dev.off()
  return(p)
}

#################
##  Load data  ##
#################

# Load the list of all individuals 
final_summary_list_path <- file.path(eqtl_input_dir,"Summary_list")
final_summary_list <- readRDS(file.path(final_summary_list_path, "Final_summary_filtering_test.rds"))

type_of_cluster_name <- c("Gradual_UP","Gradual_DOWN","Transient_UP","Transient_DOWN")
# Load the dataframes with the Hill function fits using lapply
cluster_vec <- c(3,8,9,14)
cluster_number <- paste0("Cluster_", cluster_vec)

stats_by_cluster_list <- lapply(cluster_vec, function(i) {
  read.table(file.path(eqtl_input_dir, paste0("Hill_fits_w_constrains_Cluster_", i, ".txt")),header = TRUE, sep = "\t")
})

stats_by_cluster_list <- setNames(stats_by_cluster_list, cluster_number)
stats_by_cluster_list <- read.table(file.path(eqtl_input_dir, paste0("Hill_fits_wo_constrains_Cluster_3.txt")),header = TRUE, sep = "\t")

# Extract the dataframes from the final_summary_list
# names_of_list      <- names(final_summary_list)
# r2_combination_vec <- names(final_summary_list[[names_of_list[1]]])
# 
# r2_selected_list      <- final_summary_list[[names_of_list[1]]][[r2_combination_vec[10]]]
# summary_list_selected <- r2_selected_list$summary
# df_shapiro            <- summary_list_selected$shapiro_results
# df_fit                <- stats_by_cluster_list[[1]]

# Create output directory
random_pvalues_path_before <- file.path(output_dir, "001_Figures", "009_Random_inspection_pvalues","Before")
dir.create(random_pvalues_path_before, showWarnings = FALSE,recursive = T)
random_pvalues_path_after  <- file.path(output_dir, "001_Figures", "009_Random_inspection_pvalues","After")
dir.create(random_pvalues_path_after, showWarnings = FALSE,recursive = T)

# Parameters to explore
parameters   <- c("Amp", "h", "Tao")
pvalue_names <- c("pval_Amp", "pval_h", "pval_Tao")
# Input parameters
summary_list <- final_summary_list
stats_df   <- final_summary_list$Gradual_UP_Cluster_3$R2_0.9_median_0.9$summary$df_stat_long

summary(stats_df)

### Plot the distribution of the parameters
# Create the output directory for the plots
distribution_plots_path <- file.path(output_dir, "001_Figures", "009_Random_inspection_pvalues", "Distribution_plots")
before_path <- file.path(distribution_plots_path, "Before_filtering")
dir.create(before_path, showWarnings = FALSE, recursive = TRUE)


# Generate distribution plots for each parameter
# Loop through each parameter and create the distribution plot
for (param in parameters) {
  plot_distribution(stats_df, param, before_path,before_after = "before",width = 100, height = 8)
}

lower_Tao <- 0
upper_Tao <- 80

Tao_filt <- stats_df %>%
  filter(Tao > lower_Tao, Tao < upper_Tao)%>% 
  mutate(Setup=paste0("Fit_",Gene,"_",Individuo))
  
q_amp <- quantile(Tao_filt$Amp, 0.99, na.rm = TRUE)
q_h <- quantile(Tao_filt$h, 0.99, na.rm = TRUE)

Tao_filt <- Tao_filt %>%
  filter(Amp <= q_amp, h <= q_h, h>0) %>% 
  mutate(Setup=paste0("Fit_",Gene,"_",Individuo))

summary(Tao_filt)
after_path <- file.path(distribution_plots_path, "After_filtering")
dir.create(after_path, showWarnings = FALSE, recursive = TRUE)

for (param in parameters) {
  plot_distribution(Tao_filt, param, after_path,before_after = "after_filter_Tau_amp_h_w_constrains",width = 100, height = 8)
}

#Plot histogram of the parameters after filtering
params <- c("Amp", "Tao", "h")
colors <- c("blue", "green", "red")

for (i in seq_along(params)) {
  
  param <- params[i]
  color <- colors[i]
  
  p <- ggplot(Tao_filt, aes(x = .data[[param]])) +
    geom_histogram(bins = 30, fill = color, alpha = 0.5) +
    labs(title = paste("Distribution of", param, "after filtering"), x = param, y = "Frequency") +
    theme_minimal()
  
  # Save the plot
  pdf(file = file.path(after_path, paste0("Histogram_", param, "_after_filtering_w_constrains.pdf")), width = 12, height = 8)
  print(p)
  dev.off()
}

# Plot correlation between Tao and h
cor_plt_tao_h <- ggplot(Tao_filt, aes(x=Tao,y= h))+
         geom_point() +
         geom_smooth(method = "lm", se = FALSE, color = "red") +
         labs(title = "Correlation between Tao and h after filtering",
              x = "Tao", y = "h") +
         theme_minimal()

pdf(file = file.path(after_path, "Correlation_Tao_h_after_filtering_w_constrains.pdf"), width = 12, height = 8)
print(cor_plt_tao_h)
dev.off()

### Cluster the data by h
set.seed(123) # Set seed for reproducibility
h_clusters <- kmeans(Tao_filt$h, centers = 3, nstart = 20)
# Add cluster information to the data frame
Tao_filt$h_cluster <- as.factor(h_clusters$cluster)

# Plot the distribution of h clusters

hist_h <- ggplot(Tao_filt, aes(x = h, fill = h_cluster)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.5) +
  labs(title = "Distribution of h clusters after filtering",
       x = "h", y = "Frequency") +
  theme_minimal()

pdf(file = file.path(distribution_plots_path, "h_cluster_hist_after_filtering_w_constrains.pdf"), width = 12, height = 8)
print(hist_h)
dev.off()


# Plot geom point of h clusters

cor_plt_tao_h <- ggplot(Tao_filt, aes(x = Tao, y = h, color = h_cluster)) +
  geom_point() +
  labs(title = "Correlation_Tao_h_after_filtering",
       x = "Tao", y = "h") +
  theme_minimal()

h_cluster_plot <- ggplot(Tao_filt, aes(x = Gene, y = h, color = h_cluster)) +
  geom_point() +
  labs(title = "h_clusters_distribution_after_filtering",
       x = "Gene", y = "h") +
  theme_minimal()
amp_plt <- ggplot(Tao_filt, aes(x = Gene, y = Amp, color = h_cluster)) +
  geom_point(alpha=0.5) +
  labs(title = "Amp per Gene after filtering",
       x = "Gene", y = "Amp") +
  theme_minimal()

tao_plt <- ggplot(Tao_filt, aes(x = Gene, y = Tao, color = h_cluster)) +
  geom_point() +
  labs(title = "Tao per Gene after filtering",
       x = "Gene", y = "Tao") +
  theme_minimal()
# Save the plots

pdf(file = file.path(distribution_plots_path, "h_histogram_after_filtering_w_constrains.pdf"), width = 12, height = 8)
print(hist_h)
dev.off()
pdf(file = file.path(distribution_plots_path, "Correlation_Tao_h_after_filtering_w_constrains.pdf"), width = 12, height = 8)
print(cor_plt_tao_h)
dev.off()
pdf(file = file.path(distribution_plots_path, "h_clusters_distribution_after_filtering_w_constrains.pdf"), width = 12, height = 8)
print(h_cluster_plot)
dev.off()
pdf(file = file.path(distribution_plots_path, "Amp_per_Gene_after_filtering_clustering_w_constrains.pdf"), width = 12, height = 8)
print(amp_plt)
dev.off()

### Lets inspect the clusters

length(which(Tao_filt$h_cluster == 1)) # 1 cluster)
# 3080
length(which(Tao_filt$h_cluster == 2)) # 2 cluster
# 11962
length(which(Tao_filt$h_cluster == 3)) # 3 cluster
# 1651

clus1_df <- Tao_filt[Tao_filt$h_cluster == 2,]
clus_df <- clus1_df

### Plot the distribution of the parameters for the cluster 1 colored by individuals

clus1_plot <- ggplot(clus_df, aes(x = Gene, y = Amp, color = Individuo)) +
  geom_point() +
  labs(title = "Cluster 1 - Amp per Gene colored by Individuals",
       x = "Gene", y = "Amp") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


write.table(clus_df, file = file.path(eqtl_input_dir, "Cluster_2_h_stats.txt"), sep = "\t")
write.table(Tao_filt, file = file.path(eqtl_input_dir, "Tao_filt_stats.txt"), sep = "\t")


# Call the function to generate plots
random_inspection_plot(final_summary_list, stats_df)
 
