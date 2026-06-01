#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(openxlsx)
  library(broom)
  library(patchwork)
  library(ggpubr)
}

### Set dir
main_dir <- getwd()
setwd(main_dir)
input_dir      <- "Analyses/Inputs"
script_id      <- "005_Test_filtering_options"
output_dir     <- file.path("Analyses/Outputs",script_id)
eqtl_input_dir <- file.path(input_dir,"006_EQTL_data")

dir.create(file.path(output_dir),recursive = T,showWarnings = F)

#################
##  Functions  ##
#################
# Define a publication-ready theme
theme_pub <- theme_bw(base_size = 16) +
  theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black", size = 12)
  )
# Create the function to test different strategies
#Inputs
# 1. The dataframe with the Hill function fits
# 2. The dataframe with the response to data
# 3. The dataframe with predicted values
# 4. Threshold for the R^2 value
# 5. Threshold for the median
# 6. Cluster id

test_filtering_options <- function(df_fit,df_response,df_predicted,
                                   th_R2,th_median,cluster_id){
  
  # Filter the data frame by th_R2 and th_median
  df_fit_filt <- df_fit %>% 
    filter(R2 > th_R2,
           Amp > 0)

  # Transform the data frame to long format to compute the median and select 
  # the genes
  df_wide <- df_fit_filt %>%
    select(Gene, Individuo, R2) %>%
    pivot_wider(names_from = Individuo, values_from = R2) %>% 
    mutate(R2_Median = apply(select(., -Gene), 1, median, na.rm = TRUE)) %>%
    filter(R2_Median > th_median)
  
  df_long <- df_wide %>%
    pivot_longer(cols = -c(Gene,R2_Median), names_to = "Individuo", values_to = "R2") %>%
    filter(!is.na(R2))
  
  # Select the genes that pass the filtering
  surviving_genes <- df_wide$Gene
  
  # Filter parameter data based on R2 thresholds
  param_filt <- df_fit_filt %>% filter(Gene %in% surviving_genes)
  
  ### check how many individuals by gene we have
  df_ind_by_gene <- param_filt %>% 
    group_by(Gene) %>% 
    summarise(Individuos = n_distinct(Individuo)) %>% 
    arrange(Individuos)
  
  ### check how many genes by individuals we have
  df_gen_by_ind <- param_filt %>% 
    group_by(Individuo) %>% 
    summarise(Genes = n_distinct(Gene)) %>% 
    arrange(Genes)
  
  gen_ind_path <- file.path(eqtl_input_dir,"001_Count_df_after_filtering")
  dir.create(gen_ind_path,recursive=T,showWarnings = FALSE)
  
  write.table(df_ind_by_gene,
              file.path(gen_ind_path,paste0("001_individuals_per_gene_",cluster_id,"_R2_",th_R2,"_median_",th_median,".txt")))
  write.xlsx(df_ind_by_gene,
             file.path(gen_ind_path,paste0("001_individuals_per_gene_",cluster_id,"_R2_",th_R2,"_median_",th_median,".xlsx")))
  write.table(df_gen_by_ind,
              file.path(gen_ind_path,paste0("002_genes_per_individual_",cluster_id,"_R2_",th_R2,"_median_",th_median,".txt")))
  write.xlsx(df_gen_by_ind,
             file.path(gen_ind_path,paste0("002_genes_per_individual_",cluster_id,"_R2_",th_R2,"_median_",th_median,".xlsx")))
  
  q25_ind <- quantile(as.numeric(df_ind_by_gene$Individuos), 0.25)
  q50_ind <- quantile(as.numeric(df_ind_by_gene$Individuos), 0.50)
  q75_ind <- quantile(as.numeric(df_ind_by_gene$Individuos), 0.75)
  
  ind_by_gene_plt <- ggplot(df_ind_by_gene,
                            aes(x = reorder(Gene, -Individuos), y = Individuos)) +
    
    geom_col(fill = "#4C72B0", width = 0.7, alpha = 0.9) +
    
    # Quantiles lines
    geom_hline(yintercept = q75_ind, color = "#D55E00", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = q50_ind, color = "#0072B2", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = q25_ind, color = "#009E73", linetype = "dashed", linewidth = 1) +
    
    # Add labels only once, nicely aligned
    annotate("text", x = 0.5, y = q75_ind, 
             label = sprintf("Q3 = %.1f", q75_ind),
             hjust = 0, vjust = -0.5, size = 4) +
    annotate("text", x = 0.5, y = q50_ind, 
             label = sprintf("Median = %.1f", q50_ind),
             hjust = 0, vjust = -0.5, size = 4) +
    annotate("text", x = 0.5, y = q25_ind, 
             label = sprintf("Q1 = %.1f", q25_ind),
             hjust = 0, vjust = -0.5, size = 4) +
    
    labs(title = "Number of Individuals per Gene",
         x = "",
         y = "Individuals") +
    
    theme_pub +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  q25_gene <- quantile(as.numeric(df_gen_by_ind$Genes), 0.25)
  q50_gene <- quantile(as.numeric(df_gen_by_ind$Genes), 0.50)
  q75_gene <- quantile(as.numeric(df_gen_by_ind$Genes), 0.75)
  
  gen_by_ind_plt <- ggplot(df_gen_by_ind,
                           aes(x = reorder(Individuo, -Genes), y = Genes)) +
    
    geom_col(fill = "#4C72B0", width = 0.7, alpha = 0.9) +
    
    geom_hline(yintercept = q75_gene, color = "#D55E00", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = q50_gene, color = "#0072B2", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = q25_gene, color = "#009E73", linetype = "dashed", linewidth = 1) +
    
    annotate("text", x = 0.5, y = q75_gene, label = sprintf("Q3 = %.1f", q75_gene),
             hjust = 0, vjust = -0.5, size = 4) +
    annotate("text", x = 0.5, y = q50_gene, label = sprintf("Median = %.1f", q50_gene),
             hjust = 0, vjust = -0.5, size = 4) +
    annotate("text", x = 0.5, y = q25_gene, label = sprintf("Q1 = %.1f", q25_gene),
             hjust = 0, vjust = -0.5, size = 4) +
    
    labs(title = "Number of Genes per Individual",
         x = "Individual",
         y = "Number of Genes") +
    
    theme_pub +
    theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 10))
  
  # Save the plots
  count_path <- file.path(output_dir,"001_Figures","008_Count_plots_R2")
  dir.create(count_path,recursive = T, showWarnings = FALSE)
  
  pdf(file = file.path(count_path,paste0(cluster_id,"_R2_",th_R2,"_median_",th_median,".pdf")), width = 10,height = 8)
  print(ind_by_gene_plt)
  print(gen_by_ind_plt)
  dev.off()
  
  # Get the matrices for the Hill function parameters
  df_wide_amp <- param_filt %>%
    select(Gene, Individuo, Amp) %>%
    pivot_wider(
      names_from = Individuo,
      values_from = Amp
    )
  df_wide_tao <- param_filt %>%
    select(Gene, Individuo, Tao) %>%
    pivot_wider(
      names_from = Individuo,
      values_from = Tao
    )
  
  df_wide_h <- param_filt %>%
    select(Gene, Individuo, h) %>%
    pivot_wider(
      names_from = Individuo,
      values_from = h
    )
  
  # Return summary
  summary_list <- list(
    th_R2      = th_R2,
    th_median  = th_median,
    n_genes    = length(unique(param_filt$Gene)),
    n_obs      = nrow(param_filt),
    df_stat_long = param_filt,
    df_amp     = df_wide_amp,
    df_h       = df_wide_h,
    df_tao     = df_wide_tao
  )
  return(list(
    summary = summary_list,
    df_ind_by_gene = df_ind_by_gene,
    df_gen_by_ind = df_gen_by_ind,
    ind_by_gene_plt = ind_by_gene_plt,
    gen_by_ind_plt = gen_by_ind_plt
  ))
  
  summary_path <- file.path(eqtl_input_dir,"summary_filtering_test")
  dir.create(summary_path, showWarnings = FALSE)
  saveRDS(summary_list,
          file.path(summary_path,paste0("003_summary_filtering_test_",cluster_id,"_R2_",th_R2,"_median_",th_median,".rds")))
}
#################
##  Load data  ##
#################

# Load the files named list from previous analysis
list_of_files <- list.files(eqtl_input_dir, pattern = "Hill_fits_wo_constrains_Cluster_.*\\.txt$")

# Load a list of data frames by cluster
hill_fit_list <- lapply(list_of_files, function(file) {
  read.table(file.path(eqtl_input_dir, file), header = TRUE, sep = "\t")
})
names(hill_fit_list) <- gsub("Hill_fits_wo_constrains_(Cluster_.*)\\.txt", "\\1", list_of_files)

# Create a table with the th_R2 and th_median values
R2     <- c(0.6, 0.7, 0.8, 0.9)
median <- c(0.7, 0.8, 0.9, 0.95)

# R2     <- c(0.9)
# median <- c(0.9)
# Create a data frame with the thresholds

threshold <- data.frame(
  th_R2 = rep(R2,each = length(median)),
  th_median = median
)

# Filter the data frame to keep only the combinations where th_median is 
# greater than or equal to th_R2
threshold <- threshold %>% 
  filter(th_median >= th_R2) %>%
  arrange(th_R2, th_median)
threshold

# Run the function for each cluster and each threshold
final_list <- list()
cluster_vec <- names(hill_fit_list)
vec <- seq_along(cluster_vec)[2]

type_of_cluster_name <- c("Gradual_UP","Gradual_DOWN","Transient_UP","Transient_DOWN")

# Do a for loop to process all clusters
for (i in vec) {
  
  # Select clusters to analyze
  df <- hill_fit_list[[cluster_vec[i]]] # Select cluster #3 who is at position 2
  # Let´s add a setup column to the data frame
  df$Setup  <- paste0("Fit_",df$Gene,"_",df$Individuo)
  
  type_of_cluster <- type_of_cluster_name[1]
  cluster_number  <- cluster_vec[i]
  
  print(paste0("Processing ", type_of_cluster, " ", cluster_number))
  
  cluster_id <- paste0(type_of_cluster,"_",cluster_number)
  
  # df_fit = stats_by_cluster_list[[i]]
  df_fit = df
  
  nan_index  <- which(apply(df_fit, 1, function(x) any(is.na(x))))
  df_fit_nan  <- df_fit[nan_index,c(1:3)]
  
  print(paste0("We have ",length(unique(df_fit_nan$Gene))," genes with at least one NaN values in the fit parameters."))
  print(paste0("We have ",length(unique(df_fit_nan$Individuo))," individuals with at least one NaN values in the fit parameters."))
  
  df_response  = df_fit[,grep("^exp_",colnames(df_fit))]
  df_predicted = df_fit[,grep("^fit_",colnames(df_fit))]
    
  for (j in seq_len(nrow(threshold))) {
    
    th_R2      <- threshold[j,1]
    th_median  <- threshold[j,2]
    
    result <- test_filtering_options(
      df_fit       = df_fit,
      df_response  = df_response,
      df_predicted = df_predicted,
      th_R2        = th_R2,
      th_median    = th_median,
      cluster_id   = cluster_id
    )
    
    final_list[[cluster_id]][[paste0("R2_",th_R2,"_median_",th_median)]] <- result
    
    result_path <- file.path(eqtl_input_dir,"Summary_list")
    dir.create(result_path, showWarnings = FALSE)
    # Save the results
    saveRDS(result,
            file.path(result_path,paste0("Summary_filtering_test_",cluster_id,"_R2_",th_R2,"_median_",th_median,".rds")))
  }
}

# Save the final list
saveRDS(final_list,
        file.path(eqtl_input_dir,"Summary_list","Final_summary_filtering_test.rds"))
