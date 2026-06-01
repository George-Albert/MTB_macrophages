
#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(openxlsx)
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
##  Load Data  ##
#################

# Load the delta expression data
delta_expr_df <- read.table(file.path(eqtl_input_dir, "delta_expr_df.txt"),header = TRUE, sep = "\t")

# Load the gene parameters for Tau between 0 and 72 hours
# clus2_df <- read.table(file.path(eqtl_input_dir, "Cluster_2_h_stats.txt"), sep = "\t")
stat_df  <- read.table(file.path(eqtl_input_dir, "Tao_filt_stats.txt"), sep = "\t")

# Plot h for each h_cluster
my_colors <- c("#DF535B",
               "#61D03F",
               "#2296E6")

plt_h <- ggplot(stat_df, aes(x = as.factor(h_cluster), y = h, colour=as.factor(h_cluster))) +
  geom_point() +
  labs(title = "Distribution of h values across clusters",
       x = "h Cluster",
       y = "h Value") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_color_manual(values = my_colors) +
  theme_minimal()
plt_h

# Save the plot
cor_plt_dir <- file.path(output_dir,"001_Figures", "011_Correlation_parameters")
dir.create(cor_plt_dir, recursive = TRUE, showWarnings = FALSE)

pdf(file.path(cor_plt_dir, "h_values_distribution_across_clusters.pdf"), width = 8, height = 6)
print(plt_h)
dev.off()

#  Plot h for one cluster at a time

for (i in unique(stat_df$h_cluster)) {
  
  plt_h_i <- ggplot(stat_df[stat_df$h_cluster == i, ], aes(x = h)) +
    geom_histogram(bins = 30, fill = my_colors[i],color = "black", alpha = 0.7) +
    labs(title = paste("Distribution of h values for cluster", i),
         x = "h Value",
         y = "Frequency") +
    theme_minimal()
  
  pdf(file.path(cor_plt_dir, paste0("h_values_distribution_cluster_", i, ".pdf")), width = 8, height = 6)
  print(plt_h_i)
  dev.off()
}

# Plot Hill values for each h_cluster

hill_function <- function(time, Amp, Tao, h) {
  Exp=Amp * time^h / (Tao^h + time^h)
  return(Exp)
}

# Create a sequence of time points from 0 to 72 hours
time_points <- seq(0, 72, by = 1)
# # Create an empty list to store the plots
# list_hill_plots <- list()
# Loop through each h_cluster and create the Hill plot
clus_vec <- unique(stat_df$h_cluster)
for (i in clus_vec) {
  
  message(paste("Processing h_cluster:", i))
  # Create a data frame to store the Hill values for each gene and individual
  hill_df_by_cluster <- data.frame()
  
  # Filter the data for the current h_cluster
  cluster_data <- stat_df[stat_df$h_cluster == i, ]
  all_genes <- unique(cluster_data$Gene)
  # start_pos <- match("TTYH2", all_genes)
  # all_genes[start_pos:length(all_genes)]
 
  for(gene_id in all_genes){
    
    message(paste("Processing gene:", gene_id))
    # Filter the data for the current gene
    gene_data <- cluster_data %>% 
      filter(Gene %in% gene_id) 
      
    all_ind <- unique(gene_data$Individuo)
    
    for (ind_id in all_ind) {
      
      message(paste("Processing individual:", ind_id))
      
      ind_data <- gene_data %>% 
        filter(Individuo %in% ind_id)
      
      if (nrow(ind_data) == 0) {
        next
      }
      
      # Calculate the Hill values for the current gene and individual
      hill_df <- data.frame(
        gene = gene_id,
        ind = ind_id,
        time = time_points,
        h = ind_data$h[1],      # Assuming h is constant for the individual
        Amp = ind_data$Amp[1],  # Assuming Amp is constant for the individual
        Tao = ind_data$Tao[1]   # Assuming Tao is constant for the individual
      )
      
      hill_df$Hill <- hill_function(time_points, ind_data$Amp, ind_data$Tao, ind_data$h)

      hill_df_by_cluster <- rbind(hill_df_by_cluster, hill_df)
      
    }
    #Save the Hill values to a file
    write.table(hill_df_by_cluster, 
                file = file.path(eqtl_input_dir, paste0("Hill_values_cluster_", i,".txt")), 
                sep = "\t", quote = FALSE)
  }
  
  # Create the Hill plot for the current gene and individual
  hill_plot <- ggplot(hill_df_by_cluster, aes(x = time, y = Hill,group = interaction(gene, ind), color = as.factor(ind))) +
    geom_line(linewidth = 1,alpha=0.5) +
    labs(title = paste("Hill Function for Cluster", i),
         x = "Time (hours)",
         y = "Hill Value",
         color = "Individual") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      # panel.background = element_rect(fill = "black", color = NA),
      panel.grid.major = element_line(color = "gray30"),
      panel.grid.minor = element_line(color = "gray20"),
      # legend.background = element_rect(fill = "black", color = NA),
      # legend.key = element_rect(fill = "black", color = NA),
      text = element_text(color = "white"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      plot.title = element_text(color = "white")
    )
  
  hill_plot
  # Save the plot as a PDF
  pdf(file.path(cor_plt_dir, paste0("Hill_function_cluster_", i,"_black.pdf")), width = 20, height = 10)
  print(hill_plot)
  dev.off()
  
  message(paste("Hill plot for cluster", i, "saved."))
  
}
  

####################
##  Plot h vs R2  ##
####################

# Plot h vs R2

plt_h_r2 <- ggplot(stat_df, aes(x = h, y = R2, color = as.factor(h_cluster))) +
  geom_point() +
  labs(title = "R2 vs h values",
       x = "h Value",
       y = "R2 Value") +
  scale_color_manual(values = my_colors) +
  theme_minimal()
plt_h_r2

# Save the plot

pdf(file.path(cor_plt_dir, "h_vs_R2.pdf"), width = 10, height = 8)
print(plt_h_r2)
dev.off()


############################
##  Correlation Tau vs h  ##
############################

# Correlation between Tau and h values colored by R2

plt_tau_h <- ggplot(stat_df, aes(x = Tao, y = h, color = R2)) +
  geom_point() +
  labs(title = "Correlation between Tau and h values",
       x = "Tau Value",
       y = "h Value") +
  scale_color_gradient(low = "yellow", high = "red3") +
  theme_minimal()

plt_tau_h

# Save the plot

pdf(file.path(cor_plt_dir, "Tau_vs_h_colored_by_R2.pdf"), width = 10, height = 8)
print(plt_tau_h)
dev.off()

# Now Select 6 groups and plot the correlation between Tau and h values 
# colored by R2
# 
# plot( stat_df$Tao, stat_df$h,
#      pch=19,
#      xlab="h values", 
#      ylab="Tau values", 
#      main="Correlation between Tau and h values colored by R2")
# 
# 
# rownames(stat_df) <- c(1:nrow(stat_df))
# identify(stat_df$Tao, stat_df$h, 
#          labels = rownames(stat_df), 
#          pch = 19)

selected_data_storage <- read.table(file.path(eqtl_input_dir, "grupos_seleccionados_Tau_vs_h.txt"), header = TRUE, sep = "\t")  
selected_data_storage <- selected_data_storage %>%
  distinct(Gene, Individuo, .keep_all = TRUE)

# Extrae los grupos que quieres comparar
# group_6 <- selected_data_storage %>% 
#   filter(Grupo == "Grupo_6") %>% 
#   distinct() 
# group_1 <- selected_data_storage %>% 
#   filter(Grupo == "Grupo_1") %>% 
#   distinct()
# group_2 <- selected_data_storage %>%
#   filter(Grupo == "Grupo_2") %>% 
#   distinct()
# group_3 <- selected_data_storage %>%
#   filter(Grupo == "Grupo_3") %>% 
#   distinct()
# group_4 <- selected_data_storage %>%
#   filter(Grupo == "Grupo_4") %>% 
#   distinct()
# group_5 <- selected_data_storage %>%
#   filter(Grupo == "Grupo_5") %>% 
#   distinct()
#   
group_345 <- selected_data_storage %>% filter(Grupo %in% c("Grupo_3", "Grupo_4", "Grupo_5"))
# group_12345 <- bind_rows(group_1, group_2, group_3, group_4, group_5) %>%
#   distinct()

# Subset de grupo 6 excluyendo coincidencias exactas con grupos 3, 4 y 5
group_6_cleaned <- anti_join(group_6, group_345, by = colnames(selected_data_storage)[-ncol(selected_data_storage)])

selected_data_storage <- selected_data_storage %>% 
  filter(Grupo != "Grupo_6")


combined_groups <- bind_rows(group_6_cleaned, selected_data_storage)

write.table(combined_groups, 
            file = file.path(eqtl_input_dir, "grupos_seleccionados_Tau_vs_h_cleaned.txt"), 
            sep = "\t", quote = FALSE)

# Plot the selected groups colored by R2

plot_group_correlations <- function(data, group_id, save = FALSE, out_dir = ".") {
  
  
  group_data <- data %>% filter(Grupo == group_id)
  
  # Calcular correlaciones
  cor_h_tao <- cor(group_data$h, group_data$Tao, use = "complete.obs", method = "pearson")
  cor_h_r2  <- cor(group_data$h, group_data$R2,  use = "complete.obs", method = "pearson")
  
  # Plot 1: Tao vs h
  plt1 <- ggplot(group_data, aes(x = Tao, y = h, color = R2)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = paste0(group_id, ": Tao vs h  (Pearson r = ", round(cor_h_tao, 3), ")"),
      x = "Tau", y = "h Value", color = "R²"
    ) +
    scale_color_gradient(low = "yellow", high = "red3") +
    theme_minimal()
  
  # Plot 2: h vs R2
  plt2 <- ggplot(group_data, aes(x = h, y = R2, color = R2)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = paste0(group_id, ": h vs R²  (Pearson r = ", round(cor_h_r2, 3), ")"),
      x = "h Value", y = "R²", color = "R²"
    ) +
    scale_color_gradient(low = "yellow", high = "red3") +
    theme_minimal()
  
  # Mostrar plots
  print(plt1)
  print(plt2)
  
  # Guardar si se solicita
  if (save) {
    
    pdf(file.path(out_dir, paste0("Correlation_plots_", group_id, "_h_vs_Tau.pdf")), width = 12, height = 6)
    print(plt1)
    dev.off()
    pdf(file.path(out_dir, paste0("Correlation_plots_", group_id, "_h_vs_R2.pdf")), width = 12, height = 6)
    print(plt2)
    dev.off()

  }
}

# Loop through each group and plot the correlations
unique_groups <- factor(unique(combined_groups$Grupo),
                        levels = c("Grupo_1","Grupo_2","Grupo_3",
                                   "Grupo_4", "Grupo_5", "Grupo_6"))
for (group in unique_groups) {
  
  message(paste("Processing group:", group))
  
  plot_group_correlations(combined_groups, group, save = TRUE, out_dir = cor_plt_dir)
  
}


plt1 <- ggplot(combined_groups[combined_groups$Grupo=="Grupo_1",], aes(x = Amp, y = Tao, color = R2)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    # title = paste0(group_id, ": Tao vs h  (Pearson r = ", round(cor_h_tao, 3), ")"),
    x = "Amp", y = "h Value", color = "R²"
  ) +
  scale_color_gradient(low = "yellow", high = "red3") +
  theme_minimal()
plt1
