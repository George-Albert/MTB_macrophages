
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
script_id      <- "004_Inspect_R2_negative" 
input_dir      <- "Analyses/Inputs"
output_dir     <- file.path("Analyses/Outputs",script_id)
eqtl_input_dir <- file.path(input_dir,"006_EQTL_data")

dir.create(file.path(output_dir),recursive = T,showWarnings = F)

#################
##  Functions  ##
#################


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

# Select clusters to analyze
df <- hill_fit_list[["Cluster_3"]]
# Let´s add a setup column to the data frame
df$Setup  <- paste0("Fit_",df$Gene,"_",df$Individuo)

# Filter all the obs that have R2 less than 0
df_fit_r2_neg <- df %>%
  filter(R2 < 0) 

# Filter the predicted vallues dataframe by df_fit  
match_ind <- df_fit_r2_neg$Setup

df_long <- df_fit_r2_neg %>%
  pivot_longer(
    cols = starts_with("exp_") | starts_with("fit_"),
    names_to = c(".value", "Time"),
    names_pattern = "(exp|fit)_(\\d+)"
  ) %>%
  mutate(
    Time = as.numeric(Time)
  )
##########################################################################
##  Plot the data in a for loop and use wrap to create multiples plots  ##
##########################################################################

# Create output directory for plots
cluster_vec <- unique(df$Cluster)
hill_out_dir <- file.path(output_dir,"001_Figures",paste0("008_Hill_fits_wo_constrains_fits_r2_neg_Cluster",cluster_vec))
dir.create(hill_out_dir, recursive = TRUE, showWarnings = FALSE)

# Loop through each individual and create plots
id <- unique(df_long$Setup)
id_gene <- unique(df_long$Gene)

# Desired plot dimensions in inches
plot_width  <- 3.5
plot_height <- 3.5

# time grid for full predicted curve
time_seq <- seq(0, 72, by = 1)

for (g in id_gene) {
  
  df_to_plot <- df_long %>% filter(Gene == g)
  
  # Create a data frame for the fitted curves
  df_params <- df_to_plot %>%
    distinct(Setup, Gene, Individuo, Amp, h, Tao)
  
  # Generate fitted curves using the Hill equation
  df_curve <- df_params %>%
    mutate(curve = pmap(
      list(Amp, h, Tao, Setup, Gene, Individuo),
      function(Amp, h, Tao, Setup, Gene, Individuo) {
        
        fit_vals <- (Amp * time_seq^h) / (time_seq^h + Tao^h)
        
        # Detect NaNs in the Hill curve
        if (any(is.nan(fit_vals))) {
          message(" NaN detected for: ",
                  " Gene=", Gene,
                  " | Individuo=", Individuo,
                  " | Setup=", Setup)
        }
        
        tibble(Time = time_seq, fit_full = fit_vals)
      }
    )) %>%
    select(Setup, Gene, Individuo, curve) %>%
    unnest(curve)
  
  
  setups_gene <- unique(df_to_plot$Setup)
  n_panels    <- length(setups_gene)
  
  # Dinamically determine nrow and ncol based on number of panels
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels/ncol)
  
  total_width  <- ncol * plot_width
  total_height <- nrow * plot_height
  
  message(sprintf("Plotting gene %s with layout %dx%d", g, ncol, nrow))
  
  df_zero <- data.frame(
    Time = 0,
    Expression = 0,  # o 0, o cualquier valor de referencia
    Setup = unique(df_to_plot$Setup) # si usas facet_wrap por Setup
  )
  
  p <- ggplot() +
    # Model with blue line
    geom_line(data = df_curve, aes(x = Time, y = fit_full),
              color = "#1f78b4", linewidth = 1.2) +
    # observed data
    geom_point(data = df_to_plot, aes(x = Time, y = exp),
               color = "black", size = 2) +
    # Add zero orange point
    geom_point(data = df_zero, aes(x = Time, y = Expression), 
               color = "#e66101", size = 3) +
    facet_wrap(~Setup, ncol = ncol, scales = "fixed") +
    # labels
    labs(title = paste("Gene:", g),
         x = "Time (h)", 
         y = "Expression level") +
    # axes with defined ticks
    scale_x_continuous(breaks = c(0, 2, 20, 48, 72)) +
    # Clean theme for publications
    theme_bw(base_size = 14) +  # font size
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text  = element_text(size = 12, color = "black"),
      strip.text = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, linetype = "dashed", color = "grey80"),
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(hill_out_dir, sprintf("%s.pdf", g)),
    plot   = p,
    width  = total_width,
    height = total_height,
    units  = "in"
  )
}
























