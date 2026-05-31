#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(minpack.lm)
  library(broom)
}

### Set dir
main_dir <- getwd()
setwd(main_dir)
script_id      <- "002_Fitting_Hill_eq"
input_dir      <- "Analyses/Inputs"
output_dir     <- file.path("Analyses/Outputs",script_id)
eqtl_input_dir <- file.path(input_dir,"006_EQTL_data")

dir.create(output_dir,recursive = T,showWarnings = F)

#################
##  Load data  ##
#################

# Load the list of all individuals responses data
# all_individuals_data <- readRDS(file.path(eqtl_input_dir, "all_individuals_exp_data.rds"))
# Load the delta expression data
delta_expr_df <- read.table(file.path(eqtl_input_dir, "delta_expr_df.txt"),header = TRUE, sep = "\t")

# Create the vectors for the differents type of clusters
up_vec   <- c(3,5,7,13)            #  Gradual UP 
down_vec <- c(1,6,8,10)            #  Gradual DOWN

transient_up_vec    <- c(15,16,4)   #  Transient UP
transient_dowwn_vec <- c(12,9,2)    #  Transient DOWN

type_of_cluster_name <- c("Gradual_UP","Gradual_DOWN","Transient_UP","Transient_DOWN")

# Select the clusters of interest or do a loop over all clusters
cluster_list <- list(up_vec,down_vec,transient_up_vec,transient_dowwn_vec)  
# cluster_list <- cluster_list[[1]][[1]] # Select only the first cluster list for testing
cluster_list <- list(up_vec)

system.time({
  for (j in seq_along(cluster_list)) {
    
    cluster_vec <- cluster_list[[j]]
    message(paste0("Fitting Hill equation for ", type_of_cluster_name[j], " clusters"))
    # Define the results list with the length of the cluster vector
    results <- vector("list", length = length(cluster_vec))
    
    # Define the loop over the clusters vector
    vec <- seq_along(cluster_vec)
    
    # Define the loop over the clusters
    
    for (i in vec) {
      
      clus <- cluster_vec[i]
      print(clus)
      message(paste0("Fitting Hill equation for Cluster ", clus))
      # Define the name of the cluster in the results list
      cluster_name <- paste0("Cluster_", clus)
      names(results)[[i]] <- cluster_name
      # Define the name of the data frame in the results list
      df_name <- paste0("Expression_df")
      # Add the data frame to the results list
      df_cluster <- delta_expr_df[delta_expr_df$Cluster == clus, ]
      results[[cluster_name]][[1]] <- df_cluster
      names(results[[cluster_name]])[1] <- df_name
      
      # Obtain the genes for the cluster
      gene_vec <- unique(df_cluster$Gene)
      
      # Create a data frame to store the fits and stats for each cluster
      df_cluster_fits <- data.frame()
      df_all_genes <- data.frame()
      
      # Create a list to store the Predicted_Values
      results[[cluster_name]][["Predicted_Values"]] <- list()
      
      for (gen in gene_vec) {
        # testing command
        # if (gen == "ABHD13") {
        #   stop("ABHD13 is the last gene")
        # }
        message(paste0("Fitting Hill equation for Gene ", gen))
        # Filter the data for the specific gene and order by gene name
        df_gene  <- df_cluster[df_cluster$Gene == gen, ]
        # Obtain the individuals for the gene
        ind_vec  <- unique(df_gene$Individuo)
        
        ### Set the constrains for the Hill equation
        # lower_amp <- 0
        # upper_amp <- max(df_gene[,2:5])* 1.5
        # 
        # lower_h   <- 0.1
        # upper_h   <- 50
        # 
        # upper_Tao <- 80
        # lower_Tao <- 2
        
        for (ind in ind_vec) {
          
          # print(paste0("Analyzing Gene: ", gen, " Individual: ", ind))
          message(paste0("Fitting Hill equation for Individual ", ind))
          
          # Filter the data for the specific individual
          df_ind <- df_gene[df_gene$Individuo == ind, ]
          
          # Transform the data to long format
          df_long <- df_ind %>%
            pivot_longer(
              cols         = starts_with("Exp_"),
              names_to     = "Time",
              names_prefix = "Exp_",
              values_to    = "Expression"
            ) %>%
            mutate(Time    = as.numeric(gsub("Exp_|h", "", Time)))
          
          #############################
          ##  Fit the Hill Equation  ##
          #############################
          
          # Compute the half maximum
          # half_max <- max(df_long$Expression, na.rm = TRUE) / 2
          
          # Set initial values for the parameters
          init_vals <- list(
            Amp  = as.numeric(df_long$Expression[which.max(abs(df_long$Expression))]), # max expression value
            h    = 1,
            Tao  = 72/2 # half of the time course (36h)
          )
          
          fit <- tryCatch({nlsLM(
            Expression ~ (Amp * Time^h) / (Time^h + Tao^h),
            data = df_long,
            start = list(
              Amp = init_vals$Amp,
              h   = init_vals$h,
              Tao = init_vals$Tao),
            control = nls.lm.control(maxiter = 200)
          ) }, error = function(e) return(NULL))
          
          if (!is.null(fit)) {
            # If the fit is successful, extract the coefficients
            data_plot <- df_long %>%
              mutate(Fit = predict(fit, newdata = df_long))
            
            # Compute R-squared
            # SS_res <- sum((data_plot$Expression - data_plot$Fit)^2)
            
            Mean   <- mean(data_plot$Expression)
            SS_res <- sum(resid(fit)^2)
            SS_tot <- sum((data_plot$Expression - Mean)^2)
            R2     <- 1 - SS_res / SS_tot
            
            # Create a fit data frame from 0 to 72h
            new_time <- data.frame(Time = seq(0, 72, by = 1))
            new_time$fit <- predict(fit, newdata = new_time)
            
            # Create a key for the plot
            plot_key <- paste0(gen, "_", ind)
            colnames(new_time)[2] <- paste0("Fit","_",plot_key)
            
            # Save the predicted values in the results list
            results[[cluster_name]][["Predicted_Values"]][[plot_key]] <- new_time
            
            # Extract the coefficients and summary
            coef_vals   <- coef(fit)
            fit_summary <- summary(fit)
            
            # Create a data frame to store the fit results
            df_fit <- data.frame(
              Cluster   = clus,
              Gene      = gen,
              Individuo = ind,
              Amp = coef_vals["Amp"],
              h   = coef_vals["h"],
              Tao = coef_vals["Tao"],
              RSS = deviance(fit),
              R2  = R2,
              Mean_Expression = Mean,
              StdError_Amp = fit_summary$coefficients["Amp", "Std. Error"],
              StdError_h   = fit_summary$coefficients[ "h" , "Std. Error"],
              StdError_Tao = fit_summary$coefficients["Tao", "Std. Error"],
              p_value_Amp  = fit_summary$coefficients["Amp", "Pr(>|t|)"],
              p_value_h    = fit_summary$coefficients[ "h" , "Pr(>|t|)"],
              p_value_Tao  = fit_summary$coefficients["Tao", "Pr(>|t|)"],
              exp_2  = data_plot$Expression[1],
              exp_20 = data_plot$Expression[2],
              exp_48 = data_plot$Expression[3],
              exp_72 = data_plot$Expression[4],
              fit_2  = data_plot$Fit[1],
              fit_20 = data_plot$Fit[2],
              fit_48 = data_plot$Fit[3],
              fit_72 = data_plot$Fit[4]
              )
            
          } else {
            # If the fit fails, create a data frame with NA values
            df_fit <- data.frame(
              Cluster   = clus,
              Gene      = gen,
              Individuo = ind,
              Amp = NA,
              h   = NA,
              Tao = NA,
              RSS = NA,
              R2  = NA,
              Mean_Expression = NA,
              StdError_Amp = NA,
              StdError_h   = NA,
              StdError_Tao = NA,
              p_value_Amp  = NA,
              p_value_h    = NA,
              p_value_Tao  = NA,
              exp_2  = data_plot$Expression[1],
              exp_20 = data_plot$Expression[2],
              exp_48 = data_plot$Expression[3],
              exp_72 = data_plot$Expression[4],
              fit_2  = NA,
              fit_20 = NA,
              fit_48 = NA,
              fit_72 = NA)
          }
          # Create a data frame to store the fits for each cluster
          df_cluster_fits <- rbind(df_cluster_fits, df_fit)
          message("Done, passing to the next gene:", gen)
        }
      }
      # Add the fits data frame to the results list
      results[[cluster_name]][["Fits_df"]] <- df_cluster_fits
      pred_list <- results[[cluster_name]][["Predicted_Values"]]
      Time <- new_time$Time
      # Create a data frame to store the predicted values
      fits_df <- do.call(cbind, lapply(names(pred_list), function(k) {
        df <- pred_list[[k]][2]
      }))
      # Combine the time and predicted values into a single data frame
      final_pred_df <- cbind(Time = Time, fits_df)
      results[[cluster_name]][["Predicted_Values"]] <- final_pred_df
      write.table(df_cluster_fits, file.path(eqtl_input_dir, paste0("Hill_fits_wo_constrains_", cluster_name, ".txt")), sep = "\t")
    }
    # Save the results list
    saveRDS(results, file.path(eqtl_input_dir, paste0("Hill_fits_wo_constrains_",type_of_cluster_name[j],"_responses.rds")))
  }
  
})

