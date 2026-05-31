
#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(brms)
  library(xlsx)
  library(ggpubr)
  library(bayesplot)
  library(cmdstanr)
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

# Function to create prior strings
make_priors_from_gene <- function(param_list, 
                                  dist_amp = "normal",
                                  dist_tao = "normal",
                                  dist_h   = "normal",
                                  dist_sd_tao = "normal",
                                  dist_sd_amp = "auto",
                                  dist_sd_h   = "auto",
                                  sigma_prior = "exponential(2)",
                                  include_sd = TRUE,
                                  include_cor   = FALSE,
                                  include_sigma = TRUE,
                                  max_frac_sd_amp = 0.5,
                                  max_frac_sd_tao = 0.5,
                                  max_frac_sd_h   = 0.5,
                                  prior_sd_amp = NULL,
                                  prior_sd_tao = NULL,
                                  prior_sd_h   = NULL,
                                  lb_tao = NA,
                                  ub_tao = NA,
                                  lb_h   = NA,
                                  lb_amp = NA) {
  
  # Extraer vectores
  amp_vec         <- (as.numeric(param_list[["Amp"]]))
  tao_vec         <- as.numeric(param_list[["Tao"]])
  h_vec           <- (as.numeric(param_list[["h"]]))
  amp_vec_global  <- (as.numeric(param_list[["Amp_global"]]))
  tao_vec_global  <- as.numeric(param_list[["Tao_global"]])
  h_vec_global    <- (as.numeric(param_list[["h_global"]]))
  
  # ---------------------------
  # AMP
  amp_est_global  <- mean(amp_vec_global, na.rm = TRUE)
  amp_sd_global   <- sd(amp_vec_global, na.rm = TRUE)
  amp_est         <- mean(amp_vec, na.rm = TRUE)
  amp_sd          <- sd(amp_vec, na.rm = TRUE)
  amp_sd_random   <- min(sd(amp_vec - amp_est, na.rm = TRUE), max_frac_sd_amp * amp_sd_global)
  
  # ---------------------------
  # TAO
  tao_est_global  <- mean(tao_vec_global, na.rm = TRUE)
  tao_sd_global   <- sd(tao_vec_global, na.rm = TRUE)
  tao_est         <- mean(tao_vec, na.rm = TRUE)
  tao_sd          <- sd(tao_vec, na.rm = TRUE)
  tao_sd_random   <- min(sd(tao_vec - tao_est, na.rm = TRUE), max_frac_sd_tao * tao_sd_global)
  
  # ---------------------------
  # h
  h_est_global    <- mean(h_vec_global, na.rm = TRUE)
  h_sd_global     <- sd(h_vec_global, na.rm = TRUE)
  h_est           <- mean(h_vec, na.rm = TRUE)
  h_sd            <- sd(h_vec, na.rm = TRUE)
  h_sd_random     <- min(sd(h_vec - h_est, na.rm = TRUE), max_frac_sd_h * h_sd_global)
  
  # ---------------------------
  # Convertir a strings
  amp_string <- switch(dist_amp,
                       "normal"    = paste0("normal(", round(amp_est_global, 2), ", ", round(amp_sd_global, 2), ")"),
                       "student_t" = paste0("student_t(3, ", round(amp_est_global, 2), ", ", round(amp_sd_global, 2), ")"))
  
  tao_string <- switch(dist_tao,
                       "normal"    = paste0("normal(", round(tao_est_global, 2), ", ", round(tao_sd_global, 2), ")"),
                       "student_t" = paste0("student_t(3, ", round(tao_est_global, 2), ", ", round(tao_sd_global, 2), ")"))
  
  h_string   <- switch(dist_h,
                       "normal"    = paste0("normal(", round(h_est_global, 2), ", ", round(h_sd_global, 2), ")"),
                       "student_t" = paste0("student_t(3, ", round(h_est_global, 2), ", ", round(h_sd_global, 2), ")"),
                       "gamma"     = paste0("gamma(", round(h_est_global, 2), ", ", round(h_sd_global, 2), ")"))
  
  sd_amp_string <- if (!is.null(prior_sd_amp)) {
    prior_sd_amp
  } else {
    switch(dist_sd_amp,
           "normal"    = paste0("normal(0, ", round(amp_sd_random, 2), ")"),
           "student_t" = paste0("student_t(3, 0, ", round(amp_sd_random, 2), ")"),
           "auto"      = paste0("normal(0, ", round(amp_sd_random, 2), ")"),
           "exponential" = paste0("exponential(",round(1/amp_sd_random, 2), ")"))
  }
  
  sd_tao_string <- if (!is.null(prior_sd_tao)) {
    prior_sd_tao
  } else {
    switch(dist_sd_tao,
           "normal"    = paste0("normal(0, ", round(tao_sd_random, 2), ")"),
           "student_t" = paste0("student_t(3, 0, ", round(tao_sd_random, 2), ")"),
           "auto"      = paste0("normal(0, ", round(tao_sd_random, 2), ")"),
           "exponential" = paste0("exponential(",round(1/tao_sd_random, 2), ")"))
  }
  
  sd_h_string <- if (!is.null(prior_sd_h)) {
    prior_sd_h
  } else {
    switch(dist_sd_h,
           "normal"    = paste0("normal(0, ", round(h_sd_random, 2), ")"),
           "student_t" = paste0("student_t(3, 0, ", round(h_sd_random, 2), ")"),
           "auto"      = paste0("normal(0, ", round(h_sd_random, 2), ")"),
           "exponential" = paste0("exponential(",round(1/h_sd_random, 2), ")"))
  }
  # ---------------------------
  # Crear priors para brms
  if (include_sd){
    priors <- c(
      prior_string(amp_string, nlpar = "Amp", lb = lb_amp),
      prior_string(tao_string, nlpar = "Tao", lb = lb_tao, ub =ub_tao),
      prior_string(h_string,   nlpar = "h",   lb = lb_h),
      
    # Standard deviations
      prior_string(sd_amp_string, class = "sd", group = "Individuo", nlpar = "Amp"),
      prior_string(sd_tao_string, class = "sd", group = "Individuo", nlpar = "Tao"),
      prior_string(sd_h_string,   class = "sd", group = "Individuo", nlpar = "h")
    )
  } else {
    priors <- c(
      prior_string(amp_string, nlpar = "Amp", lb = lb_amp),
      prior_string(tao_string, nlpar = "Tao", lb = lb_tao, ub =ub_tao),
      prior_string(h_string,   nlpar = "h",   lb = lb_h)
    )
  }
  if (include_sigma){
    priors <- c(priors, 
                prior_string(sigma_prior, class = "sigma"))
  }
  
  if (include_cor) {
    priors <- c(priors, set_prior("lkj(2)", class = "cor", group = "Individuo"))
  }
  
  return(priors)
}

# Function to initialize parameters for the model
init_func <- function(param_list, chain_id = 1, base_h = 1, gene_id = NULL, log_seed = TRUE,
                      fixed_amp = NULL, fixed_tao = NULL, fixed_h = NULL) {
  
  Amp <- param_list[["Amp"]]
  Tao <- param_list[["Tao"]]
  h   <- param_list[["h"]]
  
  half_max <- max(Amp, na.rm = TRUE) / 2
  amp_base <- as.numeric(quantile(Amp, 0.9, na.rm = TRUE))
  tao_base <- mean(Tao, na.rm = TRUE)
  
  Amp_val <- if (!is.null(fixed_amp)) fixed_amp else amp_base
  Tao_val <- if (!is.null(fixed_tao)) fixed_tao else tao_base
  h_val   <- if (!is.null(fixed_h))   fixed_h   else base_h
  sd_amp <- sd(Amp, na.rm = TRUE)
  sd_tao <- sd(Tao, na.rm = TRUE)
  sd_h   <- sd(h,   na.rm = TRUE)
  
  if (!is.null(gene_id)) {
    
    combined_seed <- as.integer(abs(sum(utf8ToInt(gene_id))) + chain_id * 100)
    
    if (chain_id != 1) {
      set.seed(combined_seed)
      # Randomize parameters based on the base values
      if (is.null(fixed_amp)) {
        Amp_val <- runif(1, amp_base * 0.9, amp_base * 1.1)
      }
      if (is.null(fixed_tao)) {
        Tao_val <- runif(1, tao_base * 0.9, tao_base * 1.1)
      }
      if (is.null(fixed_h)) {
        h_val <- runif(1, base_h * 0.9, base_h * 1.1)
      }
      # Add randomization for standard deviations
      sd_amp <- runif(1, sd_amp * 0.9, sd_amp * 1.1)
      sd_tao <- runif(1, sd_tao * 0.9, sd_tao * 1.1)
      sd_h   <- runif(1, sd_h   * 0.9, sd_h   * 1.1)
    }
    
    if (log_seed) {
      new_entry <- tibble(Gene = gene_id, Chain = chain_id, Seed = combined_seed)
      assign("init_seeds_log", bind_rows(get("init_seeds_log", envir = .GlobalEnv), new_entry), envir = .GlobalEnv)
    }
  }
  
  list(
    b_Amp = Amp_val,
    b_Tao = Tao_val,
    b_h   = h_val,
    sd_Individuo__Amp = sd_amp,
    sd_Individuo__Tao = sd_tao,
    sd_Individuo__h   = sd_h
    
  )
}
#################
##  Load data  ##
#################

# Load the delta expression data
delta_expr_df <- read.table(file.path(eqtl_input_dir, "delta_expr_df.txt"),header = TRUE, sep = "\t")
# Load gene with minimum amount of individuals
min_ind_df <- read.table(file.path(eqtl_input_dir, "genes_r2_threshold_and_min_individuals.txt"), header = TRUE, sep = "\t")
min_ind_df_list <- readRDS(file.path(eqtl_input_dir, "genes_r2_threshold_and_min_individuals.rds"))

selected_r2_df <-min_ind_df_list[[11]]
summary(selected_r2_df)
# 
# lower_Tao <- min(selected_r2_df$Tao, na.rm = TRUE)
# upper_Tao <- max(selected_r2_df$Tao, na.rm = TRUE)
# lower_h   <- min(selected_r2_df$h, na.rm = TRUE)
# upper_h   <- max(selected_r2_df$h, na.rm = TRUE)
# epsilon <- 0.5
# 
# selected_r2_df <- selected_r2_df %>%
#   filter(
#     Tao > lower_Tao + epsilon,
#     Tao < upper_Tao - epsilon,
#     h   > lower_h   + epsilon,
#     h   < upper_h   - epsilon
#   )

ggplot(selected_r2_df, aes(x = Gene, y = h)) +
  geom_point() +
  labs(title = "Scatter plot of Tao vs Amp") +
  theme_minimal()

set.seed(123) # Set seed for reproducibility

h_clusters <- kmeans(selected_r2_df$h, centers = 3, nstart = 20)
# Add cluster information to the data frame
selected_r2_df$h_cluster <- as.factor(h_clusters$cluster)

# Plot the distribution of h clusters

hist_h <- ggplot(selected_r2_df, aes(x = h, fill = h_cluster)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.5) +
  labs(title = "Distribution of h clusters after filtering",
       x = "h", y = "Frequency") +
  theme_minimal()

### Lets inspect the clusters

length(which(selected_r2_df$h_cluster == 1)) # 1 cluster)
# 2228
length(which(selected_r2_df$h_cluster == 2)) # 2 cluster
# 32939
length(which(selected_r2_df$h_cluster == 3)) # 3 cluster
# 5273

clus2_df <- selected_r2_df[selected_r2_df$h_cluster == 2,]
# Load the gene parameters
clus2_df <- read.table(file.path(eqtl_input_dir, "Cluster_2_h_stats.txt"), sep = "\t")
stat_df  <- read.table(file.path(eqtl_input_dir, "Tao_filt_stats.txt"), sep = "\t")

df_cluster <- clus2_df
# df_cluster <- selected_r2_df
summary(df_cluster)

# Extract the Amp, Tao, and h dataframes from the list

df_amp <- df_cluster %>% 
  select(Gene, Individuo, Amp)
df_tao <- df_cluster %>%
  select(Gene, Individuo, Tao)
df_h   <- df_cluster %>%
  select(Gene, Individuo, h)

# Extract the cluster of interest in the expression data
number_of_cluster <- 3 
data <- delta_expr_df %>% 
  filter(Cluster == number_of_cluster) 

# Put it in long format
df_delta_expr_long <- data %>%
  pivot_longer(
    cols = starts_with("Exp_"),
    names_to = "Time",
    values_to = "Expression"
  ) %>%
  mutate(
    time      = as.numeric(gsub("Exp_|h", "", Time)),
    Individuo = factor(Individuo),
    Gene      = factor(Gene)
  ) %>%
  select(Gene, Individuo, time, Expression)
# Check the structure of the data
df_delta_expr_long <- df_delta_expr_long %>% 
  mutate(
    Individuo = as.factor(Individuo),
    Gene = as.factor(Gene)
  )
str(df_delta_expr_long)

# Create a list to store results for each gene
results_by_gene <- list()

all_genes <- as.character(unique(data$Gene))
# param_name <- c("Amp", "Tao", "h")
# Define hill equation function

hill_formula <- brmsformula(
  Expression ~ Amp * time^h / (Tao^h + time^h),
  Amp ~  0 + 1 + (1 |p| Individuo),      # here we can add (1 | Individuo)
  Tao ~  0 + 1 + (1 |p| Individuo),
  h   ~  0 + 1 + (1 | Individuo),
  # Amp+Tao+h  ~ 0 + 1 + (1 | Individuo), # This is to allow the model to estimate the parameters for each individual
  nl = TRUE
)
# Amp ~ 0 + 1 + (1 | Individuo),      # here we can add (1 | Individuo)
# Tao ~ 0 + 1 + (1 | Individuo),
# h   ~ 0 + 1 + (1 | Individuo),

# Initialization log for seeds
init_seeds_log <- tibble()

already_fitted <- list.files(path = file.path(eqtl_input_dir, "brms_fits"),
                             pattern = "^fit_.*\\.rds$") |>
  gsub("fit_|\\.rds", "", x = _)

for (gene_id in all_genes) {
  
  if (gene_id %in% already_fitted) {
    message("Skipping already fitted gene: ", gene_id)
    next
  }
  message("Fitting gene: ", gene_id)
  
  # Filter data for the current gene
  df_gene <- df_delta_expr_long %>% filter(Gene == gene_id)
  
  # Next if all values are NA (can´t fit a model)
  if (all(is.na(df_gene$Expression))) {
    warning("All NA for gene ", gene_id, " – skipping.")
    next
  }
  
  df_amp_filt <- df_amp %>%
    filter(Gene == gene_id) %>%
    # mutate(logAmp = log(Amp / max_expr)) %>% 
    arrange(Gene)
  
  # ggdensity(df_amp_filt$logAmp, main = paste("Density of Amp for", gene_id))
  
  df_tao_filt <- df_tao %>%
    filter(Gene == gene_id) %>%
    # mutate(logTao = log(Tao /10)) %>%
    arrange(Gene)
 
   # ggdensity(df_tao_filt$Tao, main = paste("Density of Tao for", gene_id))
  
  df_h_filt <- df_h %>%
    filter(Gene == gene_id) %>%
    mutate(logh = log(h)) %>%
    arrange(Gene)
 
   # ggdensity(df_h_filt$h, main = paste("Density of h for", gene_id))
  
  cat("There are ", 
      length(which(unique(df_amp_filt$Gene) != unique(df_tao_filt$Gene))),
      " differences between Amp and Tao values.")
  cat("There are ", 
      length(which(unique(df_amp_filt$Gene) != unique(df_h_filt$Gene))),
      " differences between Amp and h values.")
  
  param_list <- list(
    Amp_global= df_amp_filt$Amp,
    Tao_global= df_tao_filt$Tao,
    h_global  = df_h_filt$h,
    
    Amp = df_amp_filt$Amp,
    Tao = df_tao_filt$Tao,
    h   = df_h_filt$h
  )
  
  # Generar priors from the parameters
  priors <- make_priors_from_gene(param_list = param_list,
                                  dist_amp = "normal",
                                  dist_tao = "normal",
                                  dist_h   = "normal",
                                  dist_sd_tao = "normal",
                                  dist_sd_amp = "normal",
                                  dist_sd_h   = "normal",
                                  sigma_prior = "exponential(1)",
                                  include_sd  = T,
                                  include_cor   = F,
                                  include_sigma = T,
                                  max_frac_sd_amp = 1,
                                  max_frac_sd_tao = 1,
                                  max_frac_sd_h   = 1,
                                  prior_sd_amp = NULL,
                                  prior_sd_tao = NULL,
                                  prior_sd_h   = NULL,
                                  lb_tao = 2,
                                  ub_tao = 80,
                                  lb_h   = 0.1,
                                  lb_amp = 0)
  
  # Set initial values for the parameters
  # Create a list of initial values for each chain
  n_chains <- 4
  init_list <- lapply(1:n_chains, function(c) init_func(param_list, 
                                                        chain_id = c,
                                                        base_h   = 1,
                                                        log_seed = TRUE,
                                                        gene_id  = gene_id,
                                                        fixed_amp = NULL, 
                                                        fixed_tao = NULL, 
                                                        fixed_h   = NULL))
  
  {
  chains <- n_chains
  iter   <- 6000
  warmup <- 3000
  cores  <- 6
  adapt_delta   <- 0.9999
  max_treedepth <- 20
  stepsize <- 0.001
  }
  
  total_transitions <- chains * (iter - warmup)

  # Fit the model using brms
  fit <- tryCatch({
    brm(
      formula = hill_formula,
      data    = df_gene,
      family  = gaussian(),
      prior   = priors,
      init    = init_list,
      # sample_prior = "yes",
      chains  = chains,
      iter    = iter,
      cores   = cores,
      warmup  = warmup,
      control = list(adapt_delta   = adapt_delta,
                     max_treedepth = max_treedepth,
                     step_size      = stepsize),
      seed    = 123,
      backend = "cmdstanr"
    )
  }, error = function(e) {
    message("Error fitting ", gene_id, ": ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) {
    results_by_gene[[gene_id]] <- list(model = NULL, imputed = NULL)
    next
  }
  # Extract the fitted values
  ind_coefs <- coef(fit)$Individuo
  coefs <- data.frame(
    Gene= gene_id,
    Individuo = rownames(ind_coefs),
    Amp    = ind_coefs[, "Estimate" , "Amp_Intercept"],
    Amp_sd = ind_coefs[, "Est.Error", "Amp_Intercept"],
    Tao    = ind_coefs[, "Estimate" , "Tao_Intercept"],
    Tao_sd = ind_coefs[, "Est.Error", "Tao_Intercept"],
    h      = ind_coefs[, "Estimate" , "h_Intercept"],
    h_sd   = ind_coefs[, "Est.Error", "h_Intercept"]
    
  )
  
  results_by_gene[[gene_id]] <- coefs  
  
  # Store the results in the list
  saveRDS(fit, file = file.path(eqtl_input_dir,"brms_fits", paste0("fit_", gene_id, ".rds")))
  
}
# Combine all results into a single data frame

final_results <- bind_rows(results_by_gene)

write.table(final_results, 
            file.path(eqtl_input_dir , "001_bmrs_fit_results_scaled.txt"), 
            sep = "\t", quote = FALSE)
write.csv(final_results, 
          file.path(eqtl_input_dir , "001_bmrs_fit_results_scaled.csv"), 
          row.names = FALSE)

# Carpeta con los fits guardados
fit_dir <- file.path(eqtl_input_dir, "brms_fits")

# Obtener los nombres de archivos
fit_files <- list.files(fit_dir, pattern = "^fit_.*\\.rds$", full.names = TRUE)

# Inicializar tabla
divergences_table <- data.frame(Gene = character(), Divergences = integer())

for (file in fit_files) {
  fit <- readRDS(file)
  gene_id <- gsub("fit_|\\.rds", "", basename(file))
  
  # Extraer divergencias
  np <- nuts_params(fit)
  divs <- sum(np$Parameter == "divergent__" & np$Value == 1)
  
  divergences_table <- rbind(divergences_table,
                             data.frame(Gene = gene_id, Divergences = divs))
}

# Save the divergences table

write.table(divergences_table, 
            file.path(eqtl_input_dir, "002_bmrs_divergences.txt"), 
            sep = "\t")
write.csv(divergences_table, 
          file.path(eqtl_input_dir, "002_bmrs_divergences.csv"), 
          row.names = FALSE)
  
  
  
  
  
  
  
  
  

  posterior_summary(fit, variable = c("b_Amp_Intercept", "b_Tao_Intercept", "b_h_Intercept"))
  
  # Impute the missing data
  df_na <- df_gene %>% filter(is.na(Expression))
  if (nrow(df_na) > 0) {
    imputed_samples <- posterior_predict(fit, newdata = df_na)
    df_na$Expression_imputed <- colMeans(imputed_samples)
  } else {
    df_na <- NULL
  }
  
  # Extract the fitted values
  ind_coefs <- coef(fit)$Individuo
  coefs <- data.frame(
    Gene= gene_id,
    Individuo = rownames(ind_coefs),
    Amp    = ind_coefs[, "Estimate", "Amp_Intercept"],
    Amp_sd = ind_coefs[, "Est.Error", "Amp_Intercept"],
    Tao    = ind_coefs[, "Estimate", "Tao_Intercept"],
    Tao_sd = ind_coefs[, "Est.Error", "Tao_Intercept"],
    h      = ind_coefs[, "Estimate", "h_Intercept"],
    h_sd   = ind_coefs[, "Est.Error", "h_Intercept"]
    
  )
  
  # Save the results
  results_by_gene[[gene_id]] <- list(
    model   = fit,
    imputed = df_na,
    coefs   = coefs
    )
}




df_gene_id <- df_cluster %>%
  filter(Gene == gene_id) %>%
  arrange(Gene)
df_inspect <- df_gene_id %>%
  left_join(coefs, by = c("Gene","Individuo"), suffix = c(".fit", ".bmr"))

# Correlation between fitted and brms values

cor <- cor.test(df_inspect$Amp.fit, df_inspect$Amp.bmr, method = "pearson")
cor_plot <- ggplot(df_inspect, aes(x = Amp.fit, y = Amp.bmr)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(
    title = paste0("R = ", round(cor$estimate, 3),
                   ", p.value = ", signif(cor$p.value, 3),
                   ", ",nrow(df_inspect)," Individuals"),
    x = "Fitted Amp",
    y = "brms Amp") +
  theme_minimal()

# Plot the Fits
hill_out_dir <- file.path(output_dir,"001_Figures","010_Hill_bmr_fit")
dir.create(hill_out_dir, showWarnings = FALSE)

plot_width <- 4
plot_height <- 2.5
id_gene <- unique(df_gene$Gene)

fit_stats <- list()  # Para guardar resultados por individuo

for (g in id_gene) {
  
  # Subconjuntos por gen
  df_real_gene <- df_gene %>% filter(Gene == g)
  df_param_gene <- coefs %>% filter(Gene == g)
  
  if (nrow(df_real_gene) == 0 || nrow(df_param_gene) == 0) next
  
  # --- Predicciones por individuo ---
  df_pred_gene <- df_param_gene %>%
    rowwise() %>%
    mutate(data = list({
      t <- df_real_gene %>% filter(Individuo == Individuo) %>% pull(time)
      y_obs <- df_real_gene %>% filter(Individuo == Individuo) %>% pull(Expression)
      y_pred <- (Amp * t^h) / (Tao^h + t^h)
      
      rss <- sum((y_obs - y_pred)^2, na.rm = TRUE)
      tss <- sum((y_obs - mean(y_obs))^2, na.rm = TRUE)
      r2  <- 1 - rss / tss
      
      fit_stats[[paste(g, Individuo, sep = "_")]] <<- tibble(
        Gene = g,
        Individuo = Individuo,
        RSS = rss,
        R2 = r2,
        Amp = Amp,
        Tao = Tao,
        h = h
      )
      
      tibble(Time = t, Predicted = y_pred,RSS = rss, R2 = r2)
    })) %>%
    unnest(cols = c(data)) %>%
    mutate(Gene = g, Setup = Individuo)
  
  fit_stats_df <- bind_rows(fit_stats)
  
  df_real_gene <- df_real_gene %>%
    mutate(Setup = Individuo)
  
  # Layout del facet
  setups_gene <- unique(df_pred_gene$Setup)
  n_panels <- length(setups_gene)
  ncol <- ceiling(sqrt(n_panels))
  nrow <- ceiling(n_panels / ncol)
  total_width <- ncol * plot_width
  total_height <- nrow * plot_height
  
  message(sprintf("Plotting gene %s with layout %dx%d", g, ncol, nrow))
  
  # --- GRAFICO ---
  p <- ggplot() +
    geom_line(data = df_pred_gene, aes(x = Time, y = Predicted), color = "blue") +
    geom_point(data = df_real_gene, aes(x = time, y = Expression), color = "black", size = 1.5) +
    geom_point(data = df_real_gene %>% filter(time == 0), 
               aes(x = time, y = Expression), color = "orange", size = 2.5) +
    facet_wrap(~Setup, ncol = ncol, scales = "fixed") +
    labs(title = paste("Gene:", g),
         x = "Time (h)", y = "Expression") +
    scale_x_continuous(breaks = c(0, 2, 20, 48, 72)) +
    theme_bw() +
    theme(strip.text = element_text(size = 10),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold"))
  
  # Guardar PDF
  ggsave(
    filename = file.path(hill_out_dir, sprintf("%s.pdf", g)),
    plot = p,
    width = total_width,
    height = total_height,
    units = "in"
  )
}

# --- GUARDAR ESTADISTICAS ---
fit_stats_df <- bind_rows(fit_stats)
write.csv(fit_stats_df, file = file.path(hill_out_dir, "fit_stats_by_individual.csv"), row.names = FALSE)

