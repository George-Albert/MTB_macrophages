#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(ggplot2)
  library(brms)
  library(openxlsx)
  library(ggpubr)
  library(bayesplot)
  library(cmdstanr)
}
# cmdstanr::rebuild_cmdstan()
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
           "exponential" = paste0("exponential(2)"))
  }
  
  sd_tao_string <- if (!is.null(prior_sd_tao)) {
    prior_sd_tao
  } else {
    switch(dist_sd_tao,
           "normal"    = paste0("normal(0, ", round(tao_sd_random, 2), ")"),
           "student_t" = paste0("student_t(3, 0, ", round(tao_sd_random, 2), ")"),
           "auto"      = paste0("normal(0, ", round(tao_sd_random, 2), ")"),
           "exponential" = paste0("exponential(2)"))
  }
  
  sd_h_string <- if (!is.null(prior_sd_h)) {
    prior_sd_h
  } else {
    switch(dist_sd_h,
           "normal"    = paste0("normal(0, ", round(h_sd_random, 2), ")"),
           "student_t" = paste0("student_t(3, 0, ", round(h_sd_random, 2), ")"),
           "auto"      = paste0("normal(0, ", round(h_sd_random, 2), ")"),
           "exponential" = paste0("exponential(2)"))
  }
  # ---------------------------
  # Crear priors para brms
  priors <- c(
    prior_string(amp_string, nlpar = "Amp", lb = lb_amp),
    prior_string(tao_string, nlpar = "Tao", lb = lb_tao, ub =ub_tao),
    prior_string(h_string,   nlpar = "h",   lb = lb_h),
    
    prior_string(sd_amp_string, class = "sd", group = "Individuo", nlpar = "Amp"),
    prior_string(sd_tao_string, class = "sd", group = "Individuo", nlpar = "Tao"),
    prior_string(sd_h_string,   class = "sd", group = "Individuo", nlpar = "h")
  )
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
init_func <- function(df, chain_id = 1, base_h = 1, gene_id = NULL, log_seed = TRUE,
                      fixed_amp = NULL, fixed_tao = NULL, fixed_h = NULL) {
  
  df <- df %>% filter(Gene == gene_id) %>% arrange(Gene)
  
  half_max <- max(df$Expression, na.rm = TRUE) / 2
  amp_base <- as.numeric(quantile(df$Expression, 0.9, na.rm = TRUE))
  tao_base <- df$time[which.min(abs(df$Expression - half_max))]
  
  Amp_val <- if (!is.null(fixed_amp)) fixed_amp else amp_base
  Tao_val <- if (!is.null(fixed_tao)) fixed_tao else tao_base
  h_val   <- if (!is.null(fixed_h))   fixed_h   else base_h
  
  if (!is.null(gene_id)) {
    
    combined_seed <- as.integer(abs(sum(utf8ToInt(gene_id))) + chain_id * 100)
    
    if (chain_id != 1) {
      set.seed(combined_seed)
      
      if (is.null(fixed_amp)) {
        Amp_val <- runif(1, amp_base * 0.9, amp_base * 1.1)
      }
      if (is.null(fixed_tao)) {
        Tao_val <- runif(1, tao_base * 0.9, tao_base * 1.1)
      }
      if (is.null(fixed_h)) {
        h_val <- runif(1, base_h * 0.9, base_h * 1.1)
      }
    }
    
    if (log_seed) {
      new_entry <- tibble(Gene = gene_id, Chain = chain_id, Seed = combined_seed)
      assign("init_seeds_log", bind_rows(get("init_seeds_log", envir = .GlobalEnv), new_entry), envir = .GlobalEnv)
    }
  }
  
  list(
    Amp   = Amp_val,
    Tao   = Tao_val,
    h     = h_val,
    b_Amp = Amp_val,
    b_Tao = Tao_val,
    b_h   = h_val
  )
}

# Hill equation formula for brms
hill_formula <- brmsformula(
  Expression ~ Amp * time^h / (Tao^h + time^h),
  Amp+Tao+h ~ 0 + 1 + (1 || Individuo), # This is to allow the model to estimate the parameters for each individual
  nl = TRUE
)


priors <- make_priors_from_gene(param_list = param_list,
                                dist_amp = "normal",
                                dist_tao = "normal",
                                dist_h   = "normal",
                                dist_sd_tao = "normal",
                                dist_sd_amp = "auto",
                                dist_sd_h   = "auto",
                                sigma_prior = "exponential(1)",
                                include_cor   = F,
                                include_sigma = F,
                                max_frac_sd_amp = 1,
                                max_frac_sd_tao = 1,
                                max_frac_sd_h   = 1,
                                prior_sd_amp = "normal(0,1)",
                                prior_sd_tao = "normal(0,1)",
                                prior_sd_h   = "normal(0,1)",
                                lb_tao = NA,
                                ub_tao = NA,
                                lb_h   = NA,
                                lb_amp = NA)

# Set initial values for the parameters
# Create a list of initial values for each chain
n_chains <- 4
init_list <- lapply(1:n_chains, function(c) init_func(df_gene, 
                                                      chain_id = c,
                                                      base_h   = 1,
                                                      log_seed = TRUE,
                                                      gene_id  = gene_id,
                                                      fixed_amp = NULL, 
                                                      fixed_tao = NULL, 
                                                      fixed_h   = NULL))

{chains <- 4
  iter <- 6000
  warmup <- 1000
  cores <- 5
  adapt_delta <- .9999
  max_treedepth <- 20
  stepsize <- 0.001}

total_transitions <- chains * (iter - warmup)

# Fit the model using brms
fit <- tryCatch({
  brm(
    formula = hill_formula,
    data    = df_gene,
    family  = gaussian(),
    prior   = priors,
    init    = init_list,
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
