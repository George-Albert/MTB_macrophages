library(bayesplot)

safe <- fit

posterior_array <- as.array(fit)  # array 3D: iteraciones x chains x parámetros

library(brms)
library(bayesplot)
library(patchwork)
library(ggplot2)

# Extraer muestras posteriores
posterior_array <- as.array(fit)

# Parámetros a mostrar
params <- c(
  "b_Amp_Intercept",
  "b_Tao_Intercept",
  "b_h_Intercept",
  "sd_Individuo__Amp_Intercept",
  "sd_Individuo__Tao_Intercept",
  "sd_Individuo__h_Intercept"
)
friendly_names <- c(
  "Amplitude (Population Mean)",
  "Tao (Population Mean)",
  "Hill coefficient h (Population Mean)",
  "Amplitude (SD across Individuals)",
  "Tao (SD across Individuals)",
  "Hill coefficient h (SD across Individuals)"
)

# Crear lista de gráficos
plots <- lapply(1:length(params), function(i) {
  p <- params[i]
  name <- friendly_names[i]
  
  # Densidad posterior
  dens <- mcmc_areas(posterior_array, pars = p, prob = 0.8, fill = "#69b3a2") +
    ggtitle(name) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    ) +
    labs(x = "Value", y = "Density")
  
  # Traza MCMC
  trace <- mcmc_trace(posterior_array, pars = p, color = "#404080") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    ) +
    labs(x = "Iteration", y = "Parameter Value")
  
  dens / trace  # Combina densidad arriba, traza abajo
})

# Combinar todos en layout 2x3
poster_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]])

# Mostrar
poster_plot

# Guardar en alta resolución para póster
ggsave("posterior_distributions_poster.png", poster_plot, width = 20, height = 12, dpi = 300)

