bootstrap_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
repeat {
  setup_candidate <- file.path(bootstrap_root, "Codes", "_shared", "project_setup.R")
  if (file.exists(setup_candidate)) {
    source(setup_candidate)
    break
  }
  parent <- dirname(bootstrap_root)
  if (identical(parent, bootstrap_root)) {
    stop("Could not locate Codes/_shared/project_setup.R")
  }
  bootstrap_root <- parent
}

project_root <- find_project_root()
setwd(project_root)

required_packages <- c(
  "tidyverse", "ggplot2", "ggrepel", "limma", "edgeR", "qvalue",
  "cowplot", "openxlsx", "ngram", "RColorBrewer", "dendextend",
  "reshape2", "xlsx", "minpack.lm", "broom", "patchwork", "ggforce",
  "ggpubr", "shiny", "brms", "bayesplot", "cmdstanr", "Rtsne",
  "HGNChelper", "ggnewscale", "biomaRt", "rtracklayer", "tidyr",
  "dplyr", "gridExtra", "ggthemes", "ggdendro", "dendsort", "umap",
  "circlize", "ComplexHeatmap", "ConsensusClusterPlus", "factoextra",
  "NbClust", "mclust", "clustertend", "stringr", "babelgene",
  "org.Hs.eg.db", "igraph", "ggsci"
)

assert_required_packages(required_packages)

message("Preflight OK: all required packages are available.")
