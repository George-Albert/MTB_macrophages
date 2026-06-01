

####################
##  Dependencies  ##
####################

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
input_dir <- resolve_inputs_dir(project_root)
outputs_root <- resolve_outputs_dir(project_root)

required_packages <- c("tidyverse", "ggpubr", "openxlsx", "readxl")
load_required_packages(required_packages)

#################
##  Functions  ##
#################

dcols <- function(df){data.frame(colnames(df))} 

read_excel_table <- function(path, sheet = 1) {
  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "xls")) {
    return(as.data.frame(readxl::read_excel(path, sheet = sheet)))
  }

  openxlsx::read.xlsx(path, sheet = sheet)
}

#################
##  Load Data  ##
#################

script_id  <- "010_cor_pvalues_GO.R"
output_dir <- file.path(outputs_root, script_id)

# Load the files from ClueGO
clue_GO_dir      <- file.path(outputs_root, "006_clueGO")
clueGO_tab_early <- read_excel_table(file.path(clue_GO_dir,"Early.Up-1.xls"), sheet = 1)
clueGO_tab_mid   <- read_excel_table(file.path(clue_GO_dir,"Mid.Up-1.xls"), sheet = 1)
clueGO_tab_late  <- read_excel_table(file.path(clue_GO_dir,"Late.Up-1.xls"), sheet = 1)

clueGO_list <- list(clueGO_tab_early=clueGO_tab_early,clueGO_tab_mid=clueGO_tab_mid,
                    clueGO_tab_late=clueGO_tab_late)

### We selected the 1st configuration 
config_dir         <- "1st_Configuration_by_stage"
previous_script_id <- "008_Enrichment_Analysis_by_cluster"
output_dir_prev    <- file.path(outputs_root, previous_script_id)

# Load the files from our results
my_GO_dir   <- file.path(output_dir_prev,"005_Enrichment_GO",config_dir)
early_up_1  <- read_excel_table(file.path(my_GO_dir,"Early.Up","Enrichment_GO_Early.Up.xlsx"), sheet = 1)
mid_up_1    <- read_excel_table(file.path(my_GO_dir,"Mid.Up","Enrichment_GO_Mid.Up.xlsx"), sheet = 1)
late_up_1   <- read_excel_table(file.path(my_GO_dir,"Late.Up","Enrichment_GO_Late.Up.xlsx"), sheet = 1)

colnames_vec <- c("GO_ID","description","OR","p","fdr","genes_in_query")

# We selected the columns from colnames_vec
early_up <- early_up_1[,colnames_vec]
mid_up   <- mid_up_1[,colnames_vec]
late_up  <- late_up_1[,colnames_vec]

list_my_GO <- list(early_up=early_up,mid_up=mid_up,late_up=late_up)

### Save list 

list_out   <- list()
vec <- c(1,2,3)

for (i in vec) {
  
  
  clueGO_tab <- clueGO_list[[i]] 
  clueGO_tab <- clueGO_tab[,c(1,2,4,5,12)]
  colnames(clueGO_tab) <- c("GO_ID","description","p","fdr",
                            "Associated.Genes.Found")
  # fdr_threshold <- 0.01
  x <- clueGO_tab
  y <- list_my_GO[[i]]
  
  merged_df <- merge(x, y, by = "description", suffixes = c("_cluego", "_mi_GO"))
  merged_df <- na.omit(merged_df)
  
  merged_df <- merged_df %>% 
    distinct(description, .keep_all = TRUE)
  
  merged_df$log_p_cluego <- -log10(merged_df$p_cluego)
  merged_df$log_p_mi_GO  <- -log10(merged_df$p_mi_GO)
  
  list_out[[i]] <-  merged_df
  names(list_out)[i] <- names(list_my_GO)[i]
  
  cor_test <- cor.test(merged_df$log_p_cluego, merged_df$log_p_mi_GO, method = "spearman")  
  print(cor_test)
  cor_pear_test <- cor.test(merged_df$log_p_cluego, merged_df$log_p_mi_GO, method = "pearson")  
  print(cor_pear_test)
  
  title_vec      <- c("Early.Up","Mid.Up","Late.Up")
  color_for_line <- c("#E41A1C","#377EB8","#4DAF4A")  # Colorblind-friendly
  
  cor_plt <- ggplot(merged_df, aes(x = log_p_cluego, y = log_p_mi_GO)) +
    geom_point(alpha = 0.6, color = color_for_line[i], size = 2) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", linetype = "dashed", se = TRUE) +
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
             method = "spearman", label.x = 3, label.y = max(merged_df$log_p_mi_GO) * 0.9, size = 4) +
    labs(title = title_vec[i],
         subtitle = paste("Num. of common GO terms:", dim(merged_df)[1]),
         x = expression(-log[10](p-value)~ClueGO),
         y = expression(-log[10](p-value)~mi_GO)) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.title = element_text(face = "bold"))
  cor_plt
  
  dir.create(file.path(output_dir,"001_Figures","006_Cor_GO"),recursive=T,showWarnings =F)
  filename <- file.path(output_dir,"001_Figures","006_Cor_GO",paste0("Correlation_plt_",title_vec[i],".pdf"))
  
  pdf(file =filename ,width=8,height=6)
  print(cor_plt)
  dev.off()
  

}
