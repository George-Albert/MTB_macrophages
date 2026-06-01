#########################
## 0.Load Dependencies ##
#########################

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

# Standardized package loading keeps execution reproducible across environments.
required_packages <- c(
  "tidyverse", "ggplot2", "ggrepel", "limma", "edgeR", "qvalue", "cowplot",
  "stats", "openxlsx", "ngram", "RColorBrewer", "dendextend", "reshape2"
)
load_required_packages(required_packages)


#######################
##  Set directories  ##
#######################
script_id <- "001_Setting_exp_data"
output_dir <- file.path(outputs_root, script_id)


########################
##  Loading the data  ##
########################
reads     = read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
cols      = read.table(file.path(input_dir,"002_Processed","ready_for_DE","metadata.txt"))
cols_whole= read.table(file.path(input_dir,"002_Processed","whole","metadata_whole.txt"))

##############
##  Step#1  ##
##############
{ 
  ### Inspecting the data to arrive to conlcusions about the samples and timepoints we have
  
  total_number_of_individuals <- unique(cols$Individual)
  # 102 individuals in the in the processed metadata
  length(which(colnames(reads)!=cols$Sample))
  
  # Just comparing with the raw metadata
  df <- data.frame(setup=colnames(cols_whole))
  df <- df %>%
    extract(col = setup, into = c("Individuo", "Condition", "TimePoint"),
            regex = "^(.*)_(.*)_(.*)$")
  length(unique(df$Individuo))
  length(which(sort(unique(df$Individuo)) != sort(total_number_of_individuals)))
  # 102 Individuals in the raw metadata
  
  # Differential expression outputs are produced by script 002_Differential_expression.
  DE_out_dir <- file.path(outputs_root, "002_Differential_expression")
  results <- read.table(file.path(DE_out_dir,"003_DE_def","resultados.txt"))
  DE_tab  <- read.table(file.path(DE_out_dir,"003_DE_def","DE_table.txt"))
  
  ### Individuals per condition TB
  length(which(grepl("TB_72h",colnames(reads))))
  #87
  length(which(grepl("TB_48h",colnames(reads))))
  #93
  length(which(grepl("TB_20h",colnames(reads))))
  #95
  length(which(grepl("TB_02h",colnames(reads))))
  #99
  
  ### Individuals per condition NI
  length(which(grepl("NI_72h",colnames(reads))))
  #87
  length(which(grepl("NI_48h",colnames(reads))))
  #93
  length(which(grepl("NI_20h",colnames(reads))))
  #98
  length(which(grepl("NI_02h",colnames(reads))))
  #98
  
  ### Sum of individuals per condition
  sum(c(length(which(grepl("TB_72h",colnames(reads)))),
        length(which(grepl("TB_48h",colnames(reads)))),
        length(which(grepl("TB_20h",colnames(reads)))),
        length(which(grepl("TB_02h",colnames(reads))))))
  # 374
  sum(c(length(which(grepl("NI_72h",colnames(reads)))),
        length(which(grepl("NI_48h",colnames(reads)))),
        length(which(grepl("NI_20h",colnames(reads)))),
        length(which(grepl("NI_02h",colnames(reads))))))
  # 376     
  
  # Count how many samples there are per individual, treatment and time
  sample_exp <- cols %>%
    group_by(TimePoint) %>%
    summarise(n = n(), .groups = "drop")
  
  print(sample_exp)
  
  # TimePoint     n
  # <chr>     <int>
  # 1 02h         197
  # 2 20h         193
  # 3 48h         186
  # 4 72h         174
  
  # Count how many samples there are per timepoint, treatment 
  sample_exp1 <- cols %>%
    group_by(TimePoint,Treatment) %>%
    summarise(n = n(), .groups = "drop")
  
  print(sample_exp1)
  # TimePoint Treatment     n
  # <chr>     <chr>     <int>
  # 1 02h       NI           98
  # 2 02h       TB           99
  # 3 20h       NI           98
  # 4 20h       TB           95 # Unpair by Time Point
  # 5 48h       NI           93 
  # 6 48h       TB           93
  # 7 72h       NI           87
  # 8 72h       TB           87
  
  # Create a wide format data frame with TB and NI counts
  sample_pairs <- cols %>%
    group_by(TimePoint,Treatment,Individual) %>%
    summarise(n = n(), .groups = "drop") %>% 
    pivot_wider(names_from = Treatment, values_from = n, values_fill = 0) %>% 
    arrange(Individual, TimePoint) 
  
  # Keep only the rows with both TB and NI counts
  valid_pairs <- sample_pairs %>% filter(TB > 0 & NI > 0)
  # Identify the inconsistent pairs
  incomplete_pairs <- sample_pairs %>% filter(TB == 0 | NI == 0)
  
  ### Inspect the incomplete pairs data
  table(incomplete_pairs$TimePoint)
  # 02h 20h 48h 72h 
  # 1   5   8  16 
  table(incomplete_pairs$NI)
  # 0  1 
  # 14 16
  table(incomplete_pairs$TB)
  # 0  1 
  # 16 14
  table(incomplete_pairs$Individual)
  # I_15   I_22    I_3   I_33    I_4   I_42   I_46   I_48 I_C053  I_C12  I_C13  I_C14  I_C16  I_C32   I_C4 
  # 1      1      1      1      2      2      1      1      1      1      3      1      1      2      1 
  # I_C44   I_C5  I_C50  I_C59  I_C60 
  # 1      3      1      2      3 
  
  
  #Count how many combinations each individual has
  ind_time_trt_counts <- cols %>%
    distinct(Individual, TimePoint, Treatment) %>%
    count(Individual) %>% 
    arrange(n)
  
  #Filter individuals with at least 3 combinations in NI and 4 in TB
  sample_pairs2 <- sample_pairs %>%
    pivot_longer(cols = c("NI", "TB"), names_to = "Condition", values_to = "Count") %>%
    # filter(Count > 0) %>%  # only keep rows with non-zero counts
    group_by(Individual, Condition) %>%
    summarise(Samples = sum(Count), .groups = "drop") %>% 
    pivot_wider(names_from = Condition, values_from = Samples, values_fill = 0) %>% 
    filter((TB == 4 & NI >= 3) | (TB >= 3 & NI == 4)) %>% 
    arrange(NI)
  dim(sample_pairs2)[1]
  
  # 72 individuals with at least 4 combinations in TB and 3 in NI
  
  ind_to_impute <- sample_pairs2 %>%
    filter(NI == 3) %>%
    select(Individual)
  # 8 Individuals with 3 combinations in NI to impute
  
  # Filter sample_pairs2 to keep only the individuals to impute
  sample_pairs_to_impute <- sample_pairs %>%
    filter(Individual %in% ind_to_impute$Individual &
             NI == 0) %>%
    # select(Individual, TimePoint) %>%
    arrange(Individual, TimePoint)
  
  # Save the data frames
  write.table(ind_time_trt_counts, file.path(input_dir,"006_EQTL_data","Number_of_comb_per_ind.txt"),
              sep = "\t")
  write.xlsx(ind_time_trt_counts, file.path(input_dir,"006_EQTL_data","Number_of_comb_per_ind.xlsx"))
  write.xlsx(sample_pairs2, file.path(input_dir,"006_EQTL_data","Number_of_comb_per_ind_filter.xlsx"))
  write.xlsx(sample_pairs_to_impute, file.path(input_dir,"006_EQTL_data","Number_of_comb_per_ind_to_impute.xlsx"))
  
  # Filter individuals with less than 8 combinations (4 timepoints × 2 treatments)
  ind_incomplete <- ind_time_trt_counts %>%
    filter(n < 8)
  
  length(unique(ind_incomplete$Individual))
  # 30 individuals
  
  # All the possible combinations of individuals, timepoints and treatments
  combinaciones_completas <- expand.grid(
    Individual = unique(cols$Individual),
    TimePoint = unique(cols$TimePoint),
    Treatment = c("TB", "NI")
  )
  
  # Compare with the data we have
  combinaciones_faltantes <- combinaciones_completas %>%
    anti_join(cols %>% distinct(Individual, TimePoint, Treatment),
              by = c("Individual", "TimePoint", "Treatment"))
  
  # Filter all the individuals that have all the combinations
  ind_completos <- ind_time_trt_counts %>% filter(n == 8) %>% pull(Individual)
  length(which(ind_time_trt_counts$n==8))
  length(ind_completos)
  # 72
  
  #############################################################
  ##    ##
  #############################################################
  
  # Filter cols df (metadata)
  cols_filtrado <- cols %>% filter(Individual %in% ind_completos)
  # cols_filtrado_to_impute <- cols %>% filter(Individual %in% sample_pairs2$Individual)
  
  # We are going doble check the batch again just make sure there are not duplicates
  inspect_batch <- cols_filtrado %>%
    group_by(Individual, Treatment,Batch) %>%
    summarise(n = n(), .groups = "drop") %>% 
    pivot_wider(names_from = Treatment, values_from = n, values_fill = 0) %>% 
    arrange(Batch, Individual)
  
  ind_per_batch <- cols %>%
    select(Individual, Batch) %>%
    group_by(Batch) %>% 
    unique() %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(Batch)
  
  write.table(ind_per_batch, file.path(input_dir,"006_EQTL_data","ind_per_batch.txt"))
  write.xlsx(ind_per_batch, file.path(input_dir,"006_EQTL_data","ind_per_batch.xlsx"))
  
  
}

##############
##  Step#2  ##
##############
{
  ### Do the DE analysis in order to get the response to infection data
  
  # Do the normalization 
  design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols_filtrado)
  
  # design=model.matrix(~1)
  dge <- DGEList(counts=reads[,cols_filtrado$Sample])
  dge <- calcNormFactors(dge,method="TMM")
  v <- voom(dge,plot=T)
  exp <- v$E
  exp <- as.data.frame(exp)
  
  ###########################################################
  ##  Log Fold Changes. Response to infection across time  ##
  ###########################################################
  
  cols <- cols_filtrado
  ###Create a matrix to store the results
  all_individuals_data <- list()
  time_vec <- unique(cols$TimePoint)
  
  ind_data <- data.frame(Gene = rownames(reads),
                         Exp_02h = 0,
                         Exp_20h = 0,
                         Exp_48h = 0,
                         Exp_72h = 0)
  ### Compute the difference between TB and NI for each individual and time point
  for (ind in unique(cols$Individual)) {
    # Subset by  Individual
    print(paste("The individual is: ", ind))
    
    for (time in time_vec) {
      
      # Get the column names for TB and NI for the current individual and time point
      tb_col <- cols$Sample[cols$Individual == ind & cols$TimePoint == time &
                              cols$Treatment == "TB"]
      print(paste("The TB column is: ", tb_col))
      ni_col <- cols$Sample[cols$Individual == ind & cols$TimePoint == time &
                              cols$Treatment == "NI"]
      print(paste("The NI column is: ", ni_col))
      # Check if the columns exist
      if (length(tb_col) == 0 || length(ni_col) == 0) {
        print(paste("No TB or NI columns found for individual: ", ind, " and time: ", time,". you need to impute"))
        exp[, grepl(ind, colnames(exp))]
        # ind_data[[paste0("Exp_", time)]] <- NA
        next
      } else if (length(tb_col) == 1 & length(ni_col) == 1) {
        
        print(paste("Resting TB and NI columns for individual: ", ind, " and time: ", time))
        print(paste("Executing the rest", tb_col,"-", ni_col))
        
        ind_data[[paste0("Exp_",time)]] <- exp[, tb_col] - exp[, ni_col]
        head(ind_data,n=5)
      }
    }
    
    ind_data$Individuo <- ind
    all_individuals_data[[ind]] <- ind_data
  }
  
}

##############
##  Step#3  ##
##############
{
  
  ### Create the response to infection dataframe with the cluster information we have and the DE genes
  # We use do.call(rbind, all_individuals_data) to combine the list of data 
  # frames into a single data frame
  delta_expr_df <- do.call(rbind, all_individuals_data)
  
  # Load lfc with clustering information
  cluster_input_dir <- file.path(input_dir, "003_Clustering_tabs")
  cluster_out <- read.table(file.path(cluster_input_dir,"lfc_data_with_clusters_information.txt"))
  
  # We selected the k=16 cluster results
  cluster_out1 <- cluster_out[,c(1:4,6)]
  
  # # Reorder columns
  # delta_expr_df1 <- delta_expr_df[, c("Gene", unique(cols$TimePoint))]
  
  # Add the cluster information to the delta_expr_df
  delta_expr_df$Cluster <- cluster_out1$k_16[match(delta_expr_df$Gene, rownames(cluster_out1))]
  
  delta_expr_df <- delta_expr_df %>% 
    filter(!is.na(Cluster))
  
  dim(delta_expr_df)[1]/dim(cluster_out1)[1]
  # 72 Idividuals
  
  # Save the data frame and the list to a file
  eqtl_input_dir <- file.path(input_dir,"006_EQTL_data")
  dir.create(eqtl_input_dir, showWarnings = FALSE)
  
  saveRDS(all_individuals_data, file.path(eqtl_input_dir, "all_individuals_exp_data.rds"))
  write.table(delta_expr_df, file.path(eqtl_input_dir, "delta_expr_df.txt"), sep = "\t", row.names = FALSE,
              quote = FALSE)
  
  
}


