#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(stringr)
  library(biomaRt)
  library(parallel)
  
}


#################
##  Load data  ##
#################
main_wd        <- getwd()
input_dir      <- file.path("Analyses","Inputs")
output_dir     <- file.path("Analyses","Outputs")
GO_data_dir    <- file.path(input_dir,"004_GO_data")
enrichment_dir <- file.path(GO_data_dir,"002_Enrichment_lists")

### load the tables of BP, CC and MF extracted from ClueGO
BiologicalProcess_GO_table   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_BP_06.02.2025.txt"))
CellularComponent_GO_table   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_CC_06.02.2025.txt"))
Molecular_Functions_GO_table  <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_MF_06.02.2025.txt"))

human_ncbi           <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_ncbi.txt"))
human_ncbi           <- human_ncbi[,c("NCBI.GeneID","Symbol")] 

########################
##  Define Functions  ##
########################

#
#This function takes as input the background table that the previous function 
#outputs, and also the raw table of GO annotations. It then
#reformats the table in order to put the levels as a vector, and also 
#calculates the size of each enrichment term depending on 

metadata_GO_builder <- function(GO_table_input, GO_term_list, Ontology_group) {
  names_list <- list()
  levels_list <- list()
  GO_table <- GO_table_input
  
  for (i in 1:length(GO_table[,1])){
    
    #For each term, to rewrite the levels. We split each cell and store the 
    #resulting vector
    
    levels_vec <- GO_table[i,2]
    levels_vec <- str_split(levels_vec, ",")
    levels_vec <- levels_vec[[1]]
    
    levels_list[[i]] <- as.integer(levels_vec)
  }
  
  
  #We store the levels of the term
  
  GO_table$Levels <- levels_list
  GO_table <- GO_table[,c(1,3,5)]
  
  #And finally, we add a column for the sizes of the terms in the table
  
  size_vector <- c()
  for (i in 1:nrow(GO_table)){
    col_name <-
      term_name <- GO_table[i,3]
      size <- length(GO_term_list[[i]])
      size_vector <- c(size_vector, size)
  }
  
  GO_table$Size <- size_vector
  GO_table$Ontology_group <- Ontology_group
  
  colnames(GO_table) <- c("Ontology", "Category", "Levels", "Size", "Ontology_group")
  return(GO_table)
}

biological_process  <- readRDS(file.path(enrichment_dir,"BP_list_transformed.RDS"))
cellular_component  <- readRDS(file.path(enrichment_dir,"CC_list_transformed.RDS"))
molecular_functions <- readRDS(file.path(enrichment_dir,"MF_list_transformed.RDS"))

BP_metadata_matrix <- metadata_GO_builder(BiologicalProcess_GO_table, biological_process, Ontology_group = "Biological_process")
CC_metadata_matrix <- metadata_GO_builder(CellularComponent_GO_table, cellular_component, Ontology_group = "Cellular_component")
MF_metadata_matrix <- metadata_GO_builder(Molecular_Functions_GO_table, molecular_functions, Ontology_group = "Molecular_functions")

dir.create(file.path(GO_data_dir,"004_Metadata_matrix"),showWarnings = F)

saveRDS(BP_metadata_matrix, file=file.path(GO_data_dir,"004_Metadata_matrix","BP_metadata_matrix.RDS"))
saveRDS(CC_metadata_matrix, file=file.path(GO_data_dir,"004_Metadata_matrix","CC_metadata_matrix.RDS"))
saveRDS(MF_metadata_matrix, file=file.path(GO_data_dir,"004_Metadata_matrix","MF_metadata_matrix.RDS"))
