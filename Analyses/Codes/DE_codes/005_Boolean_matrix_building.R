#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(stringr)
  library(biomaRt)
  library(parallel)
  library(dplyr)
  library(babelgene)
  library(HGNChelper)
  library(xlsx)
}


#################
##  Load data  ##
#################
main_wd     <- getwd()
input_dir   <- file.path("Analyses","Inputs")
output_dir  <- file.path("Analyses","Outputs")
GO_data_dir <- file.path(input_dir,"004_GO_data")
#Date of GO tables version: 06.02.2025

dir.create(file.path(GO_data_dir,"001_GO_table_to_list"),showWarnings = F)

### Load any table with the genes names
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
genes_names <- rownames(reads)

### load the tables of BP, CC and MF extracted from ClueGO
biological_process   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_BP_06.02.2025.txt"))
cellular_component   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_CC_06.02.2025.txt"))
molecular_functions  <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_MF_06.02.2025.txt"))

### load the modified tables of BP, CC and MF terms
filtered_biological_process_transf  <- readRDS(file.path(GO_data_dir,"filtered_biological_process_transf.RDS"))
filtered_molecular_functions_transf <- readRDS(file.path(GO_data_dir,"filtered_molecular_functions_transf.RDS"))
filtered_cellular_component_transf  <- readRDS(file.path(GO_data_dir,"filtered_cellular_component_transf.RDS"))

######################################################
##  Transform (translate) from NCBI ids to Symbols  ##
######################################################

translate_term=function(term,annotations){# annotations is in entrez id
  return(annotations$SYMBOL[which(annotations$ENTREZID %in% term)])
}

system.time({
  translated_bp=sapply(filtered_biological_process_transf,translate_term,annotations=gene_types_bp)
})
system.time({
  translated_mf=sapply(filtered_molecular_functions_transf,translate_term,annotations=gene_types_mf)
})
system.time({
  translated_cc=sapply(filtered_cellular_component_transf,translate_term,annotations=gene_types_cc)
})

dir.create(file.path(GO_data_dir,"002_Enrichment_lists"),showWarnings = F)

saveRDS(translated_bp,file.path(GO_data_dir,"002_Enrichment_lists","BP_list_transformed.RDS"))
saveRDS(translated_cc, file.path(GO_data_dir,"002_Enrichment_lists", "CC_list_transformed.RDS"))
saveRDS(translated_mf,file.path(GO_data_dir,"002_Enrichment_lists","MF_list_transformed.RDS"))

#First of all, we will extract the gene names that are contained within these 
#lists. With this gene names, we will then create the boolean matrix

biological_process_genes <- unlist(translated_bp)
biological_process_genes <- unique(biological_process_genes)

cellular_component_genes <- unlist(translated_cc)
cellular_component_genes <- unique(cellular_component_genes)

molecular_functions_genes <- unlist(translated_mf)
molecular_functions_genes <- unique(molecular_functions_genes)

names_vector <- unique(c(biological_process_genes, cellular_component_genes,molecular_functions_genes))

length(names_vector)
# 19920
#There are 19001 (19678 old data) genes with annotations either in biological process or 
#cellular component. 

#We will first initialize the matrix

GO_gene_term_matrix_bp <- matrix(data=NA, nrow=length(names_vector), 
                                 ncol=length(translated_bp))

#We name the matrix
rownames(GO_gene_term_matrix_bp) <- names_vector

colnames(GO_gene_term_matrix_bp) <- names(translated_bp)

#And we mark as TRUE if the gene in row i is contained in the term that 
#corresponds to row j

for(i in 1:length(names_vector)) {
  for (j in 1:length(translated_bp)) {
    if (names_vector[i] %in% translated_bp[[j]]) {
      GO_gene_term_matrix_bp[i,j] <- TRUE
    }
    else {
      GO_gene_term_matrix_bp[i,j] <- FALSE
    }
  }
}


#We can repeat it for the cellular component

GO_gene_term_matrix_cc <- matrix(data=NA, nrow=length(names_vector), 
                                 ncol=length(translated_cc))


rownames(GO_gene_term_matrix_cc) <- names_vector

colnames(GO_gene_term_matrix_cc) <- names(translated_cc)


for(i in 1:length(names_vector)) {
  for (j in 1:length(translated_cc)) {
    if (names_vector[i] %in% translated_cc[[j]]) {
      GO_gene_term_matrix_cc[i,j] <- TRUE
    }
    else {
      GO_gene_term_matrix_cc[i,j] <- FALSE
    }
  }
}

#We can repeat it for the molecular functions

GO_gene_term_matrix_mf <- matrix(data=NA, nrow=length(names_vector), 
                                 ncol=length(translated_mf))


rownames(GO_gene_term_matrix_mf) <- names_vector

colnames(GO_gene_term_matrix_mf) <- names(translated_mf)


for(i in 1:length(names_vector)) {
  for (j in 1:length(translated_mf)) {
    if (names_vector[i] %in% translated_mf[[j]]) {
      GO_gene_term_matrix_mf[i,j] <- TRUE
    }
    else {
      GO_gene_term_matrix_mf[i,j] <- FALSE
    }
  }
}


length(unique(rownames(GO_gene_term_matrix_bp)))
length(unique(rownames(GO_gene_term_matrix_cc)))
length(unique(rownames(GO_gene_term_matrix_mf)))

#And now, we can save the boolean matrices
dir.create(file.path(GO_data_dir,"003_Boolean_matrices"),showWarnings = F)

saveRDS(GO_gene_term_matrix_bp, file=file.path(GO_data_dir,"003_Boolean_matrices/boolean_matrix_bp.RDS"))
saveRDS(GO_gene_term_matrix_cc, file=file.path(GO_data_dir,"003_Boolean_matrices/boolean_matrix_cc.RDS"))
saveRDS(GO_gene_term_matrix_mf, file=file.path(GO_data_dir,"003_Boolean_matrices/boolean_matrix_mf.RDS"))





