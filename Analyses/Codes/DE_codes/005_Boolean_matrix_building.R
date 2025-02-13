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

human_ncbi <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","human_ncbi_dataset.txt"))
human_ncbi <- human_ncbi[which(human_ncbi$Gene.Type=="PROTEIN_CODING"),]
human_ncbi <- human_ncbi[order(human_ncbi$NCBI.GeneID),]

# 1.IBA (Inferred from Biological aspect of Ancestor): Inferido a partir de un 
# aspecto biológico del ancestro.
# 2.	IC (Inferred by Curator): Inferido por un curador (experto).
# 3.	IDA (Inferred from Direct Assay): Inferido a partir de un ensayo directo
#  (experimentos realizados en el laboratorio).
# 4.	IEA (Inferred from Electronic Annotation): Inferido a partir de una 
# anotación electrónica (anotación automatizada, como la alineación de secuencias).
# 5.	IMP (Inferred from Mutant Phenotype): Inferido a partir de un fenotipo
#  mutante (observación de un fenotipo en mutantes).
# 6.	ISS (Inferred from Sequence or Structural Similarity): Inferido a partir 
# de similitud en secuencia o estructura.
# 7.	NAS (Non-traceable Author Statement): Declaración de autor no rastreable 
# (anotaciones que se basan en artículos pero no tienen una evidencia 
# experimental directa).
# 8.	TAS (Traceable Author Statement): Declaración de autor rastreable (basado
#  en estudios de autor pero con referencia precisa a la fuente).

GO_2_list <- function(GO_table) { #We take as input any table of GO term, as is downloaded by ClueGO
  
  fgsea_list <- list()
  names <- c()
  for (i in 1:length(GO_table[,1])){
    
    #We first substitute the whitespaces, to put the names as in the fgsea named lists format
    name_string <- GO_table[i,3]
    names <- c(names,name_string)              #Here we concatenate the names
    
    #And now we split and extract our genes of interest for the corresponding term
    genes_vec <- GO_table[i,4]
    genes_vec <- str_split(genes_vec,"\\|")
    genes_vec <- genes_vec[[1]]
    genes_vec <- gsub(".*:(\\d+)", "\\1", genes_vec)  #We use this expression because all gene IDs are preceded by several letters in the clueGO terms
    fgsea_list[[i]] <- unique(genes_vec)                    #Some of them are repeated, so we have to make them unique
  }
  names(fgsea_list) <- names                                #We assign the name of each term and return the list
  return(fgsea_list)
}

# Here I separate genes from all.genes.in.the.node column. They have entrez gene id code
biological_process_transf   <- GO_2_list(biological_process)
cellular_component_transf   <- GO_2_list(cellular_component)
molecular_functions_transf  <- GO_2_list(molecular_functions)

# https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/
# This is the good one. Contain all the genes I think
human_gene_info <- read.delim(file.path(GO_data_dir,"Homo_sapiens.gene_info"))

# Lets change the GeneID column by NCBI_ID 
colnames(human_gene_info)[2] <- "NCBI.GeneID"
human_gene_info <- human_gene_info[,c("NCBI.GeneID","Symbol","Synonyms","dbXrefs",
                                      "description","type_of_gene")]
# We keep just the protein coding genes
human_gene_info <- human_gene_info[which(human_gene_info$type_of_gene=="protein-coding"),]
# We separate the dbXref based on the pattern found. They are various ID identifiers 
splitted_values <- strsplit(human_gene_info$dbXrefs,"\\|")
ordered_elements <- function(fila) {
  mim <- grep("^MIM:", fila, value = TRUE)
  hgnc <- grep("^HGNC:", fila, value = TRUE)
  ensembl <- grep("^Ensembl:", fila, value = TRUE)
  alliance <- grep("^AllianceGenome:", fila, value = TRUE)
  c(
    ifelse(length(mim) > 0, mim, NA),
    ifelse(length(hgnc) > 0, hgnc, NA),
    ifelse(length(ensembl) > 0, ensembl, NA),
    ifelse(length(alliance) > 0, alliance, NA)
  )
}

df_split <- as.data.frame(do.call(rbind, lapply(splitted_values, ordered_elements)), stringsAsFactors = FALSE)
colnames(df_split) <- c("MIM","HGNC","Ensembl","AllianceGenome")
remove_before_two_points <- function(x) {
  sub(".*:(.*)$", "\\1", x)
}

# Apply the function to each column
df_split$MIM <- remove_before_two_points(df_split$MIM)
df_split$HGNC <- remove_before_two_points(df_split$HGNC)
df_split$Ensembl <- remove_before_two_points(df_split$Ensembl)
df_split$AllianceGenome <- remove_before_two_points(df_split$AllianceGenome)

# bind to the main dataset
human_gene_info <- cbind(human_gene_info,df_split)

df_duplicated_symbols <- data.frame(human_gene_info[duplicated(human_gene_info$Symbol), ])

write.xlsx(df_duplicated_symbols,file.path(GO_data_dir,"001_GO_table_to_list","duplicated_symbols_entrez_id.xlsx"))

human_gene_info <- human_gene_info[order(human_gene_info$NCBI.GeneID),]

human_gene_info <- human_gene_info %>% filter(!Symbol %in% df_duplicated_symbols$Symbol)

length(unique(human_gene_info$Symbol))
# 20582
length(unique(human_ncbi$Symbol))
# 20594
### Check which elements are not in human_gene_info

genes_out <- human_ncbi[which(human_ncbi$Symbol %in% human_gene_info$Symbol),]
genes_not_in_human_gene_info <- data.frame(genes_not_in_human_gene_info=setdiff(human_ncbi$Symbol, human_gene_info$Symbol))
print(genes_not_in_human_gene_info)
genes_not_in_human_gene_info
# 1                       MT-ATP6
# 2                       MT-ATP8
# 3                        MT-CO1
# 4                        MT-CO2
# 5                        MT-CO3
# 6                        MT-CYB
# 7                        MT-ND1
# 8                        MT-ND2
# 9                        MT-ND3
# 10                       MT-ND4
# 11                      MT-ND4L
# 12                       MT-ND5
# 13                       MT-ND6
write.xlsx(genes_not_in_human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","genes_not_in_human_gene_info.xlsx"))
write.xlsx(human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info_entrez_id.xlsx"))
write.xlsx(human_ncbi,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_ncbi.xlsx"))

write.table(human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info_entrez_id.txt"))
write.table(human_ncbi,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_ncbi.txt"))

trans_gene_id <- human_ncbi
# trans_gene_id <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info.txt"))

#We can check how many of our gene IDs are contained in the file

# Convert the list elements to characters  
genes_in_list_bp <- unique(unlist(biological_process_transf))
genes_in_list_bp <- as.character(genes_in_list_bp)

genes_in_list_mf <- unique(unlist(molecular_functions_transf))
genes_in_list_mf <- as.character(genes_in_list_mf)

genes_in_list_cc <- unique(unlist(cellular_component_transf))
genes_in_list_cc <- as.character(genes_in_list_cc)

# Covert the NCBI IDs to characters. At the end we stayed with the NCBI dataset
# and not with the human_gene_info 
target_genes <- as.character(trans_gene_id$NCBI.GeneID)

##############################################################
##  Filter alll the elements that match with trans_gene_id  ##
##############################################################
### For biological processes
length(genes_in_list_bp)
sum(genes_in_list_bp %in% target_genes)

length(which(!genes_in_list_bp %in% target_genes))
non_common_genes_bp <- genes_in_list_bp[which(!genes_in_list_bp %in% target_genes)]

# initial_sizes=sapply(biological_process_transf,length)
filter_set=function(set,ref){
  return(set[which(set %in% ref)])
}

filtered_biological_process_transf=biological_process_transf

start=Sys.time()
filtered_biological_process_transf=sapply(biological_process_transf,filter_set,ref=target_genes)
end=Sys.time()
end-start

length(unique(unlist(filtered_biological_process_transf)))

genes_in_list_bp_filtered <- unique(unlist(filtered_biological_process_transf))
genes_in_list_bp_filtered <- as.character(genes_in_list_bp_filtered)
sum(genes_in_list_bp_filtered %in% target_genes)

### For molecular functions
length(unique(unlist(molecular_functions_transf)))
sum(genes_in_list_mf %in% target_genes)
filtered_molecular_functions_transf=molecular_functions_transf

start=Sys.time()
filtered_molecular_functions_transf=sapply(molecular_functions_transf,filter_set,ref=target_genes)
end=Sys.time()
end-start

length(unique(unlist(filtered_molecular_functions_transf)))

genes_in_list_mf_filtered <- unique(unlist(filtered_molecular_functions_transf))
genes_in_list_mf_filtered <- as.character(genes_in_list_mf_filtered)
sum(genes_in_list_mf_filtered %in% target_genes)

### For Cellular Components
length(unique(unlist(cellular_component_transf)))
sum(genes_in_list_cc %in% target_genes)
filtered_cellular_component_transf=cellular_component_transf

start=Sys.time()
filtered_cellular_component_transf=sapply(cellular_component_transf,filter_set,ref=target_genes)
end=Sys.time()
end-start

length(unique(unlist(filtered_cellular_component_transf)))

genes_in_list_cc_filtered <- unique(unlist(filtered_cellular_component_transf))
genes_in_list_cc_filtered <- as.character(genes_in_list_cc_filtered)
sum(genes_in_list_cc_filtered %in% target_genes)

######################################################
##  Transform (translate) from NCBI ids to Symbols  ##
######################################################

translate_term=function(term,annotations){
  return(annotations$Symbol[which(annotations$NCBI.GeneID %in% term)])
}

system.time({
  translated_bp=sapply(filtered_biological_process_transf,translate_term,annotations=trans_gene_id)
})
system.time({
  translated_mf=sapply(filtered_molecular_functions_transf,translate_term,annotations=trans_gene_id)
})
system.time({
  translated_cc=sapply(filtered_cellular_component_transf,translate_term,annotations=trans_gene_id)
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





