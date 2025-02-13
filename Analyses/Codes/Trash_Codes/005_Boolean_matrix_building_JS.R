
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
main_wd     <- getwd()
input_dir   <- file.path("Analyses","Inputs")
output_dir  <- file.path("Analyses","Outputs")
GO_data_dir <- file.path(input_dir,"004_GO_data")
#Date of GO tables version: 11.02.2020

dir.create(file.path(GO_data_dir,"GO_table_to_list"),showWarnings = F)

### Load any table with the genes names
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
genes_names <- rownames(reads)
human_ncbi <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","human_ncbi_dataset.tsv"))
types=data.frame(summary(factor(human_ncbi$Gene.Type)))
biotypes=rownames(types)
biotypes=biotypes[c(3,5,6,7,8,9,10)]
human_ncbi=human_ncbi[which(human_ncbi$Gene.Type %in% biotypes),]
human_ncbi           <- human_ncbi[,c("NCBI.GeneID","Symbol","Gene.Type","Ensembl.GeneIDs")] 

### load the tables of BP, CC and MF extracted from ClueGO
biological_process   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo_Sapiens_GO_BP_22.11.2024.txt"))
cellular_component   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo_Sapiens_GO_CC_22.11.2024.txt"))
molecular_functions  <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo_Sapiens_GO_MF_22.11.2024.txt"))

########################
##  Define Functions  ##
########################

#The following function converts a table of GO enrichments into a named list with 
#the term names as names, and as members of each list item, 
#the Gene IDs that are included in the table.

#The reason why I produced a list it is because I found it more manageable, 
#and it can be extracted for other uses. Originally, this was the
#function that I used to be able to perform a GSEA analysis on Mtb GO terms 
#(for homo sapiens)

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

biological_process_transf   <- GO_2_list(biological_process)
cellular_component_transf   <- GO_2_list(cellular_component)
molecular_functions_transf  <- GO_2_list(molecular_functions)

# We can translate those identifiers to whichever we may find more suitable. 
# Here, for example, we will load the tables for humans
#biological processes, cellular components and molecular functions, 
#create the gsea list, and then 
#translate the terms to the locus tag identifiers (RV numbers):

# Now we are going to translate the gene symbols to entrezgene_id. We can do it 
# only to our data or to the full gene name list downloaded from ncbi in tsv format

# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# gene_names <- human_ncbi$Symbol
# # gene_names <- rownames(reads)
# 
# # Use getBM to map the gene names to another attributes
# results <- getBM(
#   attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene_id"), # attributes you wish too extract
#   filters = "external_gene_name",  # ID of your input
#   values = gene_names,             # Your gene names
#   mart = ensembl                   # connection to dataset 
# )
# 
# # Show results
# head(results,5)

#write.table(human_ncbi,file.path(GO_data_dir,"001_GO_table_to_list","sym_to_entrezgene_id.txt"))

#trans_gene_id <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","sym_to_entrezgene_id.txt"))

#We can check how many of our gene IDs are contained in the file

# Convert the list elements to the 
genes_in_list_bp <- unique(unlist(biological_process_transf))
genes_in_list_bp <- as.character(genes_in_list_bp)

genes_in_list_mf <- unique(unlist(molecular_functions_transf))
genes_in_list_mf <- as.character(genes_in_list_mf)

genes_in_list_cc <- unique(unlist(cellular_component_transf))
genes_in_list_cc <- as.character(genes_in_list_cc)

# Covert the NCBI IDs to characters 
human_ncbi$NCBI.GeneID <- as.character(human_ncbi$NCBI.GeneID)


###########################################################
##  Filter all the elements that match with NCBI.GeneID  ##
###########################################################

### For biological processes
length(genes_in_list_bp)
sum(genes_in_list_bp %in% human_ncbi$NCBI.GeneID)
bp_vec <- seq_along(biological_process_transf)

initial_sizes=sapply(biological_process_transf,length)
filtered_biological_process_transf=biological_process_transf

start=Sys.time()
for (i in bp_vec) {
  # Filter all the elements that match with trans_gene_id
  filtered_biological_process_transf[[i]] <- biological_process_transf[[i]][which(
    biological_process_transf[[i]] %in% human_ncbi$NCBI.GeneID)]
}
end=Sys.time()
end-start

filter_set=function(set,ref){
  return(set[which(set %in% ref)])
}

start=Sys.time()
filtered_biological_process_transf=sapply(biological_process_transf,filter_set,ref=human_ncbi$NCBI.GeneID)
end=Sys.time()
end-start

genes_in_list_bp_filtered <- unique(unlist(filtered_biological_process_transf))
genes_in_list_bp_filtered <- as.character(genes_in_list_bp_filtered)

### For molecular functions
length(unique(unlist(molecular_functions_transf)))
sum(genes_in_list_mf %in% target_genes)
mf_vec <- seq_along(molecular_functions_transf)

for (i in mf_vec) {
  # Filter alll the elements that match with trans_gene_id
  molecular_functions_transf[[i]] <- molecular_functions_transf[[i]][
    molecular_functions_transf[[i]] %in% trans_gene_id$NCBI.GeneID
  ]
}
length(unique(unlist(molecular_functions_transf)))
sum(genes_in_list_mf %in% target_genes)

### For Cellular Components
length(unique(unlist(cellular_component_transf)))
sum(genes_in_list_cc %in% target_genes)
cc_vec <- seq_along(cellular_component_transf)

for (i in cc_vec) {
  # Filter alll the elements that match with trans_gene_id
  cellular_component_transf[[i]] <- cellular_component_transf[[i]][
    cellular_component_transf[[i]] %in% trans_gene_id$NCBI.GeneID
  ]
}


######################################################
##  Transform (translate) from NCBI ids to Symbols  ##
######################################################

# Make sure we defined well enough the bp_vec
bp_vec <- seq_along(filtered_biological_process_transf)

#  Introduce a progress bar
pb <- txtProgressBar(min = 0, max = length(bp_vec), style = 3) # progress bar

translate_term=function(term,annotations){
  return(annotations$Symbol[which(annotations$NCBI.GeneID %in% term)])
}

translated_bp=sapply(biological_process_transf,translate_term,annotations=human_ncbi)


system.time({
  
  for (i in bp_vec) { # terminos GO
    setTxtProgressBar(pb, i) # Update progress bar
    for (j in 1:length(biological_process_transf[[i]])) { ## genes dentro de cada término
      
      match_idx <- which(trans_gene_id$NCBI.GeneID == biological_process_transf[[i]][j])
      biological_process_transf[[i]][j] <- trans_gene_id$Symbol[match_idx]
      
      # if (length(match_idx) > 0) {
      #   
      #   biological_process_transf[[i]][j] <- trans_gene_id$Symbol[match_idx]
      #   
      # }
      # # else{
      #   biological_process_transf[[i]][j] <- biological_process_transf[[i]][j]
      # }
      
    }
  }
  resultado <- sum(1:1e7)
})
close(bp)

dir.create(file.path(GO_data_dir,"002_Enrichment_lists"),showWarnings = F)
saveRDS(biological_process_transf,file.path(GO_data_dir,"002_Enrichment_lists","BP_list_transformed.RDS"))

# Make sure we defined well enough the cc_vec
cc_vec <- seq_along(cellular_component)

# Function to process an CC element
process_cc <- function(i, cellular_component_transf, trans_gene_id) {
  for (j in 1:length(cellular_component_transf[[i]])) {
    match_idx <- which(trans_gene_id$NCBI.GeneID == cellular_component_transf[[i]][j])
    if (length(match_idx) > 0) {
      cellular_component_transf[[i]][j] <- trans_gene_id$Symbol[match_idx]
    }
  }
  return(cellular_component_transf[[i]]) # Return the processed element
}

num_cores <- detectCores() - 1

# Execute in parallel
system.time({
  transformed_cc <- mclapply(
    cc_vec, 
    process_cc, 
    cellular_component_transf = cellular_component_transf, 
    trans_gene_id = trans_gene_id, 
    mc.cores = num_cores
  )
})

### Reconstruct the original list
names(transformed_cc) <- names(cellular_component_transf)
cellular_component_transf <- transformed_cc

# Save the results
saveRDS(cellular_component_transf, file.path(GO_data_dir,"002_Enrichment_lists", "CC_list_transformed.RDS"))

# Make sure we defined well enough the mf_vec
mf_vec <- seq_along(molecular_functions)

#  Introduce a progress bar
pb <- txtProgressBar(min = 0, max = length(bp_vec), style = 3) # progress bar

system.time({
  
  for (i in mf_vec) {
    setTxtProgressBar(pb, i) # Update progress bar
    for (j in 1:length(molecular_functions_transf[[i]])) {
      
      match_idx <- which(trans_gene_id$NCBI.GeneID == molecular_functions_transf[[i]][j])
      molecular_functions_transf[[i]][j] <- trans_gene_id$Symbol[match_idx]
      
      # if (length(match_idx) > 0) {
      #   
      #   molecular_functions_transf[[i]][j] <- trans_gene_id$Symbol[match_idx]
      #   
      # }
      # # else{
      #   molecular_functions_transf[[i]][j] <- molecular_functions_transf[[i]][j]
      # }
      
    }
  }
  resultado <- sum(1:1e7)
})
close(bp)

saveRDS(molecular_functions_transf,file.path(GO_data_dir,"002_Enrichment_lists","MF_list_transformed.RDS"))

#First of all, we will extract the gene names that are contained within these 
#lists. With this gene names, we will then create the boolean matrix

biological_process_genes <- unlist(biological_process_transf)
biological_process_genes <- unique(biological_process_genes)

cellular_component_genes <- unlist(cellular_component_transf)
cellular_component_genes <- unique(cellular_component_genes)

names_vector <- unique(c(biological_process_genes, cellular_component_genes))

length(names_vector)

#There are 19678 genes with annotations either in biological process or 
#cellular component. 

#We will first initialize the matrix

GO_gene_term_matrix_bp <- matrix(data=NA, nrow=length(names_vector), 
                                 ncol=length(biological_process_transf))

#We name the matrix

rownames(GO_gene_term_matrix_bp) <- names_vector

colnames(GO_gene_term_matrix_bp) <- names(biological_process_transf)

#And we mark as TRUE if the gene in row i is contained in the term that 
#corresponds to row j

for(i in 1:length(names_vector)) {
  for (j in 1:length(biological_process_transf)) {
    if (names_vector[i] %in% biological_process_transf[[j]]) {
      GO_gene_term_matrix_bp[i,j] <- TRUE
    }
    else {
      GO_gene_term_matrix_bp[i,j] <- FALSE
    }
  }
}


#We can repeat it for the cellular component

GO_gene_term_matrix_cc <- matrix(data=NA, nrow=length(names_vector), 
                                 ncol=length(cellular_component_transf))


rownames(GO_gene_term_matrix_cc) <- names_vector

colnames(GO_gene_term_matrix_cc) <- names(cellular_component_transf)


for(i in 1:length(names_vector)) {
  for (j in 1:length(cellular_component_transf)) {
    if (names_vector[i] %in% cellular_component_transf[[j]]) {
      GO_gene_term_matrix_cc[i,j] <- TRUE
    }
    else {
      GO_gene_term_matrix_cc[i,j] <- FALSE
    }
  }
}


length(unique(rownames(GO_gene_term_matrix_bp)))
length(unique(rownames(GO_gene_term_matrix_cc)))

#And now, we can save the boolean matrices
dir.create(file.path(GO_data_dir,"003_Boolean_matrices"),showWarnings = F)

saveRDS(GO_gene_term_matrix_bp, file=file.path(GO_data_dir,"003_Boolean_matrices/boolean_matrix_bp.RDS"))
saveRDS(GO_gene_term_matrix_cc, file=file.path(GO_data_dir,"003_Boolean_matrices/boolean_matrix_cc.RDS"))
saveRDS(GO_gene_term_matrix_mf, file=file.path(GO_data_dir,"003_Boolean_matrices/boolean_matrix_mf.RDS"))
