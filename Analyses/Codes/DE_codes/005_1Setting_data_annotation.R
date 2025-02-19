{
  library(tidyverse)
  library(stringr)
  library(biomaRt)
  library(parallel)
  library(dplyr)
  library(babelgene)
  library(HGNChelper)
  library(xlsx)
  library(org.Hs.eg.db)
}


#################
##  Functions  ##
#################

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

#################
##  Load data  ##
#################
main_wd     <- getwd()
input_dir   <- file.path("Analyses","Inputs")
output_dir  <- file.path("Analyses","Outputs")
GO_data_dir <- file.path(input_dir,"004_GO_data")

###Date of GO tables version: 06.02.2025
dir.create(file.path(GO_data_dir,"001_GO_table_to_list"),showWarnings = F)

### Load our read table with the genes names
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
genes_names <- rownames(reads)

### load the tables of BP, CC and MF extracted from ClueGO
biological_process   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_BP_06.02.2025.txt"))
cellular_component   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_CC_06.02.2025.txt"))
molecular_functions  <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_MF_06.02.2025.txt"))

# Here I separate genes from all.genes.in.the.node column. They have entrez gene id code
biological_process_transf   <- GO_2_list(biological_process)
cellular_component_transf   <- GO_2_list(cellular_component)
molecular_functions_transf  <- GO_2_list(molecular_functions)

#We can check how many of our gene IDs are contained in the file

# Convert the list elements to characters  
genes_in_list_bp <- unique(unlist(biological_process_transf))
genes_in_list_bp <- as.character(genes_in_list_bp)

genes_in_list_mf <- unique(unlist(molecular_functions_transf))
genes_in_list_mf <- as.character(genes_in_list_mf)

genes_in_list_cc <- unique(unlist(cellular_component_transf))
genes_in_list_cc <- as.character(genes_in_list_cc)

# Lista de tus 290 genes (Ejemplo)
gene_list <- genes_in_list_bp
gene_list_cc <- genes_in_list_cc
gene_list_mf <- genes_in_list_mf

# Obtener información sobre el tipo de gen
gene_types_bp <- select(org.Hs.eg.db, keys = gene_list,
                     columns = c("ENTREZID","SYMBOL","GENETYPE"), 
                     keytype = "ENTREZID")
gene_types_cc <- select(org.Hs.eg.db, keys = gene_list_cc,
                        columns = c("ENTREZID","SYMBOL","GENETYPE"), 
                        keytype = "ENTREZID")
gene_types_mf <- select(org.Hs.eg.db, keys = gene_list_mf,
                        columns = c("ENTREZID","SYMBOL","GENETYPE"), 
                        keytype = "ENTREZID")

length(unique(gene_types$SYMBOL))
length(unique(gene_types_cc$SYMBOL))
length(unique(gene_types_mf$SYMBOL))

class(gene_types_bp$ENTREZID)
gene_types_bp <-gene_types_bp[order(as.integer(gene_types_bp$ENTREZID)),]
gene_types_cc <-gene_types_cc[order(as.integer(gene_types_cc$ENTREZID)),]
gene_types_mf <-gene_types_mf[order(as.integer(gene_types_mf$ENTREZID)),]

genes_duplicados    <- gene_types_bp[duplicated(gene_types_bp$SYMBOL),]

# ENTREZID SYMBOL GENETYPE
# 8486  100124696    TEC  unknown
# 12230 100187828    HBD  unknown
# 9249  100505381   MMD2  unknown
# 9522  133834869   <NA>     <NA>
# 15495 139281660   <NA>     <NA>

genes_duplicados_cc <- gene_types_cc[duplicated(gene_types_cc$SYMBOL),]

# ENTREZID SYMBOL GENETYPE
# 8611  100124696    TEC  unknown
# 5242  100187828    HBD  unknown
# 353   100505381   MMD2  unknown
# 16054 133834869   <NA>     <NA>
# 5588  139281660   <NA>     <NA>

genes_duplicados_mf <- gene_types_mf[duplicated(gene_types_mf$SYMBOL),]

# ENTREZID SYMBOL GENETYPE
# 1157  100124696    TEC  unknown
# 6997  100187828    HBD  unknown
# 13425 100505381   MMD2  unknown
# 13540 133834869   <NA>     <NA>
# 843   139281660   <NA>     <NA>

##############################################################
##  Filter all the elements that match with trans_gene_id  ##
##############################################################

### For biological processes
class(genes_in_list_bp)
target_genes_bp <- gene_types_bp$ENTREZID

length(genes_in_list_bp)
sum(genes_in_list_bp %in% target_genes_bp)
length(which(!genes_in_list_bp %in% target_genes_bp))

filter_set=function(set,ref){
  return(set[which(set %in% ref)])
}

filtered_biological_process_transf=biological_process_transf

start=Sys.time()
filtered_biological_process_transf=sapply(biological_process_transf,filter_set,ref=target_genes_bp)
end=Sys.time()
end-start

length(unique(unlist(filtered_biological_process_transf)))

genes_in_list_bp_filtered <- unique(unlist(filtered_biological_process_transf))
genes_in_list_bp_filtered <- as.character(genes_in_list_bp_filtered)
sum(genes_in_list_bp_filtered %in% target_genes_bp)


### For molecular functions
class(genes_in_list_mf)
target_genes_mf <- gene_types_mf$ENTREZID

length(unique(unlist(molecular_functions_transf)))
sum(genes_in_list_mf %in% target_genes_mf)

filtered_molecular_functions_transf=molecular_functions_transf

start=Sys.time()
filtered_molecular_functions_transf=sapply(molecular_functions_transf,filter_set,ref=target_genes_mf)
end=Sys.time()
end-start

length(unique(unlist(filtered_molecular_functions_transf)))

genes_in_list_mf_filtered <- unique(unlist(filtered_molecular_functions_transf))
genes_in_list_mf_filtered <- as.character(genes_in_list_mf_filtered)
sum(genes_in_list_mf_filtered %in% target_genes_mf)

### For Cellular Components
class(genes_in_list_cc)
target_genes_cc <- gene_types_cc$ENTREZID

length(unique(unlist(cellular_component_transf)))
sum(genes_in_list_cc %in% target_genes_cc)
filtered_cellular_component_transf=cellular_component_transf

start=Sys.time()
filtered_cellular_component_transf=sapply(cellular_component_transf,filter_set,ref=target_genes_cc)
end=Sys.time()
end-start

length(unique(unlist(filtered_cellular_component_transf)))

genes_in_list_cc_filtered <- unique(unlist(filtered_cellular_component_transf))
genes_in_list_cc_filtered <- as.character(genes_in_list_cc_filtered)
sum(genes_in_list_cc_filtered %in% target_genes_cc)

saveRDS(filtered_biological_process_transf,file.path(GO_data_dir,"filtered_biological_process_transf.RDS"))
saveRDS(filtered_molecular_functions_transf,file.path(GO_data_dir,"filtered_molecular_functions_transf.RDS"))
saveRDS(filtered_cellular_component_transf,file.path(GO_data_dir,"filtered_cellular_component_transf.RDS"))



