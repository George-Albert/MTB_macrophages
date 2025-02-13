
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

dir.create(file.path(GO_data_dir,"001_GO_table_to_list"),showWarnings = F)

### Load any table with the genes names
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
genes_names <- rownames(reads)

### load the tables of BP, CC and MF extracted from ClueGO
biological_process   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_BP_06.02.2025.txt"))
cellular_component   <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_CC_06.02.2025.txt"))
molecular_functions  <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens_GO_MF_06.02.2025.txt"))

### lets inspect each metada of the genes symbols we already found
# NCBI data 
human_ncbi <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","human_ncbi_dataset.txt"))
gene2ensembl <- read.delim(file.path(GO_data_dir,"001_GO_table_to_list","Homo Sapiens.gene2ensembl_2025.02.06.txt"))

types=data.frame(Summary=summary(factor(human_ncbi$Gene.Type)))
biotypes=rownames(types)
names(biotypes) <- biotypes
biotypes=biotypes[c("ncRNA","PROTEIN_CODING","PSEUDO","rRNA","scRNA","snoRNA","snRNA")]
# biotypes=biotypes[c("PROTEIN_CODING")]
human_ncbi=human_ncbi[which(human_ncbi$Gene.Type %in% biotypes),]
human_ncbi <- human_ncbi[,c("NCBI.GeneID","Symbol","Gene.Type","Ensembl.GeneIDs")] 

### gpi file from GO database from human file
### https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/

# !gpi-version: 1.2
# !
# !The set of protein accessions included in this file is based on UniProt 
# reference proteomes, which provide one protein per gene.
# !They include the protein sequences annotated in Swiss-Prot or the longest 
# TrEMBL transcript if there is no Swiss-Prot record.
# !Protein accessions are represented in this file even if there is no 
# associated GO annotation.
# !
human_gpi <- read.delim(file.path(GO_data_dir,"goa_human.gpi"),skip = 22,
                        header = F)
col_names_vec <- c("DB","DB_Object_ID","DB_Object_Symbol","DB_Object_Name",
                   "DB_Object_Synonym","DB_Object_Type","Taxon","Parent_Object_ID",
                   "DB_Xref","Properties")
colnames(human_gpi) <- col_names_vec

###  Human .goa file from GO data base
###  https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/
###  This is the file ClueGO uses to update the GO terms 
human_GO <- read.delim(file.path(GO_data_dir,"25.H_sapiens.goa"),skip=4,
                       header = F)
col_names_vec_GO <- c("DB","DB_Object_ID","DB_Object_Symbol","Qualifier",
                      "GO ID","DB_Reference","Evidence Code","With_or_From",
                      "Aspect","DB_Object_Name","DB_Object_Synonym","DB_Object_Type",
                      "Taxon","Date","Assigned_By","Annotation_Extension",
                      "Gene_Product_Form_ID")
colnames(human_GO) <- col_names_vec_GO

length(unique(human_GO$DB_Object_Symbol))
human_GO <- human_GO[!duplicated(human_GO$DB_Object_Symbol),]

# https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/
#This is the good one. Contain all the genes i think
human_gene_info <- read.delim(file.path(GO_data_dir,"Homo_sapiens.gene_info"))
# Lets change the GeneID column by NCBI_ID 
colnames(human_gene_info)[2] <- "NCBI.GeneID"
human_gene_info <- human_gene_info[,c("NCBI.GeneID","Symbol","Synonyms","dbXrefs",
                                      "description","type_of_gene")]
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

# Aplicar la función a cada columna
df_split$MIM <- remove_before_two_points(df_split$MIM)
# df_split$HGNC<- replace(df_split$HGNC, is.na(df_split$HGNC), "HGNC:HGNC:")
df_split$HGNC <- remove_before_two_points(df_split$HGNC)
df_split$Ensembl <- remove_before_two_points(df_split$Ensembl)
df_split$AllianceGenome <- remove_before_two_points(df_split$AllianceGenome)

human_gene_info <- cbind(human_gene_info,df_split)
### load the tables of BP, CC and MF extracted from ClueGO

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
# biological processes, cellular components and molecular functions, 
# create the gsea list, and then 
# translate the terms to the Symbol identifiers:

# Now we are going to translate the gene symbols to entrezgene_id. We can do it 
# only to our data or to the full gene name list downloaded from ncbi in tsv format

# Mymarts <- listMarts()
# host="https://may2024.archive.ensembl.org"
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# datasets <- listDatasets(ensembl)
# list.filters <- listFilters(ensembl)
# gene_names <- human_ncbi$Symbol
# gene_names <- non_common_genes_bp
# # # gene_names <- rownames(reads)
# attributes_mart <- listAttributes(ensembl)
# # Use getBM to map the gene names to another attributes
# results <- getBM(
# attributes = c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),
# filters = "with_entrezgene",
# values = TRUE,  # 'TRUE' indica que queremos las entradas con NCBI Gene IDs
# mart = ensembl
# )
# 
# # Show results
# head(results,5)

write.table(human_ncbi,file.path(GO_data_dir,"001_GO_table_to_list","human_ncbi_dataset.txt"))
write.table(human_GO,file.path(GO_data_dir,"001_GO_table_to_list","25.H_sapiens.txt"))
write.table(human_gpi,file.path(GO_data_dir,"001_GO_table_to_list","goa_human_gpi.txt"))
write.table(human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info.txt"))

trans_gene_id <- human_gene_info
# trans_gene_id <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info.txt"))

#We can check how many of our gene IDs are contained in the file

# Convert the list elements to characters  
genes_in_list_bp <- unique(unlist(biological_process_transf))
genes_in_list_bp <- as.character(genes_in_list_bp)

genes_in_list_mf <- unique(unlist(molecular_functions_transf))
genes_in_list_mf <- as.character(genes_in_list_mf)

genes_in_list_cc <- unique(unlist(cellular_component_transf))
genes_in_list_cc <- as.character(genes_in_list_cc)

# Covert the NCBI IDs to characters 
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

#There are 19678 genes with annotations either in biological process or 
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
