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

### This is from the NCBI web
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

# gene.info.url = ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
# https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/
# We are going to inspect clueGO gene2 ensembl file too after this file form NCBI
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
gene2ensemble <- read.delim(file.path(GO_data_dir,"Homo Sapiens.gene2ensembl_2022.05.25.txt"))
gene2accession <- read.delim(file.path(GO_data_dir,"Homo Sapiens.gene2accession_2022.05.25.txt"))

# Downloaded from 
# https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
hgnc_complete_set <- read.delim(file.path(GO_data_dir,"hgnc_complete_set.txt"))
hgnc_subset <-hgnc_complete_set[,c("hgnc_id","symbol","entrez_id","ensembl_gene_id","locus_group")] 

### lets see the gene types of this list
gene_type_vec <- data.frame(unique(gene2ensemble$GeneType))

### lets keep all the gene types that contain protein coding in their description
# gene2ensemble <- gene2ensemble[grepl("protein_coding",gene2ensemble$GeneType),]
# gene2accession <- gene2accession[grepl("protein-coding",gene2accession$GeneType),]

### lets re-check the gene types again
gene_type_vec2 <- data.frame(unique(gene2ensemble$GeneType))

write.xlsx(genes_not_in_human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","genes_not_in_human_gene_info.xlsx"))
write.xlsx(human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info_entrez_id.xlsx"))
write.xlsx(human_ncbi,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_ncbi.xlsx"))

write.table(human_gene_info,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info_entrez_id.txt"))
write.table(human_ncbi,file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_ncbi.txt"))

####We are going to use the gene2ensemble or gene2accession files from clueGO 
####to assure the match in the annotation scheme
# trans_gene_id <- human_ncbi
trans_gene_id <- hgnc_subset
### Change the name of the UniqueID.EntrezGeneID (entrez_id) column to NCBI.GeneID
colnames(trans_gene_id)[3] <- "NCBI.GeneID"
# trans_gene_id <- read.table(file.path(GO_data_dir,"001_GO_table_to_list","Homo_sapiens_gene_info.txt"))

library(org.Hs.eg.db)

# Lista de tus 290 genes (Ejemplo)
gene_list <- as.character(non_common_genes_bp$non_common_genes_bp)  # Sustituye por tus Entrez IDs
gene_list <- genes_in_list_bp
# Obtener información sobre el tipo de gen
gene_types <- select(org.Hs.eg.db, keys = gene_list,
                     columns = c("ENTREZID","SYMBOL","GENETYPE"), 
                     keytype = "ENTREZID")

length(unique(gene_types$SYMBOL))

genes_duplicados <- gene_types[duplicated(gene_types$SYMBOL),]

# Mostrar resultados
print(pseudogenes)


library(biomaRt)

# Conectar con Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Obtener información sobre el tipo de gen
genes_info <- getBM(attributes = c('entrezgene_id', 'ensembl_gene_id', 'hgnc_symbol',
                                   'external_gene_name', 'gene_biotype'),
                    filters = 'entrezgene_id',
                    values = gene_list,
                    mart = ensembl)

length(unique(genes_info$hgnc_symbol))
length(unique(genes_info$entrezgene_id))
length(unique(genes_info$ensembl_gene_id))

genes_duplicados_1 <- genes_info[duplicated(genes_info$entrezgene_id),]


