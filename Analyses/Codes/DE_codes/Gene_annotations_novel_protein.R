
# Run first 008_Enrichment_Analysis_by_cluster.R 
# to obtain genes.names vector of Ensemble IDs


library(biomaRt)

# Conectarse a la base de datos Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Vector de Ensembl IDs
ensembl_ids <- c("ENSG00000124593", "ENSG00000188897", "ENSG00000241489", 
                 "ENSG00000243696", "ENSG00000244255", "ENSG00000250138", 
                 "ENSG00000250264", "ENSG00000250644", "ENSG00000254979", 
                 "ENSG00000255508", "ENSG00000255730", "ENSG00000256500", 
                 "ENSG00000256514", "ENSG00000257524", "ENSG00000257767", 
                 "ENSG00000258461", "ENSG00000258472", "ENSG00000259132", 
                 "ENSG00000260537", "ENSG00000261884", "ENSG00000267228", 
                 "ENSG00000267303", "ENSG00000267645", "ENSG00000267740", 
                 "ENSG00000268400", "ENSG00000269242", "ENSG00000269711", 
                 "ENSG00000270149", "ENSG00000270299", "ENSG00000271254", 
                 "ENSG00000271741", "ENSG00000272410", "ENSG00000272617", 
                 "ENSG00000273088", "ENSG00000278817", "ENSG00000280987", 
                 "ENSG00000282034", "ENSG00000282988", "ENSG00000283189")

# Obtener los nombres de genes (Gene Symbols) correspondientes
gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),  # Cambiamos los atributos
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Mostrar el resultado
print(gene_annotations)
gene_annotations <- data.frame(gene_annotations)
write.table(gene_annotations,file.path(input_dir,"002_Processed","genotype_matrix_building_files","novel_proteins_wo_gene_name.txt"))

list_marts     <- listMarts()
ensembl        <- useMart("ENSEMBL_MART_ENSEMBL")
list_dataset   <- listDatasets(ensembl)
dataset_to_use <- "hsapiens_gene_ensembl"

# ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset = dataset_to_use)

ensembl <- useDataset(dataset =  dataset_to_use, mart = ensembl)
gene_symbols <- rownames(reads)  # Cambia esto por tus propios genes

# Obtener los Ensembl IDs correspondientes
conversion <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                    filters = 'hgnc_symbol',
                    values = gene_symbols,
                    mart = ensembl)

# Ver los resultados de la conversión
print(conversion)
GOs <- getBM(attributes=c('ensembl_gene_id', 'go_id', 'name_1006', 
                          'go_linkage_type', 'namespace_1003'), 
             filters = 'ensembl_gene_id', 
             values = conversion$ensembl_gene_id, 
             mart = ensembl)
unique(GOs$go_linkage_type)
length(which(GOs$go_linkage_type=="IEA"))
length(unique(GOs$go_id))
