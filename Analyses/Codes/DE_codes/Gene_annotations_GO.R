library(GO.db)
library("org.Hs.eg.db")

annotations_ensembl <- mapIds(org.Hs.eg.db,
                              keys(org.Hs.eg.db, "GO"),
                              column = "ENSEMBL", 
                              keytype = "GO",
                              multiVals = "list")
annotations_ncbi <- mapIds(org.Hs.eg.db, 
                           keys(org.Hs.eg.db, "GO"),
                           column = "ENTREZID", 
                           keytype = "GO",
                           multiVals = "list")

x <- names(annotations_ncbi)
y <- names(annotations_ensembl)

identical(x,y)

go_descriptions <- data.frame( descriptions=mapIds(
  GO.db,  
  keys = names(annotations_ncbi),  
  column = "TERM",  
  keytype = "GOID",  
  multiVals = "first"  # Solo queremos una descripción única por GO term
))

go_descriptions$GO_ID <- rownames(go_descriptions)


