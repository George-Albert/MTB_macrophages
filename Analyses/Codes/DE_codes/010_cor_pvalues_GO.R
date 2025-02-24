

####################
##  Dependencies  ##
####################
{
  library(tidyverse)
  library(xlsx)
  }

#################
##  Functions  ##
#################

dcols <- function(df){data.frame(colnames(df))} 

#################
##  Load Data  ##
#################
output_dir <- "Analyses/Outputs"
input_dir  <- "Analyses/Inputs"

clue_GO_dir      <- "Analyses/Outputs/006_clueGO"
clueGO_tab_early <- read.xlsx(file.path(clue_GO_dir,"Early.Up-1.xls"),sheetIndex = 1)
clueGO_tab_mid   <- read.xlsx(file.path(clue_GO_dir,"Mid.Up-1.xls"),sheetIndex = 1)
clueGO_tab_late  <- read.xlsx(file.path(clue_GO_dir,"Late.Up-1.xls"),sheetIndex = 1)

clueGO_list <- list(clueGO_tab_early=clueGO_tab_early,clueGO_tab_mid=clueGO_tab_mid,
                    clueGO_tab_late=clueGO_tab_late)

config_dir <- "1st_Configuration_by_stage"
my_GO_dir <- paste0("Analyses/Outputs/005_Enrichment_GO/",config_dir)
early_up_1  <- read.xlsx(file.path(my_GO_dir,"Early.Up","Enrichment_GO_Early.Up.xlsx"),sheetIndex = 1)
mid_up_1    <- read.xlsx(file.path(my_GO_dir,"Mid.Up","Enrichment_GO_Mid.Up.xlsx"),sheetIndex = 1)
late_up_1   <- read.xlsx(file.path(my_GO_dir,"Late.Up","Enrichment_GO_Late.Up.xlsx"),sheetIndex = 1)

colnames_vec <- c("GO_ID","description","OR","p","fdr","genes_in_query")
early_up <- early_up_1[,colnames_vec]
mid_up   <- mid_up_1[,colnames_vec]
late_up  <- late_up_1[,colnames_vec]

list_my_GO <- list(early_up=early_up,mid_up=mid_up,late_up=late_up)

# ### Reformat cluego table
# df_clueGO <- clueGO_list[[i]]
# df_clueGO <- df_clueGO %>%
#   mutate(Genes_list = strsplit(str_remove_all(df_clueGO$Associated.Genes.Found, "\\[|\\]"), ", "))%>%
#   mutate(num_genes_per_term = sapply(Genes_list, length))
# df_clueGO <- df_clueGO[,c(1,2,4,5,13,14)]
# colnames(df_clueGO) <- c("GO_ID","description","p","fdr",
#                          "Genes_list","num_genes_per_term")
# df_clueGO <- df_clueGO[order(df_clueGO$GO_ID),]
# 
# ### Reformat my GO table
# df_myGO <- list_my_GO[[i]]
# df_myGO <- df_myGO %>%
#   mutate(Genes_list = strsplit(str_remove_all(genes_in_query, "\\[|\\]"), ", ")) %>%
#   mutate(num_genes_per_term = sapply(Genes_list, length))
# df_myGO <- df_myGO[,-c(6)]
# df_myGO <- df_myGO[order(df_myGO$GO_ID),]
# 
# ### Lets merge by description column
# x <- unique(df_clueGO)
# y <- df_myGO
# 
# merged_df <- merge(x, y, by = "description", suffixes = c("_cluego", "_mi_GO"))
# merged_df <- na.omit(merged_df)

# merged_df1 <- merged_df[which(merged_df$GO_ID_cluego=="GO:0009058"),]
# my_gene_query     <- merged_df1$Genes_list_cluego[[1]]
# clueGO_gene_query <- merged_df1$Genes_list_mi_GO[[1]]
# length(my_gene_query)
# length(clueGO_gene_query)
# genes_ind <- which(my_gene_query %in% clueGO_gene_query)
# genes_in <- my_gene_query[genes_ind]
# genes_not_in_clueGO <- setdiff(my_gene_query, clueGO_gene_query)

list_out   <- list()

for (i in c(1,2,3)) {
  
  
  clueGO_tab <- clueGO_list[[i]] 
  clueGO_tab <- clueGO_tab[,c(1,2,4,5,12)]
  colnames(clueGO_tab) <- c("GO_ID","description","p","fdr",
                            "Associated.Genes.Found")
  # fdr_threshold <- 0.01
  x <- clueGO_tab
  y <- list_my_GO[[i]]
  merged_df <- merge(x, y, by = "description", suffixes = c("_cluego", "_mi_GO"))
  merged_df <- na.omit(merged_df)
  
  merged_df <- merged_df %>% 
    distinct(description, .keep_all = TRUE)
  
  merged_df$log_p_cluego <- -log10(merged_df$p_cluego)
  merged_df$log_p_mi_GO  <- -log10(merged_df$p_mi_GO)
  
  list_out[[i]] <-  merged_df
  names(list_out)[i] <- names(list_my_GO)[i]
  
  cor_test <- cor.test(merged_df$log_p_cluego, merged_df$log_p_mi_GO, method = "spearman")  
  print(cor_test)
  cor_pear_test <- cor.test(merged_df$log_p_cluego, merged_df$log_p_mi_GO, method = "pearson")  
  print(cor_pear_test)
  
  title_vec      <- c("Early.Up","Mid.Up","Late.Up")
  color_for_line <- c("red","blue","green")
  
  cor_plt <- ggplot(merged_df, aes(x = log_p_cluego, y = log_p_mi_GO)) +
    geom_point(alpha = 0.7, color = color_for_line[i]) +
    geom_smooth(method=lm,formula = "y ~ x",show.legend = FALSE) +
    theme_minimal() +
    geom_vline(xintercept = 2.5)+
    geom_hline(yintercept = 0)+
  labs(title = paste0(title_vec[i],
                      "\nR-Spearman=", round(cor_test$estimate, digits = 2),
                      "\nR-Pearson=", round(cor_pear_test$estimate, digits = 2),
                      "\nNum_of_common terms=", dim(merged_df)[1]),
         
      x = "-log10(p-value) ClueGO",
      y = "-log10(p-value) mi_GO") +
    theme(plot.title = element_text(hjust = 0.5))
  cor_plt
  
  dir.create(file.path(output_dir,"001_Figures","006_Cor_GO"),recursive=T,showWarnings =F)
  filename <- file.path(output_dir,"001_Figures","006_Cor_GO",paste0("Correlation_plt_",title_vec[i],".pdf"))
  
  pdf(file =filename ,width=8,height=6)
  print(cor_plt)
  dev.off()
  

}
# 
# ### Inspect the mmerged data frame looking for the log pvalues
# df <- list_out[[2]]
# 
# df <- df[which(df$log_p_mi_GO > 7.0),]
# df <- df[which(df$log_p_cluego < 2.5),]
# 
# ### Test that GO term
# term_test <- late_up[which(late_up$GO_ID=="GO:0009058"),]
# my_gene_query <- term_test$genes_in_query
# 
# my_gene_query     <- df$genes_in_query
# clueGO_gene_query <- data.frame(clueGO_tab$Associated.Genes.Found)
# length(clueGO_gene_query[1,])
# 
# ### Split gene query by genes
# my_gene_query <- strsplit(my_gene_query, ",\\s*")[[1]]
# length(my_gene_query)
# 
# ### Load the node attribute from ClueGo output
# clueGO_node_att <- read.delim(file.path(clue_GO_dir,"Mid.Up_0","Mid.Up_0 NodeAttributeTable.txt"))
# clueGO_term_test <- clueGO_node_att[which(clueGO_node_att$ID=="GO:0009058"),]
# 
# genes_found <- clueGO_term_test$Associated.Genes.Found
# genes_found <- strsplit(gsub("\\[|\\]", "", genes_found), ",\\s*")[[1]]
# length(genes_found)
# 
# genes_associated <- clueGO_term_test$All.Associated.Genes 
# genes_associated <- strsplit(gsub("\\[|\\]", "", genes_associated), ",\\s*")[[1]]
# 
# identical(my_gene_query,genes_found)
# length(setdiff(my_gene_query,genes_found)) # See if there are any different element
# length(setdiff(genes_found,my_gene_query))
# 
# 
# x <- GO_gene_term_matrix_bp[,which(colnames(GO_gene_term_matrix_bp)=="antigen processing and presentation of endogenous peptide antigen via MHC class I")]    
# length(which(x==TRUE))
# x <- x[which(x==TRUE)]
# my_associated_genes <- names(x)
# 
# identical(my_associated_genes,genes_associated)
# setdiff(my_associated_genes,genes_associated) # See if there are any different element
# setdiff(genes_associated,my_associated_genes)
# 

