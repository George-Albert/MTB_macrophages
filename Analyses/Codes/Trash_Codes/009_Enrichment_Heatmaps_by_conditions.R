######################
### 0.Dependencies ###
######################
{
  library(qvalue)
  library(igraph)
  library(xlsx)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(ggsci)
  library(ggdendro)
  library(dendextend)
  library(dendsort)
}

###########################
### 1.Declare Functions ###
###########################

dcols=function(x){data.frame(colnames(x))}
calc_ht_size = function(ht, unit = "inch") {
        pdf(NULL)
        ht = draw(ht)
        w = ComplexHeatmap:::width(ht)
        w = convertX(w, unit, valueOnly = TRUE)
        h = ComplexHeatmap:::height(ht)
        h = convertY(h, unit, valueOnly = TRUE)
        dev.off()
        
        c(w, h)
}

#################################
### 2. Set working directory  ###
#################################

main_wd <- getwd()
setwd(main_wd)
input_folder     <- "Analyses/Inputs/005_Enrichment/GO_terms"
output_dir       <- "Analyses/Outputs"
enrichment_dir   <- file.path(output_dir,"005_Enrichment_GO")
config_dir       <- "test_Configuration_by_condition"
dir_of_gene_list <- file.path(enrichment_dir,config_dir)

#########################
### 3. Load GO lists  ###
#########################

# clusters_filenames <- list.files(path=dir_of_gene_list,"*.txt",recursive = T)
# clusters_filenames <- clusters_filenames[grep("community",clusters_filenames)]
# clusters_list <- lapply(file.path(dir_of_gene_list,clusters_filenames), read.table)
# names(clusters_list) <- gsub("/.*","",clusters_filenames)

all_file_name <- list.files(dir_of_gene_list,pattern = ".txt",recursive = T)
file_name <- all_file_name[!grepl("community", all_file_name)]
list_enrichment <- lapply(file.path(dir_of_gene_list,file_name), read.table)
names(list_enrichment) <- gsub("/.*","",file_name)

# dcols(list_enrichment$`C1_C2_Glycerol_Dextrose/Enhanced_Down/Enrichment_GO.txt`)
list_enrichment_simplyfied <- lapply(list_enrichment, function(x) {
  x[c("GO_ID","description","levels","OR","p","fdr","Min_value_CI",
  "Max_value_CI","genes_in_query")]
  }) 

### Stay with the union of the significative enrichments 
filter_list <- function(x){
  
  df <- x[which(x$OR > OR & x$fdr < fdr),]
  
  return(df)
}

OR <- 3
fdr <- 0.01

### I stay with the terms that fullfill the criteria
lista_filtered <- lapply(list_enrichment_simplyfied,filter_list)

### join the dataframes by conditions
union_of_df <- bind_rows(lista_filtered)
### remove the duplicated terms sort the data
total_terms <- unique(union_of_df$GO_ID)
total_terms <- total_terms[order(total_terms)]

# list_filtered_clean <- lista_filtered
list_filtered_clean <- lapply(list_enrichment_simplyfied,function(x){
  l <- subset(x, GO_ID %in% total_terms)
  l <-l[order(l$GO_ID),]
  return(l)})

OR_list <- lapply(list_filtered_clean, function(x)
{
  # x[["OR"]][is.infinite(x[["OR"]])] <- 1e100
  x <- x[c("description","OR","Min_value_CI","Max_value_CI","fdr")]
  # x <- x[,1:2]
  # as.matrix(x)
  return(x)
}) 

df <- cbind(OR_list[[1]],OR_list[[2]][2:5],OR_list[[3]][2:5],OR_list[[4]][2:5])
colnames(df) <- c("description_infected_20h",
                  "OR_infected_20h","lower_CI_infected_20h","upper_CI_infected_20h","FDR_infected_20h",
                  "OR_infected_2h","lower_CI_infected_2h","upper_CI_infected_2h","FDR_infected_2h",
                  "OR_infected_48h","lower_CI_infected_48h","upper_CI_infected_48h","FDR_infected_48h",
                  "OR_infected_72h","lower_CI_infected_72h","upper_CI_infected_72h","FDR_infected_72h")
rownames(df)      <- df$description

# saveRDS(lista_filtered_up,file=file.path(input_folder,paste0("lista_up_GO_enrichment_OR_",OR,"_fdr_",fdr,".rds")))
# write.table(df,file=file.path(input_folder,paste0("df_up_GO_enrichment_OR_",OR,"fdr_",fdr,".txt")))
# name_to_save_xlsx <- c("")
# write.xlsx(lista_filtered_down, file=file.path(input_folder,paste0("lista_down_GO_enrichment_OR_",OR,"_fdr_",fdr,".xlsx")),
#            sheetName = name_to_save_xlsx)

  
  ### Save the max and min values distinct from Inf
  colmns_to_loop <- (ncol(df))
  max_values_df <- apply(df[,2:colmns_to_loop], 2, function(x) max(x[is.finite(x)], na.rm = TRUE))
  min_values_df <- apply(df[,2:colmns_to_loop], 2, function(x) min(x[is.finite(x)], na.rm = TRUE))
  
  for (col in colnames(df[,2:colmns_to_loop])) {
    
    replace_value <- max_values_df[col] + 2
    df[is.infinite(df[, col]), col] <- replace_value
    
  }
  
###############
### Heatmap ###
###############

### Apply hclust function
hclust_matrix <- as.matrix(df[,c(2,6,10,14)])
# hclust_matrix[is.infinite(hclust_matrix)] <- 1e100
hclust_matrix <- apply(hclust_matrix, MARGIN=c(1, 2), FUN = log2)
max_values_hclust <- apply(hclust_matrix, 2, function(x) max(x[is.finite(x)]))
hclust_matrix[(is.infinite(hclust_matrix))] <- -1
dist_matrix <- dist(hclust_matrix)
OR_clust <- hclust(dist_matrix,method="ward.D2")
col_dend <- hclust(dist(t(hclust_matrix)))
columns_name <- names(list_enrichment)
colnames(hclust_matrix) <- columns_name
columns_name <- factor(colnames(hclust_matrix),levels=c("infected_2h",
                                                        "infected_20h",
                                                        "infected_48h",
                                                        "infected_72h"))
column_title = "Conditions"
clust_name <- "GO"

summary(hclust_matrix)
color_breaks <- c(-3, 0, 4)
# my_palette <- c( "yellow",
#                  "blue",
#                  colorRampPalette(rev(brewer.pal(8, "Spectral")))(n = length(color_breaks[2:4])-1),
#                  "red")
col_fun = colorRamp2(color_breaks, c("blue","white","red"))
# col_fun = colorRamp2(color_breaks, my_palette)
k=16
ht_plt <- Heatmap(hclust_matrix,
                  na_col = "grey2",
                  col = col_fun,
                  # split = k,
                  name="Log2(OR)",
                  column_order = columns_name,
                  show_column_names = T,
                  column_names_gp = gpar(fontsize = 6),
                  row_names_gp = gpar(fontsize = 6),
                  column_title = column_title,
                  column_title_side = "bottom",
                  # row_dend_reorder=T,
                  border_gp = gpar(col = "black", lty = 2),
                  # heatmap_height = unit(6, "cm"),
                  # heatmap_width = unit(8, "cm"),
                  width=unit(2, "cm"),
                  # show_column_dend = T,
                  # column_dend_side = "top",
                  cluster_rows = color_branches(OR_clust,k=k,groupLabels = T),
                  # cluster_columns = color_branches(col_dend),
                  row_title = clust_name,
                  row_title_gp = gpar(fontize = 2),
                  row_dend_side = "right",
                  row_names_side = "left",
                  # row_dend_width = unit(2, "cm"),
                  show_row_names = T ,
                  show_row_dend = T,
                  layer_fun = function(j, i, x, y, w, h, fill) {
                    actual_values <- pindex(hclust_matrix,as.numeric(i),as.numeric(j))
                    grid.text(
                      label = round(actual_values, digits = 2),
                      x = x, y = y,gp = gpar(fontsize = 3)
                    )},
                  # cell_fun = function(j, i, x, y, width, height, fill) {
                  #   grid.text(round(hclust_matrix[i, j],digits = 2), x, y, 
                  #             gp = gpar(fontsize = 3))},
                  heatmap_legend_param = list(title = "Log2(OR)",
                                              title_position = "leftcenter-rot",
                                              labels_gp = gpar(font = 3),
                                              title_gp = gpar( fontsize = 8)))

draw(ht_plt)

size <- calc_ht_size(ht_plt)
cluster <- dendextend:::cutree(OR_clust, k=k,order_clusters_as_data = F)
cluster_df <- data.frame(cluster)
idx <- match(rownames(hclust_matrix), names(cluster))
hcl_matrix <- cbind(hclust_matrix,cluster=cluster[idx])
hclust_df <- data.frame(hcl_matrix)
hclust_df <- hclust_df[order(hclust_df$cluster),]

OR <- 3
fdr <- 0.01

write.table(hclust_df,file = file.path(input_folder,paste0("Cluster_by_condition_OR_",OR,"_fdr_",fdr,".txt")))

dir.create(file.path(output_dir,"001_Figures","005_Enrichment"),recursive=T,showWarnings =F)
filename <- file.path(output_dir,"001_Figures","005_Enrichment",paste0("Heatmap_GO_by_condition_OR_",OR,".pdf"))

pdf(file =filename ,width=size[1]+3,height=size[2]+4)
draw(ht_plt)
dev.off()

  # ha = HeatmapAnnotation(pt = anno_points(1:61))
  # ha = rowAnnotation(pt = anno_points(1:61))
  
# ht_plt + rowAnnotation(foo = anno_block(
#   panel_fun = function(index, levels) {
#     grid.rect(gp = gpar(col = "black"))
#     ggplot(aes)
#   },
#   width = unit(2, "cm")))

#################################
### Trends by clusters Fig.3G ###
#################################
# 
# hclust_df_up <- read.table(file = file.path(input_folder,"Cluster_OR_4_fdr_0.05_up_k=6.txt"))
# 
# # Compute the mean by cluster
# result_down <- hclust_df_up %>%
#   group_by(cluster) %>%  # Group by column "cluster"
#   summarize(
#     Mean_G_D = mean(G_D),
#     Mean_G = mean(G),
#     Mean_G_LCFA = mean(G_LCFA)
#   )
# 
# ### inspect the result
# print(result_up)
# 
# write.table(result_down,file = file.path(input_folder,"table_down_mean_per_cluster.txt"))
# write.table(result_up,file = file.path(input_folder,"table_up_mean_per_cluster.txt"))
# 
# 
# data_up <- pivot_longer(result_up,cols = !cluster,names_to = "Conditions",values_to = "Mean")
# data_up$Conditions <- factor(data_up$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )
# 
# data_down <- pivot_longer(result_down,cols = !cluster,names_to = "Conditions",values_to = "Mean")
# data_down$Conditions <- factor(data_down$Conditions,levels = c("Mean_G_D", "Mean_G", "Mean_G_LCFA") )
# 
# myarrow=arrow(angle = 15, ends = "last", type = "closed")
# col <- c("#852121", "#8f4a17", "#2c682c")
# 
# list_cluster_down <- split(data_down,data_down$cluster)
# list_cluster_up <- split(data_up,data_up$cluster)
# 
# 
# plot_and_save <- function(dataframe, filename) {
#   
#   mean_plt_fun <- ggplot(data=dataframe,aes(x=Conditions,y=Mean, colour=Conditions,group=1))+
#     geom_line(color="black",arrow=myarrow)+
#     geom_point(size=4)+
#     scale_color_manual(name="Condition",values = col)+
#     ylab("Mean")+
#     ggtitle("")+
#     theme_classic()
#     
#   print(mean_plt_fun)
#     
#   pdf(file =filename)
#   print(mean_plt_fun)
#   dev.off()
# }
# 
# for (i in 1:length(list_cluster_down)) {
#   filename <- file.path(output_dir,"Figures_paper",paste("plot_mean_k_", i,"_down.pdf"))
#   plot_and_save(list_cluster_down[[i]], filename)
# }
# 
# for (i in 1:length(list_cluster_up)) {
#   filename <- file.path(output_dir,"Figures_paper",paste("plot_mean_k_", i,"_up.pdf"))
#   plot_and_save(list_cluster_up[[i]], filename)
# }

