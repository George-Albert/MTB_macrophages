######################
### 0.Dependencies ###
######################
{
  library(tidyverse)
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
  library(scales)
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
# config_dir       <- "1st_Configuration_by_stage"
# config_dir       <- "2nd_Configuration_by_stage"
config_dir       <- "3rd_Configuration_by_stage"

dir_of_gene_list <- file.path(enrichment_dir,config_dir)

#########################
### 3. Load GO lists  ###
#########################

# clusters_filenames <- list.files(path=dir_of_gene_list,"*.txt",recursive = T)
# clusters_filenames <- clusters_filenames[grep("community",clusters_filenames)]
# clusters_list <- lapply(file.path(dir_of_gene_list,clusters_filenames), read.table)
# names(clusters_list) <- gsub("/.*","",clusters_filenames)

all_file_name <- list.files(dir_of_gene_list,pattern = ".txt",recursive = T)
file_name <- all_file_name[!grepl("community", all_file_name) & 
                             !grepl("Configuration", all_file_name)]
list_enrichment <- lapply(file.path(dir_of_gene_list,file_name), read.table)
names(list_enrichment) <- gsub("/.*","",file_name)

# list_enrichment_up <- list_enrichment[c(2,4,6,8)]
# list_enrichment_down <- list_enrichment[c(1,3,5,7)]
list_enrichment_up <- list_enrichment[c(2,4,6,7)]
list_enrichment_down <- list_enrichment[c(1,3,5)]

list_enrichment_to_loop <- list(list_enrichment_up,list_enrichment_down)
list_name_to_loop       <- c("List_UP","List_Down")
level_index             <- list(UP=c("Early.Up","Mid.Up","Late.Up","Persistent.UP"),
                                      Down= c("Early.Down","Mid.Down","Late.Down"))

for (i in seq_along(list_enrichment_to_loop)) {
  
  list_enrichment <- list_enrichment_to_loop[[i]]
  list_name <- list_name_to_loop[i]
  ### Keep the columns we want and create a column with number of the cluster
  list_enrichment_simplyfied <- lapply(names(list_enrichment), function(name) {
    df <- list_enrichment[[name]]  # Extraer el data frame
    df <- df[c("GO_ID", "description", "levels", "OR", "p", "fdr", "Min_value_CI", 
               "Max_value_CI", "genes_in_query")]  # Seleccionar columnas
    df$Cluster <-  name  # Agregar columna con el nombre del data frame
    return(df)
  })
  
  ### keep the original names
  names(list_enrichment_simplyfied) <- names(list_enrichment)
  
  ### Stay with the union of the significative enrichments 
  filter_list <- function(x, th_qB=0.2, th_qOR=0.95,th_qp=0.05) {
    
    th_OR <- quantile(x$OR,th_qOR)
    th_p  <- quantile(x$p,th_qp)
    
    x$hit=0
    x$hit[which(x$OR > th_OR | x$p < th_p)] = 1
    # print(paste0(names(x)," ", length(which(x$hit==1))," th_OR=",th_OR," th_p=",th_p))
    x=x[order(-x$OR),]
    x$Rank_OR=c(1:nrow(x))
    x=x[order(x$p),]
    x$Rank_p=c(1:nrow(x))
    x$Rank_all=sqrt(x$Rank_p^2+x$Rank_OR^2)
    x=x[order(x$Rank_all),]
    return(head(x, 10))
  }
  
  # OR_values <- rep(0.5,3) 
  # fdr_values <-rep(0.01,3) 
  
  ### I stay with the terms that fullfill the criteria
  lista_filtered <- Map(filter_list, list_enrichment_simplyfied)
  lista_filtered <- Filter(function(df) nrow(df) > 0, lista_filtered)
  lista_filtered <- lapply(lista_filtered, function(df) df[order(-df$OR), ])
  
  lista_filtered <- lista_filtered[match(level_index[[i]], names(lista_filtered))]
  
  tab <- lapply(lista_filtered, function(x) (sum(nrow(x))))
  total_sum <- sum(unlist(tab))
  print(total_sum)
  
  ### join the dataframes by conditions
  union_of_df <- bind_rows(lista_filtered)
  row_order_index <- union_of_df[,c("description","Cluster")]
  row_order_index$Cluster <- factor(row_order_index$Cluster,levels = level_index[[i]])
  row_order_index <- row_order_index[order(row_order_index$Cluster),]
  rownames(row_order_index) <- row_order_index$description
  
  ### remove the duplicated terms sort the data
  total_terms <- unique(union_of_df$GO_ID)
  duplicated_terms <- union_of_df[duplicated(union_of_df$GO_ID),]
  
  # Lets create the stages and columns of the heatmap
  # list_enrichment_to_clean <- list_enrichment_simplyfied
  
  ### Filter the list by GO_ID
  list_filtered_clean <- lapply(list_enrichment_simplyfied,function(x){
    l <- subset(x, GO_ID %in% total_terms) #return the elements of GO_ID that are in total terms
    # l <-l[order((l$GO_ID)),]
    l <- l[match(total_terms, l$GO_ID), ]
    return(l)})
  
  OR_list <- lapply(list_filtered_clean, function(x)
  {
    # x[["OR"]][is.infinite(x[["OR"]])] <- 1e100
    x <- x[c("description","OR","Min_value_CI","Max_value_CI","fdr")]
    # x <- x[,1:2]
    # as.matrix(x)
    return(x)
  }) 
  
  # Combine all data frames in OR_list using cbind
  # df <- cbind(OR_list[[1]],OR_list[[2]][2:5],OR_list[[3]][2:5],OR_list[[4]][2:5])
  df <- Reduce(function(x, y) cbind(x, y[, 2:5]), OR_list)
  rownames(df)<- df$description
  
  colnames(df)[1:5] <- paste0(colnames(OR_list[[1]])[1:5],"_", names(OR_list[1]))
  
  colnames(df) <- c(
    colnames(df)[1:5], # Colnames from the first element already defined above
    unlist(lapply(2:length(OR_list), function(i) {# generate the nemes for the rest of elements
      paste0(colnames(OR_list[[i]])[2:5],"_", names(OR_list[i]))
    }))
  )
  
  
  # saveRDS(lista_filtered_up,file=file.path(input_folder,paste0("lista_up_GO_enrichment_OR_",OR,"_fdr_",fdr,".rds")))
  # write.table(df,file=file.path(input_folder,paste0("df_up_GO_enrichment_OR_",OR,"fdr_",fdr,".txt")))
  # name_to_save_xlsx <- c("")
  # write.xlsx(lista_filtered_down, file=file.path(input_folder,paste0("lista_down_GO_enrichment_OR_",OR,"_fdr_",fdr,".xlsx")),
  #            sheetName = name_to_save_xlsx)
  
  
  ### Save the max and min values distinct from Inf
  colmns_to_loop <- (ncol(df))
  max_values_df <- apply(df[,2:colmns_to_loop], 2, function(x) max(x[is.finite(x)], na.rm = TRUE))
  min_values_df <- apply(df[,2:colmns_to_loop], 2, function(x) min(x[is.finite(x)], na.rm = TRUE))
  
  # for (col in colnames(df[,2:colmns_to_loop])) {
  #   
  #   replace_value <- max_values_df[col] + 2
  #   df[is.infinite(df[, col]), col] <- replace_value
  #   
  # }
  
  ###############
  ### Heatmap ###
  ###############
  df_to_select <- df[,grepl("OR",colnames(df)[1:4])]
  dcols(df_to_select)
  df_to_select <- df_to_select[,1:4]
  # df_to_select <- df_to_select[rownames(row_order_index),]
  ### Apply hclust function
  hclust_matrix <- as.matrix(df_to_select)
  hclust_matrix <- apply(hclust_matrix, MARGIN=c(1, 2), FUN =function(x) log2(x))
  
  max_values_hclust <- apply(hclust_matrix, 2, function(x) max(x[is.finite(x)]))
  min_values_hclust <- apply(hclust_matrix, 2, function(x) min(x[is.finite(x)]))
  
  #remove any infinite value
  # hclust_matrix <- hclust_matrix[!apply(hclust_matrix, 1, function(x) any(is.infinite(x))), ]
  
  hclust_matrix[(is.infinite(hclust_matrix)) & hclust_matrix > 0] <-  100
  hclust_matrix[(is.infinite(hclust_matrix)) & hclust_matrix < 0] <- -100
  
  dist_matrix <- dist(hclust_matrix)
  OR_clust <- hclust(dist_matrix,)
  col_dend <- hclust(dist(t(hclust_matrix)))
  # columns_name <- names(list_enrichment_to_clean)
  # columns_name <- c("Early.UP","Late.UP","Mid.UP","Persistent.UP")
  columns_name <- level_index[[i]]
  colnames(hclust_matrix) <- gsub("OR_","",colnames(hclust_matrix))
  colnames(hclust_matrix)
  hclust_matrix <- hclust_matrix[, columns_name]
  colnames(hclust_matrix)
  
  # level_to_use <- c("Early.UP","Mid.UP","Late.UP","Persistent.UP")
  level_to_use <- level_index[[i]]
  
  columns_name <- factor(colnames(hclust_matrix),levels=level_to_use)
  column_title = "Clusters"
  clust_name <- ""
  
  summary(hclust_matrix)
  # color_breaks <- c(-100,-3, 0, 5, 100)
  # color_breaks <- c(-2, 0, 3.5)
  color_breaks <- seq(-2, 3.5, length.out = 11)  # Ajusta según tus datos
  my_palette <- rev(brewer.pal(11, "RdBu")) 
  
  col_fun = colorRamp2(color_breaks, my_palette)
  # col_fun = colorRamp2(color_breaks, c("blue","white","red"))
  layer_fun <- function(j, i, x, y, w, h, fill) {
    actual_values <- pindex(hclust_matrix,as.numeric(i),as.numeric(j))
    # Etiquetar Inf y -Inf
    is_inf <- actual_values == 100
    is_neg_inf <- actual_values == -100
    
    # Dibujar texto para Inf y -Inf
    if (any(is_inf)) {
      grid.rect(x[is_inf], y[is_inf], width = w * 0.9, height = h * 0.9, 
                gp = gpar(fill = "#67001F", col = NA))
      grid.text("Inf", x[is_inf], y[is_inf], gp = gpar(fontsize = 4, col = NA))
    }
    if (any(is_neg_inf)) {
      grid.rect(x[is_neg_inf], y[is_neg_inf], width = w * 0.9, height = h * 0.9, 
                gp = gpar(fill = "#053061", col = NA))
      grid.text("-Inf", x[is_neg_inf], y[is_neg_inf], gp = gpar(fontsize = 4, col = NA))
    }
    # Dibujar texto para otros valores
    grid.text(
      label = ifelse(is_inf | is_neg_inf, "", round(actual_values, digits = 2)),
      x = x, y = y,
      gp = gpar(fontsize = 3))
  }
  k=3
  ht_plt <- Heatmap(hclust_matrix,
                    na_col = "grey2",
                    col = col_fun,
                    # split = k,
                    name="Log2(OR)",
                    column_order = levels(columns_name),
                    show_column_names = T,
                    column_names_gp = gpar(fontsize = 6),
                    row_names_gp = gpar(fontsize = 6),
                    column_title = column_title,
                    column_title_side = "bottom",
                    row_order = rownames(hclust_matrix),
                    # row_split = split,
                    # row_dend_reorder=T,
                    border_gp = gpar(col = "black", lty = 2),
                    # heatmap_height = unit(6, "cm"),
                    # heatmap_width = unit(8, "cm"),
                    width=unit(2.5, "cm"),
                    # show_column_dend = T,
                    # column_dend_side = "top",
                    # cluster_rows = color_branches(OR_clust,k=k,groupLabels = T),
                    # cluster_columns = color_branches(col_dend),
                    row_title = clust_name,
                    row_title_gp = gpar(fontize = 2),
                    row_dend_side = "right",
                    row_names_side = "left",
                    # row_dend_width = unit(2, "cm"),
                    show_row_names = T ,
                    show_row_dend = T,
                    layer_fun = layer_fun,
                    # cell_fun = function(j, i, x, y, width, height, fill) {
                    #   grid.text(round(hclust_matrix[i, j],digits = 2), x, y, 
                    #             gp = gpar(fontsize = 3))},
                    heatmap_legend_param = list(title = "Log2(OR)",
                                                title_position = "leftcenter-rot",
                                                labels_gp = gpar(font = 3),
                                                title_gp = gpar( fontsize = 8)))
  
  draw(ht_plt)
  
  size       <- calc_ht_size(ht_plt)
  cluster    <- dendextend:::cutree(OR_clust, k=k,order_clusters_as_data = F)
  cluster_df <- data.frame(cluster)
  idx        <- match(rownames(hclust_matrix), names(cluster))
  hcl_matrix <- cbind(hclust_matrix,cluster=cluster[idx])
  hclust_df  <- data.frame(hcl_matrix)
  hclust_df  <- hclust_df[order(hclust_df$cluster),]
  
  #  OR <- 2
  #  fdr <- 0.01
  
  # write.table(hclust_df,file = file.path(input_folder,paste0("Cluster_by_clusters_OR_",OR,"_fdr_",fdr,".txt")))
  # write.table(hclust_df,file = file.path(input_folder,paste0("Cluster_by_stage.txt")))
  write.table(hclust_df,file = file.path(input_folder,paste0("Cluster_by_stage_sortby_OR_p_value_",list_name,"_Persistent.txt")))
  
  dir.create(file.path(output_dir,"001_Figures","005_Enrichment",config_dir),recursive=T,showWarnings =F)
  filename <- file.path(output_dir,"001_Figures","005_Enrichment",config_dir,paste0("Heatmap_GO_by_stage_",list_name,"_Persistent.pdf"))
  
  pdf(file =filename ,width=size[1]+2,height=size[2]+2)
  draw(ht_plt)
  dev.off()
}


