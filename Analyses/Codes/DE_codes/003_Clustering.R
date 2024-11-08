
#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyr)
  library(ggthemes)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(RColorBrewer)
  library(ggdendro)
  library(dendextend)
  library(dendsort)
  library(reshape2)
  library(umap)
  library(circlize)
  library(ComplexHeatmap)
  ### Clustering packages
  # library(RCKS)
  # library(M3C)
  library(ConsensusClusterPlus)
  library(factoextra)
  library(NbClust)
  library(mclust)
  library(clustertend)

}

### Set dir
main_dir <- getwd()
setwd(main_dir)
input_dir  <- "Analyses/Inputs"
output_dir <- "Analyses/Outputs"

### load data
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
cols=read.table(file.path(input_dir,"002_Processed","ready_for_DE","metadata.txt"))
cols_whole=read.table(file.path(input_dir,"002_Processed","whole","metadata_whole.txt"))
results <- read.table(file.path(output_dir,"003_DE_def","resultados.txt"))

length(which(colnames(reads)!=rownames(cols)))
length(which(paste0(cols$Individual,"_",cols$Setup)!=rownames(cols)))
length(which(rownames(reads)!= rownames(results)))

dir.create(file.path(output_dir,"004_Clustering_plots"),showWarnings = F)
cluster_dir <- file.path(output_dir,"004_Clustering_plots")

### Compute the size of the heatmap
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
##########################################
##  1. Find optimal number of Clusters  ##
##########################################

{
  tab=results
  
  grep("beta",colnames(tab))
  colnames(tab)[c(1,7,13,19,25,31,37,43,49,55,61,67,73)]
  grep("hit",colnames(tab))
  colnames(tab)[grep("hit",colnames(tab))]
  grep("SE",colnames(tab))
  colnames(tab)[grep("SE",colnames(tab))]
  
  # Select these betas
  # "beta_inf_2"         "beta_inf_20"        "beta_inf_48"        "beta_inf_72" 
  
  ### Build table with the infection logFCs of genes that respond to infection
  ###  at at least one timepoint. threshold > 0.05
  
  tab$hit_any=apply(tab[,c("hit_inf_2","hit_inf_20","hit_inf_48","hit_inf_72")],1,function(x){as.numeric(sum(x)>0)})
  tab_hit_labels=tab[which(tab$hit_any==1),c("hit_inf_2","hit_inf_20","hit_inf_48","hit_inf_72")]
  tab=tab[which(tab$hit_any==1),c( "beta_inf_2","beta_inf_20","beta_inf_48","beta_inf_72")]
  
  ### Normalize the betas of each gene time-wise. Do I have to rescale again?
  ttab=t(tab)
  norm_ttab=sapply(1:ncol(ttab),function(x){ttab[,x]/max(abs(ttab[,x]))})
  tab_norm=t(norm_ttab)
  rownames(tab_norm)=rownames(tab)
  
  ### Let us rename these objects to store them more properly/safely
  betas=tab
  norm_betas=data.frame(tab_norm)
  
  set.seed(12345)
  
  ### Elbow method (look at the knee)
  # Elbow method for kmeans
  elb_plot <- fviz_nbclust(norm_betas, kmeans, method = "wss") +
    geom_vline(xintercept = 3, linetype = 2)
  # k=3
  pdf(file.path(cluster_dir,paste0("001_elbow_plt.pdf")),width=12,height=10)
  print(elb_plot)
  dev.off()
  
  # Average silhouette for kmeans
  silo_plt <- fviz_nbclust(norm_betas, kmeans, method = "silhouette")
  #k=2
  pdf(file.path(cluster_dir,paste0("002_silhouette_plt.pdf")),width=12,height=10)
  print(silo_plt)
  dev.off()
  
  ### take a time to run
  clust_res <- NbClust(data = norm_betas, diss = NULL, distance = "euclidean", 
                       min.nc = 2,max.nc = 20, method = "kmeans", index = "all",
                       alphaBeale = 0.1)
  # *** : The Hubert index is a graphical method of determining the number of clusters.
  # In the plot of Hubert index, we seek a significant knee that corresponds to a 
  # significant increase of the value of the measure i.e the significant peak in Hubert
  # index second differences plot. 
  # 
  # *** : The D index is a graphical method of determining the number of clusters. 
  # In the plot of D index, we seek a significant knee (the significant peak in Dindex
  # second differences plot) that corresponds to a significant increase of the value of
  # the measure. 
  # 
  # ******************************************************************* 
  #   * Among all indices:                                                
  #   * 13 proposed 2 as the best number of clusters 
  # * 2 proposed 3 as the best number of clusters 
  # * 1 proposed 4 as the best number of clusters 
  # * 4 proposed 6 as the best number of clusters 
  # * 1 proposed 11 as the best number of clusters 
  # * 1 proposed 12 as the best number of clusters 
  # * 1 proposed 18 as the best number of clusters 
  # * 1 proposed 19 as the best number of clusters 
  # 
  
  clust_res$Best.nc
  ### Is not the better approach
  best.partition.list <- data.frame(clust_res$Best.partition)
  
  ### Define maximum number of clusters
  kmax <- 20
  
  ### k-means clustering with concensus
  data <- as.matrix(norm_betas,nrow=nrow(norm_betas), ncol=ncol(norm_betas))
  data <- na.omit(data)
  tdata <- t(data)
  
  rcc_kmeans = ConsensusClusterPlus(tdata,maxK=kmax,reps=100,pItem=0.8,pFeature=1,
                              title="kmeans_concensus",distance="euclidean",clusterAlg="km",
                              plot= "pngBMP")
  ### hierarchical clustering with concensus
  rcc_hc = ConsensusClusterPlus(tdata,maxK=kmax,reps=100,pItem=0.8,pFeature=1,
                             title="hclust_concensus",distance="pearson",clusterAlg="hc",
                             plot="pngBMP")
  
  set.seed(12345)
  
  #################################
  ## Mclust clustering algorithm ##
  #################################
  
  mc <- Mclust(norm_betas) # Model-based-clustering
  summary(mc) # Print a summary
  
  mc$modelName # Optimal selected model ==> "VVV"
  mc$G # Optimal number of cluster => 3
  head(mc$z, 30) # Probality to belong to a given cluster
  head(mc$classification, 30) # Cluster assignement of each observation
  
  cluster_out <- data.frame(norm_betas,Mclust_7k=mc$classification)
  
  # BIC values used for choosing the number of clusters
  BIC_plt <- fviz_mclust(mc, "BIC", palette = "jco")
  
  pdf(file.path(cluster_dir,paste0("008_bic_Mclust.pdf")),width=12,height=10)
  print(BIC_plt)
  dev.off()
  
  # Classification: plot showing the clustering
  pca_mclust <- fviz_mclust(mc, "classification", geom = "point",
              pointsize = 1.5, palette = "jco")
  pdf(file.path(cluster_dir,paste0("009_pca_Mclust.pdf")),width=12,height=10)
  print(pca_mclust)
  dev.off()
  # Classification uncertainty
  uncert_plt <- fviz_mclust(mc, "uncertainty", palette = "jco")
  pdf(file.path(cluster_dir,paste0("010_uncertainty_Mclust.pdf")),width=12,height=10)
  print(uncert_plt)
  dev.off()
  
  ##############################################
  ##  Run kmeans for kmax number of clusters  ##
  ##############################################
  ### Create a list to save the kmeans results
  kmeans_res <- list()
  cluster_list <- list()
  kmax <- 20
  ### vector to loop
  vec <- seq(kmax)
  
  for (number_clusters in vec) {
    
    kmeans_res[[number_clusters]] <- kmeans(tab[, 1:4], centers=number_clusters,
                                            iter.max = 100, 
                                            nstart = number_clusters)
    cluster_list[[number_clusters]] <-kmeans_res[[number_clusters]]$cluster
    
  }
  ### check the order of the gene names in each data
  length(which(rownames(norm_betas) != names(kmeans_res[[6]]$cluster)))
  
  cluster_out <- data.frame(cluster_out,kmeans_cluster_7k=kmeans_res[[7]]$cluster,
                          kmeans_cluster_12k=kmeans_res[[12]]$cluster,
                          kmeans_cluster_16k=kmeans_res[[16]]$cluster)
  
  ### Check the size of the cluster partitions
  kmeans_res[[7]]$size
  # 482  374   53 1667 2128  187
  kmeans_res[[12]]$size
  # [1]   72  319 1005  303  233  140    6   39  566  861   74 1273
  kmeans_res[[16]]$size
  # 133  74 513 184  39 646 244 255   6  42 872 115 214 247 742 565
  
  ### Plot and Save PCA for kmean and Hclust dendrogram for 3 k number of clusters
  
  # Define vector with the 3 k we are going to use
  k_clusters <- c(k_7=7,k_12=12,k_16=16)
  
  for (k in k_clusters) {
    
    print(k)
    ### PCA for kmean
    pca_kmean <- fviz_cluster(kmeans_res[[k]], norm_betas,stand = F,axes = c(1,2),geom="point",
                              ellipse.alpha = 0.1,ellipse.type = "norm",ggtheme = theme_classic())

    ### Hierarchical clustering
    hclust_dend <- fviz_dend(hclust(dist(norm_betas)), k = k, k_colors = "jco",
                             as.ggplot = TRUE, show_labels = FALSE)
    
    pdf(file.path(cluster_dir,paste0("003_pca_kmean_k=",k,".pdf")),width=12,height=10)
    print(pca_kmean)
    dev.off()
    
    pdf(file.path(cluster_dir,paste0("004_hclust_dendrogram_k=",k,".pdf")),width=12,height=10)
    print(hclust_dend)
    dev.off()
    
  }
  
  
  #####################################
  ##  Cluster Tendency with Hopkins  ##
  #####################################
  
  
  # If the value of Hopkins statistic is close to zero, then we can reject 
  # the null hypothesis and conclude that the data set D is significantly 
  # a clusterable data.
 
  ###Compute Hopkins statistic for log fold changes
  # set.seed(12345)
  hopkins(norm_betas, n = nrow(norm_betas)-1)
  
  # About 0.09 tho, is clusterizable
  
  ### Dissimilarity matrix
  dis_mat <- fviz_dist(dist(norm_betas), show_labels = FALSE)+
    labs(title = "TB infection over time (DM)")
 
  ggsave(
    filename = file.path(cluster_dir, "005_Dissimilarity_matrix.png"),
    plot = dis_mat,
    width = 12, 
    height = 10, 
    dpi = 300  # Adjust the dpi resolution
  )
 
  #####################################
  ## Hierarchichal k-means algorithm ##
  #####################################
 
  hckmeans_res <- list()
  # rect_colors <-colorRampPalette(brewer.pal(12, "Set3"))(16)
  
  for(k in c(1:3)){
    
    hckmeans_res[[k]] <- hkmeans(
      norm_betas,
      k=k_clusters[k],
      hc.metric = "euclidean",
      hc.method = "ward.D2",
      iter.max = 10,
      km.algorithm = "Hartigan-Wong"
    )
      ### Create a dendrogram
      rect_colors  <-  colorRampPalette(brewer.pal(12, "Set3"))(k)
      
      # dend_hckmean <-  hkmeans_tree(hckmeans_res[[k]], rect.col = rect_colors)
      # Visualize the tree
      dend_hckmean <-fviz_dend(hckmeans_res[[k]], cex = 0.6, palette = "jco",show_labels = F,
                rect = TRUE, rect_border = "jco", rect_fill = TRUE)

      pca_hckmean <- fviz_cluster(hckmeans_res[[k]], norm_betas,stand = F,axes = c(1,2),geom="point",
                                  ellipse.alpha = 0.1,ellipse.type = "norm",ggtheme = theme_classic())
      
      pdf(file.path(cluster_dir,paste0("006_pca_hckmean_k=",k_clusters[k],".pdf")),width=12,height=10)
      print(pca_hckmean)
      dev.off()
      
      pdf(file.path(cluster_dir,paste0("007_hckmean_dendrogram_k=",k_clusters[k],".pdf")),width=12,height=10)
      print(dend_hckmean)
      dev.off()
    
  }
  safe <- cluster_out
  cluster_out <- data.frame(cluster_out,hckmeans_7k=hckmeans_res[[1]]$cluster,
                            hckmeans_12k=hckmeans_res[[2]]$cluster,
                            hckmeans_16k=hckmeans_res[[3]]$cluster)
  
  cluster_input_dir <- file.path(input_dir,"003_Clustering_tabs")
  dir.create(cluster_input_dir,recursive = T,showWarnings = F)
  
  saveRDS(hckmeans_res, file = file.path(cluster_input_dir,"001_hckmeans.RDS")) 
 
  write.table(cluster_out,file.path(cluster_input_dir,"lfc_with_clustering.txt"))
  
  ###############################################################
  # for the number of clusters selected, get the k-means object##
  ###############################################################
  cluster_out <- read.table(file.path(cluster_input_dir,"lfc_with_clustering.txt"))
  
  col_fun = colorRamp2(c(-1, 0, 1), c("#FFA373", "white","#50486D"))
  col_fun(seq(-1, 1))
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white","red"))
  # cluster_colors <- c("1" = "#FF6347",   # Tomato
  #                     "2" = "#4682B4",   # SteelBlue
  #                     "3" = "#32CD32",   # LimeGreen
  #                     "4" = "#FFD700",   # Gold
  #                     "5" = "#8A2BE2",   # BlueViolet
  #                     "6" = "#FF1493",   # DeepPink
  #                     "7" = "#2E8B57")   # SeaGreen
  
  cluster_colors <- c("1" = "#FF6347",   # Tomato
                      "2" = "#4682B4",   # SteelBlue
                      "3" = "#32CD32",   # LimeGreen
                      "4" = "#FFD700",   # Gold
                      "5" = "#8A2BE2",   # BlueViolet
                      "6" = "#FF1493",   # DeepPink
                      "7" = "#2E8B57",   # SeaGreen
                      "8" = "#D2691E",   # Chocolate
                      "9" = "#A52A2A",   # Brown
                      "10" = "#00FFFF",  # Cyan
                      "11" = "#FF00FF",  # Magenta
                      "12" = "#800000",  # Maroon
                      "13" = "#808000",  # Olive
                      "14" = "#FF8C00",  # DarkOrange
                      "15" = "#C71585",  # MediumVioletRed
                      "16" = "#000080")  # Navy
  vec <- 5:ncol(cluster_out)
  
  for(i in vec){
    
    lfc_normalized <- cluster_out[,c(1:4,i)]
    
    lfc_normalized[,5]=factor(lfc_normalized[,5],levels=unique(sort(lfc_normalized[,5])))
    
    table=lfc_normalized[order(lfc_normalized[,5]),]
    hclust_matrix <- table[1:4] %>% 
      as.matrix()
    columns_name <- colnames(table)[1:4]
    column_title = "Conditions"
    clust_name <- "Genes"
    cluster_colors_ha <- cluster_colors[1:max(as.numeric(table[,5]))]
    split <- as.numeric(table[,5])
    order_rows <- rownames(table)
    length(which(rownames(hclust_matrix)!=rownames(table)))
    
    ha <- rowAnnotation(Cluster = as.factor(table[,5]),col = list(Cluster = cluster_colors_ha),
                        annotation_legend_param = list(
                          title = "Cluster", 
                          labels = unique(table[,5])  # Esto incluye los números de los clusters
                        ),
                        border = TRUE)
    
    ht_plt <- Heatmap(hclust_matrix,
                      na_col = "grey2",
                      col = col_fun,
                      split = split,
                      # row_split = split,
                      row_order = order_rows,
                      name="expression levels",
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
                      width=unit(3, "cm"),
                      # show_column_dend = T,
                      # column_dend_side = "top",
                      # cluster_rows = color_branches(OR_clust,k=k),
                      # cluster_columns = color_branches(col_dend),
                      row_title = clust_name,
                      row_title_gp = gpar(fontize = 2),
                      # row_dend_side = "right",
                      # row_names_side = "left",
                      # row_dend_width = unit(2, "cm"),
                      show_row_names = F ,
                      show_row_dend = F,
                      right_annotation = ha,
                      # layer_fun = function(j, i, x, y, width, height, fill) {
                      #   grid.text(round(hclust_matrix[i, j],digits = 2), x, y, gp = gpar(fontsize = 3))},
                      heatmap_legend_param = list(title = "Normalized Exp Levels",
                                                  title_position = "leftcenter-rot",
                                                  labels_gp = gpar(font = 3),
                                                  title_gp = gpar( fontsize = 8)))
    
    draw(ht_plt)
    
    ### Compute the size of the htmap
    size <- calc_ht_size(ht_plt)
    
    ### Save the data
    filename <- file.path(cluster_dir,paste0("011_Heatmap_normalized_",colnames(table)[5],".pdf"))
    
    pdf(file =filename ,width=size[1]+1,height=size[2]+2)
    draw(ht_plt,ht_gap=unit(0.5,"mm"))
    dev.off()
    
    
  }
  
  
  }




