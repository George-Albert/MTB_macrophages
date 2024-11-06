
#########################
## 0.Load Dependencies ##
#########################

{
  library(ggplot2)
  library(ggrepel)
  library(limma)
  library(edgeR)
  library(qvalue)
  library(cowplot)
  library(stats)
  library(openxlsx)
  library(ngram)
  #library(rtracklayer)
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
  library(M3C)
  library(ConsensusClusterPlus)
  library(factoextra)
  library(NbClust)
  library(mclust)
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
  fviz_nbclust(norm_betas, kmeans, method = "wss") +
    geom_vline(xintercept = 3, linetype = 2)
  # k=3
  
  # Average silhouette for kmeans
  fviz_nbclust(norm_betas, kmeans, method = "silhouette")
  #k=2
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
                              plot= "pdf")
  ### hierarchical clustering with concensus
  rcc_hc = ConsensusClusterPlus(tdata,maxK=kmax,reps=100,pItem=0.8,pFeature=1,
                             title="hclust_concensus",distance="pearson",clusterAlg="hc",
                             plot="pdf")
  
  ### Create a list to save the kmeans results
  kmenas_res <- list()
  cluster_list <- list()
  ### vector to loop
  vec <- seq(kmax)
  
  ### Run kmeans for kmax number of clusters
  for (number_clusters in vec) {
    
    kmeans_res[[number_clusters]] <- kmeans(tab[, 1:4], centers=number_clusters,
                                            iter.max = 100, 
                                            nstart = number_clusters)
    cluster_list[[number_clusters]] <-kmenas_res[[number_clusters]]$cluster
    
  }
  ### check the order of the gene names in each data
  length(which(rownames(norm_betas) != names(kmenas_res[[6]]$cluster)))
  
  kmeans_data <- data.frame(norm_betas,kmeans_cluster_6k=kmenas_res[[6]]$cluster,
                          kmeans_cluster_12k=kmenas_res[[12]]$cluster,
                          kmeans_cluster_16k=kmenas_res[[16]]$cluster)
  
  ### Check the size of the cluster partitions
  kmeans_res[[6]]$size
  # 482  374   53 1667 2128  187
  kmeans_res[[12]]$size
  # 303  140 1273   39   74  566  233  319 1005   72    6  861
  kmeans_res[[16]]$size
  # 133  74 513 184  39 646 244 255   6  42 872 115 214 247 742 565
  
  
  
  ### PCA for k=6
  fviz_cluster(kmeans_res[[6]], norm_betas,stand = F,axes = c(1,2),geom="point",
               ellipse.alpha = 0.1,ellipse.type = "norm",ggtheme = theme_classic())
  ### Hierarchical clustering on the random dataset
  fviz_dend(hclust(dist(norm_betas)), k = 6, k_colors = "jco",
            as.ggplot = TRUE, show_labels = FALSE)
  
  ### PCA for k=12
  fviz_cluster(kmeans_res[[12]], norm_betas,stand = F,axes = c(1,2),geom="point",
               ellipse.alpha = 0.1,ellipse.type = "norm",ggtheme = theme_classic()) 
  ### Hierarchical clustering on the random dataset
  fviz_dend(hclust(dist(norm_betas)), k = 12, k_colors = "jco",
            as.ggplot = TRUE, show_labels = FALSE)
 
  ### PCA for k=16
  fviz_cluster(kmeans_res[[16]], norm_betas,stand = F,axes = c(1,2),geom="point",
               ellipse.alpha = 0.1,ellipse.type = "norm",ggtheme = theme_classic()) 
  ### Hierarchical clustering on the random dataset
  fviz_dend(hclust(dist(norm_betas)), k = 16, k_colors = "jco",
            as.ggplot = TRUE, show_labels = FALSE)
  
  
  #####################################
  ##  Cluster Tendency with Hopkins  ##
  #####################################
  
  library(clustertend)
  
  # If the value of Hopkins statistic is close to zero, then we can reject 
  # the null hypothesis and conclude that the data set D is significantly 
  # a clusterable data.
 
  ###Compute Hopkins statistic for log fold changes
  set.seed(123)
  hopkins(norm_betas, n = nrow(norm_betas)-1)
  
  fviz_dist(dist(norm_betas), show_labels = FALSE)+
    labs(title = "TB infection over time (DM)")
  ###############################################################
  # for the number of clusters selected, get the k-means object##
  ###############################################################
  kmeans_out = kmeans(norm_betas,centers=nclus)
  print(names(kmeans_out))
  
  lfc_normalized <- norm_betas
  
  kmeans_clusters <- kmeans(lfc_normalized, centers=number_clusters,iter.max = 100, nstart = number_clusters)
  
  
  
  
  
  
  
  
  lfc_normalized$Nbclust <- best.partition.list$clust_res.Best.partition
  lfc_normalized$kmeans=kmeans_out$cluster
  lfc_normalized$ordered_genes=factor(rownames(lfc_normalized),levels=names(kmeans_out$cluster))
  # lfc_normalized$ordered_genes=factor(rownames(lfc_normalized),levels=rownames(best.partition.list))
  
  table=lfc_normalized[order(lfc_normalized$kmeans),]
  # table=lfc_normalized[order(lfc_normalized$Nbclust),]
  
  table$ordered_genes=factor(table$ordered_genes,levels=as.character(table$ordered_genes))
  tablem=melt(table[,c(1:4,6)],by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=levels(tablem$ordered_genes))
  tablem=tablem[order(tablem$ordered_genes),]
  
  colores_ref=brewer.pal(11,"RdBu")
  values_ref=seq(0,1,by=0.1)
  
  pl_tot <- ggplot(tablem) +
    geom_tile(aes(x=variable, y = ordered_genes, fill = value),colour=NA)+
    scale_fill_gradientn(colours = rev(colores_ref),values=values_ref)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, colour="black"),
          axis.text.y = element_blank())
  
  pl_tot
  cosa=table[,c(5,6)]
  
  # colnames(cosa)=c("beta_inf_48","ordered_genes")
  colnames(cosa)=c("Clusters","ordered_genes")
  
  tablem=melt(cosa,by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=levels(table$ordered_genes))
  tablem=tablem[order(tablem$value),]
  length(which(table$ordered_genes!=tablem$ordered_genes))
  
  number_clusters <- nclus
  #colores_ref=brewer.pal(number_clusters,"Set3")
  
  pl_kmeans <- ggplot(tablem) +
    geom_tile(aes(x=variable, y = ordered_genes, fill = factor(value)),colour=NA)+
    #scale_fill_manual(values = colores_ref)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, colour="black"),
          axis.text.y = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +  # Quita el espacio en blanco en el eje x
    scale_y_discrete(expand = c(0, 0))   # Quita el espacio en blanco en el eje y
  
  kmneans_plot=plot_grid(pl_kmeans,pl_tot,ncol=2,rel_widths=c(2,4))
  
  pdf(file.path(cluster_dir,paste0("001_kmneans_plot_k=",number_clusters,".pdf")),width=7,height=30)
  print(kmneans_plot)
  dev.off()
  
  
  number_clusters=12
  set.seed(123)
  kmeans_clusters <- kmeans(norm_betas[, 1:4], centers=number_clusters,iter.max = 100, nstart = number_clusters)
  
  norm_betas$kmeans=kmeans_clusters$cluster
  
  
  table=norm_betas[order(norm_betas$kmeans),]
  table$ordered_genes=factor(table$ordered_genes,levels=as.character(table$ordered_genes))
  tablem=melt(table[,c(1:4,6)],by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=levels(tablem$ordered_genes))
  tablem=tablem[order(tablem$ordered_genes),]
  
  colores_ref=brewer.pal(11,"RdBu")
  values_ref=seq(0,1,by=0.1)
  
  pl_tot <- ggplot(tablem) +
    geom_tile(aes(x=variable, y = ordered_genes, fill = value),colour=NA)+
    scale_fill_gradientn(colours = rev(colores_ref),values=values_ref)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, colour="black"),
          axis.text.y = element_blank())
  
  pl_tot
  cosa=table[,c(7,6)]
  
  # colnames(cosa)=c("beta_inf_48","ordered_genes")
  colnames(cosa)=c("Clusters","ordered_genes")
  
  tablem=melt(cosa,by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=levels(table$ordered_genes))
  tablem=tablem[order(tablem$value),]
  length(which(table$ordered_genes!=tablem$ordered_genes))
  
  #colores_ref=brewer.pal(number_clusters,"Set3")
  
  pl_kmeans <- ggplot(tablem) +
    geom_tile(aes(x=variable, y = ordered_genes, fill = factor(value)),colour=NA)+
    #scale_fill_manual(values = colores_ref)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, colour="black"),
          axis.text.y = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +  # Quita el espacio en blanco en el eje x
    scale_y_discrete(expand = c(0, 0))   # Quita el espacio en blanco en el eje y
  
  kmneans_plot=plot_grid(pl_kmeans,pl_tot,ncol=2,rel_widths=c(2,4))
  
  pdf(file.path(cluster_dir,paste0("003_kmneans_plot_k=",number_clusters,".pdf")),width=7,height=30)
  print(kmneans_plot)
  dev.off()
}

##################################################
##  2. Do clustering. Hierarchical clustering.  ##
##################################################

{
  
  ### Let´s try doing clustering using hierarchical clustering:
  
  clusters <- hclust(dist(norm_betas))
  clusters$order=rev(clusters$order)
  dend1 <- as.dendrogram(clusters)
  #dend1 <- color_branches(dend1, h=1.5)
  dend1 <- color_branches(dend1, k = 10)
  dend1 <- color_labels(dend1, k = 10)
  # plot(dend1)
  
  norm_betas$cluster_10=data.frame(cutree(dend1,k=10))[,1]
  norm_betas$ordered_genes=factor(rownames(norm_betas),levels=clusters$labels[clusters$order])
  norm_betas=norm_betas[order(norm_betas$ordered_genes),]
  clusters_as_they_appear=unique(norm_betas$cluster_10)
  
  ggd1=as.ggdend(dend1)
  arbol1=ggplot(ggd1, labels = FALSE,horiz = TRUE)
  arbol1
  
  pdf(file.path(cluster_dir,"001_arbol1_mod.pdf"),height=30,width=10)
  print(arbol1)
  dev.off()
  
  table=norm_betas
  tablem=melt(table[,c(1:4,6)],by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=(clusters$labels[clusters$order]))
  tablem=tablem[order(tablem$ordered_genes),]
  
  colores_ref=brewer.pal(11,"RdBu")
  # colores_ref=colorRampPalette(c("#0072B2", "#D55E00"))(100)
  values_ref=seq(0,1,by=0.1)
  
  pl_tot <- ggplot(tablem) +
    geom_tile(aes(x=variable, y = ordered_genes, fill = value),colour=NA,height=0.8)+
    scale_fill_gradientn(colours = rev(colores_ref),values=values_ref)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.background = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=14, colour="black"),
          axis.text.y = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +  # Quita el espacio en blanco en el eje x
    scale_y_discrete(expand = c(0, 0))   # Quita el espacio en blanco en el eje y
  pl_tot
  
  pdf(file.path(cluster_dir,"002_heatmap_exp_ok.pdf"),width=5,height=40)
  print(pl_tot)
  dev.off()
  
  ##################################################################################
  ##      This is problematic: if we run HCs cutting the tree at arbitrarily      ##
  ##    defined cutoffs (or selecting a certain number of clusters), we produce   ##
  ##      too many small clusters at some sections of the tree, or too large      ##
  ##     clusters at others... alternatively (which I did provisionally, but      ##
  ##     don´t quite like) one needs to define arbitrary cutoffs at different     ##
  ##     parts of the dendrogram to retrieve clusters that look ok, ít is too     ##
  ##                         arbitrary/hard to justify...                         ##
  ##################################################################################
  
}







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

### Heatmap function
heat_fun <- function(hclust_by_Rowv,gene_clust,col_dend,clust_name,column_title=character(0),show_column_names = F,column_order){
  
  ht_plt <- Heatmap(hclust_by_Rowv,
                    col = col_fun,
                    name="Expression Levels",
                    column_order = column_order,
                    show_column_names = show_column_names,
                    column_names_gp = gpar(fontsize = 6),
                    column_title = column_title,
                    column_title_side = "bottom",
                    border_gp = gpar(col = "black", lty = 2),
                    # heatmap_height = unit(6, "cm"),
                    # heatmap_width = unit(8, "cm"),
                    width=unit(3, "cm"),
                    # show_column_dend = T,
                    # column_dend_side = "top",
                    cluster_rows = color_branches(gene_clust),
                    # cluster_columns = color_branches(col_dend),
                    row_title = clust_name,
                    row_title_gp = gpar(fontize = 2),
                    row_dend_side = "right",
                    # row_dend_width = unit(2, "cm"),
                    show_row_names = F ,
                    show_row_dend = T,
                    heatmap_legend_param = list(title = "Expression level",
                                                title_position = "leftcenter-rot",
                                                labels_gp = gpar(font = 3),
                                                title_gp = gpar( fontsize = 8)))
  return(ht_plt)
}





