
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
  library(dendextend)
  library(reshape2)
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

#############################################################
##  1. Do clustering. Version 1: Hierarchical clustering.  ##
#############################################################

{
  ## Build table with the infection logFCs of genes that respond to infection
  ##  at at least one timepoint.
  
  tab=results
  tab$hit_any=apply(tab[,c(4,8,12,16)],1,function(x){as.numeric(sum(x)>0)})
  tab_hit_labels=tab[which(tab$hit_any==1),c(4,8,12,16)]
  tab=tab[which(tab$hit_any==1),c(1,5,9,13)]
  
  ### Normalize the betas of each gene time-wise.
  ttab=t(tab)
  norm_ttab=sapply(1:ncol(ttab),function(x){ttab[,x]/max(abs(ttab[,x]))})
  tab_norm=t(norm_ttab)
  rownames(tab_norm)=rownames(tab)
  
  # Let us rename these objects to store them more properly/safely
  betas=tab
  norm_betas=tab_norm
  
  tab=data.frame(norm_betas)
  
  ### Let´s try doing clustering using hierarchical clustering:
  
  clusters <- hclust(dist(tab))
  clusters$order=rev(clusters$order)
  dend1 <- as.dendrogram(clusters)
  #dend1 <- color_branches(dend1, h=1.5)
  dend1 <- color_branches(dend1, k = 10)
  dend1 <- color_labels(dend1, k = 10)
  #plot(dend1)
  
  tab$cluster_10=data.frame(cutree(dend1,k=10))[,1]
  tab$ordered_genes=factor(rownames(tab),levels=clusters$labels[clusters$order])
  tab=tab[order(tab$ordered_genes),]
  clusters_as_they_appear=unique(tab$cluster_10)
  
  
  ggd1=as.ggdend(dend1)
  arbol1=ggplot(ggd1, labels = FALSE,horiz = TRUE)
  arbol1
  
  pdf(file.path(cluster_dir,"001_arbol1_mod.pdf"),height=30,width=10)
  print(arbol1)
  dev.off()
  
  table=tab
  tablem=melt(table[,c(1:4,6)],by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=(clusters$labels[clusters$order]))
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
  
  pdf(file.path(cluster_dir,"002_heatmap_exp_ok.pdf"),width=5,height=30)
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

############################################
##  2. Do clustering. Version 2: k-means  ##
############################################

{
  number_clusters=18
  set.seed(123)
  kmeans_clusters <- kmeans(tab[, 1:4], centers=number_clusters,iter.max = 100, nstart = number_clusters)
  
  tab$kmeans=kmeans_clusters$cluster
  
  
  table=tab[order(tab$kmeans),]
  table$ordered_genes=factor(table$ordered_genes,levels=as.character(table$ordered_genes))
  tablem=melt(table[,c(1:4,6)],by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=levels(table$ordered_genes))
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
  
  colnames(cosa)=c("beta_inf_48","ordered_genes")
  
  tablem=melt(cosa,by=c("ordered_genes"))
  
  tablem$ordered_genes=factor(tablem$ordered_genes,levels=levels(table$ordered_genes))
  tablem=tablem[order(tablem$ordered_genes),]
  
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
          axis.text.y = element_blank())
  
  kmneans_plot=plot_grid(pl_kmeans,pl_tot,ncol=2,rel_widths=c(2,4))
  
  pdf(file.path(cluster_dir,"003_kmneans_plot.pdf"),width=7,height=30)
  print(kmneans_plot)
  dev.off()
}
