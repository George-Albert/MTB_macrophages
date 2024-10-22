
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
  library(puma)
  library(umap)
  library(circlize)
  library(ComplexHeatmap)
  
  
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
  
  grep("beta",colnames(tab))
  colnames(tab)[c(1,  6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61)]
  grep("hit",colnames(tab))
  colnames(tab)[c(5,10,15,20,25,30,35,40,45,50,55,60,65)]
  grep("SE",colnames(tab))
  colnames(tab)[c(3,8,13,18,23,28,33,38,43,48,53,58,63)]
  
  # Select these betas
  # "beta_inf_2"         "beta_inf_20"        "beta_inf_48"        "beta_inf_72" 
  
  tab$hit_any=apply(tab[,c(5,10,15,20)],1,function(x){as.numeric(sum(x)>0)})
  tab_hit_labels=tab[which(tab$hit_any==1),c(5,10,15,20)]
  tab=tab[which(tab$hit_any==1),c(1,6,11,16)]
  
  ### Normalize the betas of each gene time-wise.
  ttab=t(tab)
  norm_ttab=sapply(1:ncol(ttab),function(x){ttab[,x]/max(abs(ttab[,x]))})
  tab_norm=t(norm_ttab)
  rownames(tab_norm)=rownames(tab)
  
  ### Normalize the SE of each gene time-wise.
 
  tab_se <- results[rownames(tab_norm), c(3,8,13,18)]
  ttab_se=t(tab_se)
  norm_ttab_se = sapply(1:ncol(ttab_se), function(x) { ttab_se[,x] / max(abs(ttab[,x])) })
  norm_se <- t(norm_ttab_se)
  rownames(norm_se)=rownames(tab_se)
  length(which(rownames(tab_norm)!=rownames(norm_se)))
  
  # Let us rename these objects to store them more properly/safely
  betas=tab
  norm_betas=data.frame(tab_norm)
  se <- data.frame(norm_se)

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
          axis.text.y = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +  # Quita el espacio en blanco en el eje x
    scale_y_discrete(expand = c(0, 0))   # Quita el espacio en blanco en el eje y
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
  number_clusters=8
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


###############################
##  3. Clustering with PUMA  ##
###############################


#########################
### Declare functions ### 
#########################

### variance formula Varianza = Σ((xi - x̄)^2) / (n - 1)
var_fun <- function(x,se_included){
  
  if (se_included==T) {
    x <- x^2
    x_norm2=rowSums(x)
    print(length(x_norm2))
    out <- sum(x_norm2*(se^2))/length(x_norm2)
  }else{  x <- x^2
  x_norm2=rowSums(x)
  print(length(x_norm2))
  out <- sum(x_norm2)/length(x_norm2)
  }
  
}

### Explained variance function
exp_var_fun_clust <- function(LogFC, se,se_included = T,max_n_clusters) {
  
  e <- LogFC
  
  ### Make sure to always get the same result
  set.seed(123)
  # cl<-pumaClust(e = e,se = se ,clusters=9)
  
  ### Maximun number of clusters
  max_n_clusters <- max_n_clusters
  
  ### create 5 lists where Im going to save the results
  
  pumacl <- list()
  var_list <- list()
  clusters_centrados <- list()
  bic_list <- list()
  bic_out_list <- list()
  optF_list <- list()
  likelipergene_list <- list()
  
  
  # pb <- txtProgressBar(min=0,max=100)
  for (num_clusters in seq_len(max_n_clusters)) {
    
    ### Run pumaClustii and keep the results in the pumacl list
    # cl<-pumaClust(e = e,se = se, clusters=9)
    cl<-pumaClust(e = e,se = se, clusters=num_clusters)
    pumacl[[num_clusters]] <- cl
    
    print(paste0("The number of clusters is ",num_clusters))
    
    ### Outputs from pumaClustii
    clusters <- data.frame(clusters=cl$cluster)
    n_clusters <- num_clusters
    resumen <- summary(as.factor(cl$cluster))
    numero_de_clusters <- unique(cl$cluster)
    numero_de_clusters <- numero_de_clusters[order(numero_de_clusters)]
    
    centers <- data.frame((cl$centers))
    centers <- centers[numero_de_clusters,]
    colnames(centers) <- colnames(e)[1:num_contras]
    rownames(centers) <- paste0("cluster_",numero_de_clusters)
    variance_of_center <- data.frame(cl$centersigs)
    # colnames(variance_of_center) <- colnames(e)
    likelipergene <- data.frame(cl$likelipergene)
    colnames(likelipergene) <- paste0("cluster_",numero_de_clusters)
    
    bic_out <- cl$bic
    bic_out_list[[num_clusters]] <- bic_out
    
    ### Assign clusters and colors to the data we are using (expression level of each gene) ###
    colors <- rainbow(num_clusters)
    clusters$colors <- NA
    
    ### Obtain unique clusters numbers
    unique_numbers <- unique(clusters$clusters)
    
    ### Match the colors with the cluster 
    clusters$colors <- colors[match(clusters$clusters, unique_numbers)]
    clusters <- clusters[order(clusters$clusters),]
    
    ### Center the raw data at the origin
    LogFC <- e[rownames(clusters),]
    LogFC$clusters <- clusters$clusters
    length(which(rownames(LogFC)!=rownames(clusters)))
    
    ###Obtain the medoid and rest it from the data
    medoid <- colMeans(LogFC[,1:num_contras])
    centerd_obs <- t(apply(LogFC[,1:num_contras], 1, function(x) x - medoid))  
    centerd_obs <- as.data.frame(centerd_obs)
    
    length(which(rownames(centerd_obs)!=rownames(clusters)))
    centerd_obs$clusters <- clusters$clusters
    
    ### Separate by clusters ###
    print("Separating clusters...")
    
    separate_clusters_list <- split(centerd_obs, centerd_obs$clusters)
    names(separate_clusters_list) <- paste0("clusters_",names(separate_clusters_list))
    ### remove the column with the cluster number
    # separate_clusters_list <- lapply(separate_clusters_list, function(x) x[,1:(ncol(x)-1)])
    
    ### Center at the origin the clusters ###
    print(paste0("the number of clusters is ",num_clusters))
    center_by_clust <- centers[1:num_clusters,]
    
    vec_clust <- seq_along(separate_clusters_list)
    for(i in vec_clust){
      # Define the vector and the number of times to repeat it
      rest_vector <- as.numeric(center_by_clust[i,])
      # rest_vector <- unname(rest_vector)
      n <- nrow(separate_clusters_list[[i]])
      
      # Create the data frame with the repeated vector
      rest_dataframe <- data.frame(matrix(rep(rest_vector, n), nrow = n, byrow = TRUE))
      
      # Display the resulting data frame
      selected_cluster <- separate_clusters_list[[i]][,1:(ncol(separate_clusters_list[[i]])-1)]
      clusters_centrados[[i]]<- selected_cluster-rest_dataframe
    }
    
    names(clusters_centrados) <- names(separate_clusters_list)
    clusters_centrados_df <- do.call(rbind,clusters_centrados)
    
    # centered_clusters_df <- do.call(rbind,centered_clusters_list)
    print("Computing residual variance")
    ### Compute residuals variance ####
    res_var <- var_fun(clusters_centrados_df,se_included=se_included)
    
    print("Computing total variance")
    ### Compute total variance ####
    var_total <- var_fun(centerd_obs,se_included=se_included)
    print("Computing explained variance")
    explained_var <- 1-(res_var/var_total)
    print(paste0("The Explained Variance is ",explained_var))
    print("Saving explained total variances values to a list")
    
    var_list[[num_clusters]] <- explained_var
    names(var_list[num_clusters]) <- paste0("K = ", num_clusters)
    # setTxtProgressBar(pb, num_clusters)
    #### Compute BIC=n*ln(RSS/n)+k*ln(n) ####
    print("Computing BIC score")
    
    likelihood_per_gene <- apply(likelipergene, 2, sum)
    likelihood <- sum(log(likelihood_per_gene))
    n <- nrow(clusters_centrados_df)
    # n <- nrow(likelipergene)
    ssr <- res_var
    k <- num_clusters
    # k <- sum(likelipergene > 0)
    # bic_score <- n*log(ssr) + k*log(n)
    bic_score <- -2 * likelihood + k * log(n)
    
    bic_list[[num_clusters]] <- bic_score
    print(paste0("The BIC is ",bic_score))
    
  }
  return(list(bic_list,var_list,clusters_centrados_df,pumacl,likelipergene_list,bic_out_list))
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

### Create the vertical heatmaps
create_heatmaps <- function(n) {
  heatmaps <- list()
  # n is num_clusters
  for (i in 1:n) {
    
    col_dend <- col_den_list[[i]]
    cluster_order <- rownames(matrix_cluster[[i]])
    gene_clust <- h_clust_list[[i]]
    k <- paste0("k=",i)
    columns_name <- colnames(matrix_cluster[[num_contras]])
    # Create a heatmap
    if (i==length(matrix_cluster)) {
      
      ht_plt=heat_fun(matrix_cluster[[i]],gene_clust,col_dend,clust_name=k,column_title = "Conditions",show_column_names=T,
                      column_order=columns_name)
    }else{
      ht_plt <- heat_fun(matrix_cluster[[i]],gene_clust,col_dend,clust_name=k,column_order=columns_name)
    }
    
    # Agregate a heatmap to the list
    heatmaps[[i]] <- ht_plt
  }
  
  # Combinar los heatmaps verticalmente
  combined_heatmap <- Reduce(`%v%`, heatmaps)
  
  return(combined_heatmap)
}
##########################

### Time-wise Normalized table. We are going to use it in PUMA

puma_lfc <- norm_betas[,1:4]
# tab1 <- results[rownames(tab), c(3,8,13,18)]
puma_se <- se

set.seed(123)
library(mclust)
gmm_model <- Mclust(puma_lfc)

# Show the optimal number of clusters based on BIC
optimal_clusters <- gmm_model$G
cat("Número óptimo de clusters basado en BIC:", optimal_clusters, "\n")

# Open a PDF device
pdf(file.path(cluster_dir,"004_GMM_according_BIC.pdf"), width = 8, height = 6)
plot(gmm_model, what = "BIC")
dev.off()

# aic <- gmm_model$loglik - length(gmm_model$parameters$mean) # AIC aproximado
# cat("Valor AIC:", aic, "\n")


k_values <- 1:20  # Rango de clusters que quieres probar
bic_values <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  # RUn PumaClust to define the numeber of clusters based on BIC
  puma_result <- pumaClust(puma_lfc,se=puma_se,subset = rownames(puma_lfc), clusters = k_values[i])
  
  # Extract BIC value
  bic_values[i] <- puma_result$bic
}

pdf(file.path(cluster_dir,"005_PUMA_BIC.pdf"), width = 8, height = 6)
plot(k_values, bic_values, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters",
     ylab = "Bayesian Information Criterion (BIC)",
     main = "Selección de Número de Clusters basado en BIC")
dev.off()

# Identify the optimal number of clusters (minimum BIC value found)
optimal_k <- k_values[which.min(bic_values)]
cat("El número óptimo de clusters es:", optimal_k)

# Another way using my function
num_contras <- ncol(betas)
max_n_clusters <- 25
list_plt <- exp_var_fun_clust(puma_lfc, puma_se, se_included = F,max_n_clusters=max_n_clusters)

{ 
  bic_list <- list_plt[[1]]
  var_list <- list_plt[[2]]
  var_list[is.null(var_list)]
  centered_clusters_list <- list_plt[[3]]
  pumacl_list <-list_plt[[4]]
  likelipergene_list <- list_plt[[5]]
  bic_out_list <- list_plt[[6]]
  ### plot variance and BIC
  var_df <- data.frame(t(data.frame(var_list)))
  rownames(var_df) <- c(1:max_n_clusters)
  num_clusters <- which.min(var_df[,1])
  
  x <- as.numeric(rownames(var_df))
  y <- var_df[,1]
  
  var_plt <- ggplot(data=data.frame(var_df), aes(x=x, y=100*y)) +
    geom_line(linetype = "dashed")+
    geom_point(size=3)+
    # xlim(0,max_n_clusters)+
    # ylim(-1,100)+
    xlab("numbers of clusters")+
    ylab("Explained Variance")+
    theme_classic()
  
  var_plt
  
  pdf(file.path(cluster_dir,"006_Explained_variance"),width=9,height=6)
  print(var_plt)
  dev.off()
  
  bic_df <- data.frame(t(data.frame(bic_out_list)))
  rownames(bic_df) <- c(1:max_n_clusters)
  num_clusters <- which.min(bic_df[,1])
  
  x1 <- as.numeric(rownames(bic_df))
  y1 <- bic_df[,1]
  # bic_df[which(bic_df$t.data.frame.bic_list..==0),] <- 1250
  bic_plt <- ggplot(data=data.frame(bic_df), aes(x=x1, y=y1)) +
    geom_line(linetype = "dashed")+
    geom_point(size=3)+
    xlim(0,max_n_clusters)+
    # ylim(-1,100)+
    xlab("numbers of clusters")+
    ylab("BIC score")+
    theme_classic()
  
  bic_plt
  
  cat("The minimum BIC values is at cluster number: ",which.min(bic_df[[1]]))
  
  pdf(file.path(cluster_dir,"007_Explained_variance_BIC"),width=9,height=6)
  print(bic_plt)
  dev.off()
  
}

################################
##  Heatmap from PUMA output  ##
################################

set.seed(123)

num_contras <- ncol(betas)
e <- puma_lfc
se <-puma_se

# We defined the optimal number of clusters based on the explained Variance and BIC
k=12
# k=8

cl <- pumaClust(e=e,se=se,subset = rownames(puma_lfc), clusters = k)

### Outputs from pumaClustii
clusters <- data.frame(clusters=cl$cluster)
n_clusters <- max(clusters)
resumen <- summary(as.factor(cl$cluster))
# Para k = 8
# 1    2    3    4    5    6    7    8 
# 164  349 1266 1288  637  400  165  622 

# Para k = 12
#    1    2    3    4    5    6    7    8    9   10   11   12 
# 1055  136  155    2  731  135  495  164  383 1079  550    6 

numero_de_clusters <- unique(cl$cluster)
# numero_de_clusters <- numero_de_clusters[order(numero_de_clusters)]

e$clusters <- clusters$clusters 
e <- e[order(e$clusters),]

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white","red"))

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white","red"))
columns_name <- colnames(e[,1:num_contras])
hclust_matrix <- as.matrix(e[,1:num_contras])
column_title = "Conditions"
clust_name <- "Genes"

# Definir anotación de fila para los clusters
ha <- rowAnnotation(Cluster = as.factor(e$clusters))
e$clusters <- as.factor(e$clusters)

ht_plt <- Heatmap(hclust_matrix,
                  na_col = "grey2",
                  col = col_fun,
                  split = e$clusters,
                  # row_split = e$clusters,
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
                  width=unit(2, "cm"),
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
filename <- file.path(cluster_dir,paste0("008_Heatmap_clustering_2_normalized_k=",k))

pdf(file =filename ,width=size[1]+1,height=size[2]+2)
draw(ht_plt,ht_gap=unit(0.5,"mm"))
dev.off()


# dim(e[which(e[,4]>0 & e$clusters==1),])






