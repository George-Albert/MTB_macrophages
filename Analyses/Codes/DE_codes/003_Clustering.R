
#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyr)
  library(gridExtra)
  library(dplyr)
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
}
{
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
    
    kmeans_res[[number_clusters]] <- kmeans(norm_betas[, 1:4], centers=number_clusters,
                                            iter.max = 100, 
                                            nstart = number_clusters)
    cluster_list[[number_clusters]] <-kmeans_res[[number_clusters]]$cluster
    
  }
  
  wss <- numeric(20)
  wss[1] <- (nrow(norm_betas)-1)*sum(apply(norm_betas,2,var))
  
  # Compute the within sum of squares (WSS) for each number of clusters
  for (i in 2:20) {
    wss[i] <- sum(kmeans(norm_betas, centers = i,iter.max = 100, 
                         nstart = number_clusters)$withinss)
    
  }
  
  # Realizar el gráfico después de completar el bucle
  plot(1:20, wss, type = "b", xlab = "Number of Clusters",
       ylab = "Within groups sum of squares")
  
 ### Function to calculate the BIC and the AIC for clustering algorithms
   
  kmeansAIC_BIC = function(fit){
    
    m = ncol(fit$centers)   # Number of conditions or samples
    n = length(fit$cluster) # Number of genes or rows
    k = nrow(fit$centers)   # Number of centers or clusters
    D = fit$tot.withinss    # Total within sum of squares
    return(data.frame(AIC = D + 2*m*k,
                      BIC = D + log(n)*m*k))
  }
  
   ### Define the list to capture AIC and BIC outputs
   ic_values          <- list()
   ic_values_hckmeans <- list()
   ### Define the list to capture kmeans an hckmeans outputs
   
   kmeans_model   <- list()
   hckmeans_model <- list()  
   ### Define the max number of clusters yopu want to inspect
   k_vec <- seq(40)
   
   for (k in k_vec) {
     
     kmeans_model[[k]] <- kmeans(norm_betas, 
                                 centers = k, 
                                 nstart = k,
                                 iter.max = 100)
     hckmeans_model[[k]] <- hkmeans(norm_betas,
                                    k=k,
                                    hc.metric = "euclidean",
                                    hc.method = "ward.D2",
                                    iter.max = 100,
                                    km.algorithm = "Hartigan-Wong")
     
     ic_values[[k]]          <- kmeansAIC_BIC(kmeans_model[[k]])
     ic_values_hckmeans[[k]] <- kmeansAIC_BIC(hckmeans_model[[k]])
   }
   
   
   ### Extract the AIC and BIC values form the list
   AIC_values <- sapply(ic_values, function(x) x[[1]])
   BIC_values <- sapply(ic_values, function(x) x[[2]])
   AIC_values_hckmeans <- sapply(ic_values_hckmeans, function(x) x[[1]])
   BIC_values_hckmeans <- sapply(ic_values_hckmeans, function(x) x[[2]])
  
   ### Create the df with all the information criterion from the two models
   ic_df <- data.frame(k.clusters=k_vec, AIC.Values=AIC_values,BIC.Values=BIC_values,
                       AIC.Values.hckmeans=AIC_values_hckmeans,BIC.Values.hckmeans=BIC_values_hckmeans)
   
   ### Plot and save the graphs
   
   #ic_df         <- data frame with the number of clusters and the AIC and BIC values of each model
   #col_to_select <- string of the colname with the respective clustering model
   #ylab          <- y axis label as string
   #color         <- color of the points
   #color_opt     <- color orf the point representing the optimal number of clusters
   
   info_crit_plt <-function(ic_df,col_to_select,ylab ="IC",color ="red",color_opt="blue") {
     
     min_ic <- min(ic_df[[col_to_select]])
     cluster_selected <- ic_df[which(ic_df[col_to_select]==min_ic),1]
     cat("the optimal number of clusters is:",cluster_selected)
     
     IC_plt <- ggplot(ic_df)+
       geom_point(aes(x=k.clusters,y=.data[[col_to_select]]),size=2.5,color=color)+
       geom_point(data = subset(ic_df, ic_df[[col_to_select]] == min_ic),
                  aes(x = k.clusters, y = !!sym(col_to_select)), color = color_opt, size = 3.5) + # Punto mínimo en azul
       xlab("Number of Clusters")+
       ylab(ylab)+
       theme_classic()+
       theme(axis.text.y   = element_text(size=18),
             axis.text.x   = element_text(size=18),
             axis.title.y  = element_text(size=18),
             axis.title.x  = element_text(size=18),
             panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"))
     
     pdf(file.path(cluster_dir,paste0("001_",ylab,"_plt.pdf")),width = 12,height = 10)
     print(IC_plt)
     dev.off()
     
     list_to_return <- list(IC_plt,cluster_selected)
     return(list_to_return)
   }
   
   AIC_plt<- info_crit_plt(ic_df = ic_df,col_to_select = "AIC.Values",ylab = "AIC" )
   cluster_selected_AIC <- AIC_plt[[2]]
   # the optimal number of clusters is: 24
   
   BIC_plt <- info_crit_plt(ic_df = ic_df,col_to_select = "BIC.Values",ylab = "BIC",color = "blue",color_opt = "red")
   cluster_selected_BIC <- BIC_plt[[2]]
   # the optimal number of clusters is: 9
   
   AIC_hckmeans_plt <- info_crit_plt(ic_df = ic_df,col_to_select = "AIC.Values.hckmeans",ylab = "AIC_hckmeans" )
   cluster_selected_AIC_hckmeans <- AIC_hckmeans_plt[[2]]
   # the optimal number of clusters is: 22
   
   BIC_hckmeans_plt <- info_crit_plt(ic_df = ic_df,col_to_select = "BIC.Values.hckmeans",ylab = "BIC_hckmeans",color = "blue",color_opt = "red" )
   cluster_selected_BIC_hckmeans <- BIC_hckmeans_plt[[2]]
   # the optimal number of clusters is: 11

   
   ### check the order of the gene names in each data
   length(which(rownames(norm_betas) != names(kmeans_model[[24]]$cluster)))
   length(which(rownames(norm_betas) != names(hckmeans_model[[22]]$cluster)))
  
   names_of_df_col <- c(paste0("kmeans_cluster_k=",cluster_selected_BIC),
                        paste0("kmeans_cluster_k=",cluster_selected_AIC),
                        paste0("hckmeans_cluster_k=",cluster_selected_BIC_hckmeans),
                        paste0("hckmeans_cluster_k=",cluster_selected_AIC_hckmeans))
   
   }  

#################################
## Mclust clustering algorithm ##
#################################
mc <- Mclust(norm_betas) # Model-based-clustering
summary(mc) # Print a summary

mc$modelName # Optimal selected model ==> "VEE"
mc$G # Optimal number of cluster => 8
head(mc$z, 30) # Probality to belong to a given cluster
head(mc$classification, 30) # Cluster assignement of each observation

Mclust_name <- paste0("Mclust_k=",mc$G)
cluster_out <- data.frame(norm_betas,
                          kmeans   = kmeans_model[[cluster_selected_BIC]]$cluster,
                          kmeans2  = kmeans_model[[cluster_selected_AIC]]$cluster,
                          hckmeans = hckmeans_model[[cluster_selected_BIC_hckmeans]]$cluster,
                          hckmeans2= hckmeans_model[[cluster_selected_AIC_hckmeans]]$cluster,
                          Mclust.  = mc$classification)

colnames(cluster_out) <-c(colnames(norm_betas), names_of_df_col,Mclust_name)

# BIC values used for choosing the number of clusters
BIC_plt <- fviz_mclust(mc, "BIC", palette = "jco")

pdf(file.path(cluster_dir,paste0("001_BIC_Mclust.pdf")),width=12,height=10)
print(BIC_plt)
dev.off()

cluster_input_dir <- file.path(input_dir,"003_Clustering_tabs")
dir.create(cluster_input_dir,recursive = T,showWarnings = F)

saveRDS(kmeans_model, file = file.path(cluster_input_dir,"001_kmeans_clustering.RDS")) 
saveRDS(hckmeans_model, file = file.path(cluster_input_dir,"001_hckmeans_clustering.RDS")) 

write.table(cluster_out,file.path(cluster_input_dir,"lfc_with_clustering.txt"))

###############################################################
# for the number of clusters selected, get the k-means object##
###############################################################

{

cluster_out <- read.table(file.path(cluster_input_dir,"lfc_with_clustering.txt"))

# col_fun = colorRamp2(c(-1, 0, 1), c("#FFA373", "white","#50486D"))
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
# Colors using colorRampPalette
n_colors <- 25  # Total de colores que quieres
palette_continuous <- colorRampPalette(cluster_colors)(n_colors)

# Now we use the selected colors
cluster_colors <- palette_continuous
names(cluster_colors) <- 1:length(cluster_colors)

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
  filename <- file.path(cluster_dir,paste0("002_Heatmap_normalized_",colnames(table)[5],".pdf"))
  
  pdf(file =filename ,width=size[1]+1,height=size[2]+2)
  draw(ht_plt,ht_gap=unit(0.5,"mm"))
  dev.off()
  
  }
}


######################################
##  Compute the trends per cluster  ##
######################################
{
# Compute the mean and the sd by cluster
result_mean_list <- list() # list to save 
result_sd        <- list() # list to save the sd 
result_mean      <- list() # list to save the mean

df <- cluster_out 
vec_of_clust_columns <- c(5:ncol(df))

for (i in vec_of_clust_columns){
  
  folder_to_save <- colnames(df)[i]
  df_por_clusters <- split(df, df[,i])
  names(df_por_clusters) <- paste0("k=",names(df_por_clusters))
  
  #j is going to from k=0 up the max k clusters of the list
  plots_list       <- list() # list to save plt to create a grid
  
  for (j in seq(df_por_clusters)) {
    
    result_mean_list[[j]] <- df_por_clusters[[j]][1:4] %>%
      summarize(
        Mean_beta_inf_2 = mean(beta_inf_2),
        Mean_beta_inf_20 = mean(beta_inf_20),
        Mean_beta_inf_48 = mean(beta_inf_48),
        Mean_beta_inf_72 = mean(beta_inf_72),
        sd_beta_inf_2 = sd(beta_inf_2),
        sd_beta_inf_20 = sd(beta_inf_20),
        sd_beta_inf_48 = sd(beta_inf_48),
        sd_beta_inf_72 = sd(beta_inf_72)
      )
    
    result_sd[[j]] <- result_mean_list[[j]][5:8]
    result_mean[[j]] <- result_mean_list[[j]][1:4]
    
    data_mean <- pivot_longer(result_mean[[j]],cols = c(Mean_beta_inf_2,Mean_beta_inf_20,Mean_beta_inf_48,Mean_beta_inf_72),
                              names_to = "Conditions",values_to = "Mean")
    data_mean$Conditions <- factor(data_mean$Conditions,levels = c("Mean_beta_inf_2","Mean_beta_inf_20","Mean_beta_inf_48","Mean_beta_inf_72") )
    
    data_sd <- pivot_longer(result_sd[[j]],cols = c(sd_beta_inf_2,sd_beta_inf_20,sd_beta_inf_48,sd_beta_inf_72),
                            names_to = "Conditions",values_to = "sd")
    data_sd$Conditions <- factor(data_sd$Conditions,levels = c("sd_beta_inf_2","sd_beta_inf_20","sd_beta_inf_48","sd_beta_inf_72"))
    
    dataframe <- data.frame(data_mean,sd=data_sd$sd)
    dataframe$x_num <- 1:nrow(dataframe)
    
    myarrow=arrow(angle = 15, ends = "last", type = "closed")
    col <- c("red", "green", "blue","cyan")
    # Define the line color and SD shading
    line_color <- "black" # Change this value to the desired line color
    ribbon_color <- "grey30" # Change this value to the desired shading color
    
    mean_plt_fun <- ggplot(data = dataframe, aes(x = x_num, y = Mean)) +
      geom_line(aes(group = 1), color = line_color, arrow = myarrow, lwd = 2) +
      geom_point(color = col, size = 4) +
      geom_ribbon(aes(ymin = Mean - sd, ymax = Mean + sd), fill = ribbon_color, alpha = 0.2) +
      scale_x_continuous(breaks = dataframe$x_num, labels = dataframe$Conditions) +
      ylab("Mean") +
      xlab("Conditions") +
      ggtitle(paste0("k=",j)) +
      theme_classic() +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 90),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", 
                                    fill = NA, 
                                    linewidth = 1, 
                                    linetype = "solid"))
    
    # print(mean_plt_fun)
    
    # Save the plot in the list
    plots_list[[paste0("k=", j)]] <- mean_plt_fun
    
    trend_dir <- file.path(cluster_dir,"Trends_per_cluster",folder_to_save)
    dir.create(trend_dir,showWarnings = F,recursive = T)
    
    pdf(file.path(trend_dir,paste0("k=",j,"_plt.pdf")),width = 10,height = 10)
    print(mean_plt_fun)
    dev.off()
    
    # Calculate the number of plots
    num_plots <- length(plots_list)
    
    # Calculate the number of columns and rows to automatically adjust the layout
    ncol_plots <- ceiling(sqrt(num_plots))  # Square root approximation for plot distribution
    nrow_plots <- ceiling(num_plots / ncol_plots)  # Adjust based on the number of columns
    
    pdf(file.path(trend_dir,"grid_plt.pdf"),width = 40,height = 38)
    # Use do.call with grid.arrange to display the plots with the appropriate number of rows and columns
    do.call(grid.arrange, c(plots_list, ncol = ncol_plots, nrow = nrow_plots))
    dev.off()
  }
  
}
}










































