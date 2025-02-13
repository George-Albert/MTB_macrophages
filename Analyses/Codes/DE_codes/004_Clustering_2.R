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
cluster_input_dir <- file.path(input_dir,"003_Clustering_tabs")
cluster_dir <- file.path(output_dir,"004_Clustering_plots")

### load data
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
cols=read.table(file.path(input_dir,"002_Processed","ready_for_DE","metadata.txt"))
cols_whole=read.table(file.path(input_dir,"002_Processed","whole","metadata_whole.txt"))
results <- read.table(file.path(output_dir,"003_DE_def","resultados.txt"))

length(which(colnames(reads)!=rownames(cols)))
length(which(paste0(cols$Individual,"_",cols$Setup)!=rownames(cols)))
length(which(rownames(reads)!= rownames(results)))

cluster_out1 <- read.table(file.path(cluster_input_dir,"lfc_data_with_clusters_information.txt"))
### We selected the k=16 cluster results
cluster_out1 <- cluster_out1[,c(1:4,7)]# Select the k=16

### load hckmeans results
hckmeans_model <- readRDS(file.path(cluster_input_dir,"001_hckmeans_clustering.RDS"))

### We are going to select the k=16 clustering result
hckmeans_model <- hckmeans_model[[16]]
table(hckmeans_model$cluster)
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
# 729 622 598 158 509 460 500  65 275 102 405  78 202  75  84  29

### Load the BIC and the AIC values
AIC_hckmeans <- read.table(file.path(cluster_input_dir,"AIC_hckmeans.txt"))
BIC_hckmeans <- read.table(file.path(cluster_input_dir,"BIC_hckmeans.txt"))

ic_df <- cbind(AIC_hckmeans,BIC_hckmeans)
ic_df$k_clusters <- as.numeric(seq(ic_df$x))
colnames(ic_df) <- c("AIC","BIC","k_clusters")

#################
##  Functions  ##
#################

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

#########################
##  AIC and BIC plots  ##
#########################
{ ### plot the BIC and AIC values togheter 

  min_aic     <- min(ic_df$AIC)
  min_bic     <- min(ic_df$BIC)
  col_min_bic <- "#FFA373"
  col_bic     <- "#50486D"
  col_min_aic <- "blue"
  col_aic     <- "red"
  
  IC_plt <- ggplot(ic_df[1:nrow(ic_df),])+
    geom_point(aes(x=k_clusters,y=AIC,color="AIC"),size=3)+
    geom_point(data = subset(ic_df, AIC == min_aic),
               aes(x = k_clusters, y = AIC, color = "Min AIC"),size = 5) + 
    geom_point(aes(x=k_clusters,y=BIC,color="BIC"),size=3)+
    geom_point(data = subset(ic_df, BIC == min_bic),
               aes(x = k_clusters, y = BIC, color = "Min BIC"),size = 5)+
    scale_color_manual(
      name = "Information Criteria",
      values = c(
        "AIC" = col_aic,
        "BIC" = col_bic,
        "Min AIC" = col_min_aic,
        "Min BIC" = col_min_bic
      ),
      labels = c("AIC", "BIC", "Min AIC", "Min BIC")
    ) +
    geom_vline(xintercept = 16, linetype = "dashed", color = "black", 
               linewidth = 1.2) +
    xlab("Number of Clusters")+
    ylab("AIC and BIC")+
    theme_classic()+
    theme(axis.text.y   = element_text(size=18),
          axis.text.x   = element_text(size=18),
          axis.title.y  = element_text(size=18),
          axis.title.x  = element_text(size=18),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black",
                                      fill=NA, linewidth=1,
                                      linetype="solid"))
  IC_plt
  
  pdf(file.path(cluster_dir,paste0("001_AIC_and_BIC_k_16_plt.pdf")),width = 15,
      height = 10)
  print(IC_plt)
  dev.off()
}


####################
##  Heatmap k=14  ##
####################

{
  
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
  
  lfc_normalized <- cluster_out1
  
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
  
  ### Compute the size of the heatmap
  size <- calc_ht_size(ht_plt)
  
  ### Save the data
  filename <- file.path(cluster_dir,paste0("003_Heatmap_hckmeans_",colnames(table)[5],".pdf"))
  
  pdf(file =filename ,width=size[1]+1,height=size[2]+2)
  draw(ht_plt,ht_gap=unit(0.5,"mm"))
  dev.off()
  }
###################
##  Trend plots  ##
###################

result_mean_to_save <- readRDS(file.path(cluster_input_dir,"002_df_Mean_sd_by_cluster_k_16.RDS"))

  # Compute the mean and the sd by cluster

  {
  df <- cluster_out1
  folder_to_save <- paste0("hckmeans_",colnames(df)[5])
  df_por_clusters <- split(df, df[,5])
  names(df_por_clusters) <- paste0("k=",names(df_por_clusters))
  
  # list to save plt to create a grid
  plots_list <- list() 
  
  for (j in seq(df_por_clusters)) {
    
    dataframe <- result_mean_to_save[[j]]

    myarrow=arrow(angle = 15, ends = "last", type = "closed")
    col <- c("black","red", "green", "blue","cyan")
    # Define the line color and SD shading
    line_color   <- "black" # Change this value to the desired line color
    ribbon_color <- "grey30" # Change this value to the desired shading color
    data_ribbon  <- data.frame(x = c(0, dataframe$x_num[1]), 
                              y = c(0, dataframe$Mean[1]), 
                              ymin = c(0, dataframe$Mean[1] - dataframe$sd[1]), 
                              ymax = c(0, dataframe$Mean[1] + dataframe$sd[1]))
    
    mean_plt_fun <- ggplot(data = dataframe, aes(x = x_num, y = Mean)) +
      annotate("segment", x = 0, y = 0, xend = dataframe$x_num[1], yend = dataframe$Mean[1], 
               color = line_color, lwd = 3) +
      geom_line(aes(group = 1), color = line_color, arrow = myarrow, lwd = 3) +
      geom_point(color = col[2:5], size = 5) + geom_point(color = col[1], size = 0.1)+
      geom_ribbon(aes(ymin = Mean - sd, ymax = Mean + sd), fill = ribbon_color,
                  alpha = 0.2) +
      geom_ribbon(data = data_ribbon, 
                  aes(x = x, y = y, ymin = ymin, ymax = ymax), 
                  fill = ribbon_color, alpha = 0.2) +
      scale_x_continuous(breaks = dataframe$x_num, labels = dataframe$Conditions,expand = c(0, 0.1)) +
      geom_hline(yintercept=0, linetype="dashed", color="black", linewidth=1.5)+
      ylab("Mean") +
      xlab("Conditions") +
      ylim(-1.1,1.1) +
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
                                        linetype = "solid"),
            plot.margin = unit(c(1, 1, 1, 1), "cm"))
    
    # print(mean_plt_fun)
    
    # Save the plot in the list
    plots_list[[paste0("k=", j)]] <- mean_plt_fun
    
    trend_dir <- file.path(cluster_dir,"Trends_per_cluster",folder_to_save)
    dir.create(trend_dir,showWarnings = F,recursive = T)
    
    pdf(file.path(trend_dir,paste0("k=",j,"_plt.pdf")),width = 10,height = 10)
    print(mean_plt_fun)
    dev.off()
    
  }
  
  
  # Calculate the number of plots
  num_plots <- length(plots_list)
  
  # Calculate the number of columns and rows to automatically adjust the layout
  ncol_plots <- ceiling(sqrt(num_plots))  # Square root approximation for plot distribution
  nrow_plots <- ceiling(num_plots / ncol_plots)  # Adjust based on the number of columns
  
  # Arrange all plots in a grid and remove x-axis labels for all but bottom row
  for (i in seq_along(plots_list)) {
    if (i <= (num_plots - ncol_plots)) {
      plots_list[[i]] <- plots_list[[i]] + theme(axis.text.x = element_blank(), 
                                                 axis.title.x = element_blank())
    }
  }
  
  rel_heights <- c(rep(1, nrow_plots - 1), 1.4)
  
  pdf(file.path(trend_dir,"grid_16_plt.pdf"),width = 32,height = 28)
  # Use do.call with grid.arrange to display the plots with the appropriate number of rows and columns
  # do.call(plot_grid, c(plotlist =plots_list, nrow = nrow_plots,ncol = ncol_plots,
  #                      rel_heights=rel_heights))
  plot_grid(plotlist = plots_list, align = "v", nrow = nrow_plots,
            ncol = ncol_plots,rel_heights=rel_heights)
  
  dev.off()
  
  
}


