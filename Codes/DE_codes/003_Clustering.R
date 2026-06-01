
#########################
## 0.Load Dependencies ##
#########################

bootstrap_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
repeat {
  setup_candidate <- file.path(bootstrap_root, "Codes", "_shared", "project_setup.R")
  if (file.exists(setup_candidate)) {
    source(setup_candidate)
    break
  }
  parent <- dirname(bootstrap_root)
  if (identical(parent, bootstrap_root)) {
    stop("Could not locate Codes/_shared/project_setup.R")
  }
  bootstrap_root <- parent
}

paths <- bootstrap_project()
project_root <- paths$project_root
input_dir <- paths$input_dir
outputs_root <- paths$outputs_root

required_packages <- c(
  "tidyr", "gridExtra", "dplyr", "ggthemes", "ggplot2", "ggrepel",
  "cowplot", "RColorBrewer", "ggdendro", "dendextend", "dendsort",
  "reshape2", "umap", "circlize", "ComplexHeatmap",
  "ConsensusClusterPlus", "factoextra", "NbClust", "mclust", "clustertend"
)
load_required_packages(required_packages)

{
  ### Optional clustering packages used in earlier exploration:
  # library(RCKS)
  # library(M3C)
}

### Configure project-relative paths
script_id  <- "003_Clustering"
output_dir <- file.path(outputs_root, script_id)
data_to_plot_dir <- file.path(output_dir,"002_Data_to_plot")
cluster_dir <- file.path(output_dir,"004_Clustering_plots")

dir.create(data_to_plot_dir,recursive = T,showWarnings = F)
dir.create(file.path(input_dir),recursive = T,showWarnings = F)
dir.create(file.path(output_dir,"001_Figures"),recursive = T,showWarnings = F)
dir.create(file.path(output_dir,"004_Clustering_plots"),showWarnings = F)

######################
##  Load Functions  ##
######################
{
  ### Function to calculate the BIC and the AIC for clustering algorithms
  kmeansAIC_BIC = function(fit){
    
    m = ncol(fit$centers)   # Number of conditions or samples
    n = length(fit$cluster) # Number of genes or rows
    k = nrow(fit$centers)   # Number of centers or clusters
    D = fit$tot.withinss    # Total within sum of squares
    return(data.frame(AIC = D + 2*m*k,
                      BIC = D + log(n)*m*k))
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
  
  ### Function to plot the information criterion values
  #ic_df         <- data frame with the number of clusters and the AIC and BIC values of each model
  #col_to_select <- string of the colname with the respective clustering model
  #ylab          <- y axis label as string
  #color         <- color of the points
  #color_opt     <- color orf the point representing the optimal number of clusters
  
  info_crit_plt <-function(ic_df,col_to_select,ylab ="IC",color_optimal = "#f03b20",color_points = "#2c7fb8" ) {
    
    min_ic <- min(ic_df[[col_to_select]])
    cluster_selected <- ic_df[which(ic_df[col_to_select]==min_ic),1]
    cat("the optimal number of clusters is:",cluster_selected)
    
    # Crear plot mejorado
    IC_plt <- ggplot(ic_df, aes(x = k.clusters, y = .data[[col_to_select]])) +
      geom_point(size = 3, color = color_points) +                       # todos los puntos
      geom_point(data = subset(ic_df, ic_df[[col_to_select]] == min_ic),
                 aes(x = k.clusters, y = !!sym(col_to_select)),
                 color = color_optimal, size = 5) +                      # punto óptimo
      geom_text(data = subset(ic_df, ic_df[[col_to_select]] == min_ic),
                aes(x = k.clusters, y = !!sym(col_to_select),
                    label = paste0("Optimal: ", k.clusters)),
                vjust = -1, size = 6, fontface = "bold", color = color_optimal) +
      xlab("Number of Clusters") +
      ylab(ylab) +
      theme_classic(base_size = 16) +
      theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black"),
        plot.background = element_blank()
      )
    
    # Guardar en PDF de alta resolución
    pdf(file.path(cluster_dir, paste0("001_", ylab, "_plt.pdf")), width = 12, height = 10)
    print(IC_plt)
    dev.off()
    
    list_to_return <- list(IC_plt,cluster_selected)
    return(list_to_return)
  }
  
  
}

##############
##  Step#1  ##
##############
{
  
  
  
  ### load data
  # The results from DE analysis are in the 002_Differential_expression output folder
  resul_dir <- file.path(outputs_root, "002_Differential_expression")
  results <- read.table(file.path(resul_dir,"003_DE_def","resultados.txt"))
  
  ################################################
  ##  1. Select Conditions stats and normalize  ##
  ################################################
  
  tab=results
  
  ### Select and inspect the data
  grep("beta",colnames(tab))
  colnames(tab)[c(1,7,13,19,25,31,37,43,49,55,61,67,73)]
  grep("hit",colnames(tab))
  colnames(tab)[grep("hit",colnames(tab))]
  grep("SE",colnames(tab))
  colnames(tab)[grep("SE",colnames(tab))]
  
  # Select these betas
  # "beta_inf_2"         "beta_inf_20"        "beta_inf_48"        "beta_inf_72"
  
  ### Build table with the infection logFCs of genes that respond to infection
  ###  at at least one timepoint. FDR< 0.05 and |logFC| >0.2
  
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
  ##########################################
  ##  1. Find optimal number of Clusters  ##
  ##########################################
  
  ##############################################
  ##  Run kmeans for kmax number of clusters  ##
  ##############################################
  ### Create a list to save the kmeans results
  kmeans_res <- list()
  cluster_list <- list()
  kmax <- 25
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
  
  # Plot the WSS to visualize the elbow
  plot(1:20, wss, type = "b", xlab = "Number of Clusters",
       ylab = "Within groups sum of squares")
  
  ##############################################
  
  ##  Compute AIC and BIC for kmeans and hckmeans models  ##
  
  ### Define the list to capture AIC and BIC outputs
  ic_values          <- list()
  ic_values_hckmeans <- list()
  
  ### Define the list to capture kmeans an hckmeans outputs
  kmeans_model   <- list()
  hckmeans_model <- list()
  
  ### Define the max number of clusters you want to inspect
  k_vec <- seq(40)
  
  for (k in k_vec) {
    
    kmeans_model[[k]]   <- kmeans(norm_betas,
                                  centers = k,
                                  nstart = k,
                                  iter.max = 500)
    hckmeans_model[[k]] <- hkmeans(norm_betas,
                                   k=k,
                                   hc.metric = "euclidean",
                                   hc.method = "ward.D2",
                                   iter.max = 500,
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
  
  ################################
  ##  Plot and save the graphs  ##
  ################################
  
  AIC_plt<- info_crit_plt(ic_df = ic_df,col_to_select = "AIC.Values",ylab = "AIC" )
  cluster_selected_AIC <- AIC_plt[[2]]
  # the optimal number of clusters is: 33
  
  BIC_plt <- info_crit_plt(ic_df = ic_df,col_to_select = "BIC.Values",ylab = "BIC",color_optimal = "#f03b20",color_points = "#2c7fb8")
  cluster_selected_BIC <- BIC_plt[[2]]
  # the optimal number of clusters is: 14
  
  AIC_hckmeans_plt <- info_crit_plt(ic_df = ic_df,col_to_select = "AIC.Values.hckmeans",ylab = "AIC_hckmeans" )
  cluster_selected_AIC_hckmeans <- AIC_hckmeans_plt[[2]]
  # the optimal number of clusters is: 29
  
  BIC_hckmeans_plt <- info_crit_plt(ic_df = ic_df,col_to_select = "BIC.Values.hckmeans",ylab = "BIC_hckmeans",color_optimal = "#f03b20",color_points = "#2c7fb8")
  cluster_selected_BIC_hckmeans <- BIC_hckmeans_plt[[2]]
  # the optimal number of clusters is: 15
  
  ### Save data to plot later
  write.table(ic_df,file.path(data_to_plot_dir,"001_AIC_BIC_values_per_k.txt"))
  
  ### check the order of the gene names in each data
  length(which(rownames(norm_betas) != names(kmeans_model[[24]]$cluster)))
  length(which(rownames(norm_betas) != names(hckmeans_model[[22]]$cluster)))
  
  names_of_df_col <- c(paste0("kmeans_cluster_k=",cluster_selected_BIC),
                       paste0("kmeans_cluster_k=",cluster_selected_AIC),
                       paste0("hckmeans_cluster_k=",cluster_selected_BIC_hckmeans),
                       paste0("hckmeans_cluster_k=",cluster_selected_AIC_hckmeans))
  
  ### Create the final table with the clustering information
  
  
  cluster_out <- data.frame(norm_betas,
                            kmeans   = kmeans_model[[cluster_selected_BIC]]$cluster,
                            kmeans2  = kmeans_model[[cluster_selected_AIC]]$cluster,
                            hckmeans = hckmeans_model[[cluster_selected_BIC_hckmeans]]$cluster,
                            hckmeans2= hckmeans_model[[cluster_selected_AIC_hckmeans]]$cluster
  )
  
  colnames(cluster_out) <-c(colnames(norm_betas), names_of_df_col)
  
  cluster_input_dir <- file.path(input_dir,"003_Clustering_tabs")
  dir.create(cluster_input_dir,recursive = T,showWarnings = F)
  
  # Save the kmeans and hckmeans objects
  saveRDS(kmeans_model, file = file.path(cluster_input_dir,"001_kmeans_clustering.RDS"))
  saveRDS(hckmeans_model, file = file.path(cluster_input_dir,"001_hckmeans_clustering.RDS"))
  
  write.table(cluster_out,file.path(cluster_input_dir,"lfc_with_clustering.txt"))
  
}

##############
##  Step#2  ##
##############
{
  ##########################################################
  # Plot the heatmaps for k clusters based on AIC and BIC ##
  ##########################################################
  
  # cluster_input_dir <- file.path(input_dir,"003_Clustering_tabs")
  # cluster_out <- read.table(file.path(cluster_input_dir,"lfc_with_clustering.txt"))
  
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
  n_colors <- 35  # Total de colores que quieres
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
    clust_name   = "Genes"
    cluster_colors_ha <- cluster_colors[1:max(as.numeric(table[,5]))]
    split <- as.numeric(table[,5])
    order_rows <- rownames(table)
    length(which(rownames(hclust_matrix)!=rownames(table)))
    
    ha <- rowAnnotation(Cluster = as.factor(table[,5]),col = list(Cluster = cluster_colors_ha),
                        annotation_legend_param = list(
                          title = "Cluster",
                          labels = unique(table[,5])  # This includes the clsuter numbers present in the data
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
                      # cluster_rows = T,
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

##############
##  Step#3  ##
##############
{
  ######################################
  ##  Compute the trends per cluster  ##
  ######################################
  
  # Compute the mean and the sd by cluster
  result_mean_list <- list() # list to save
  result_sd        <- list() # list to save the sd
  result_mean      <- list() # list to save the mean
  result_mean_to_save   <- list()
  
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
      
      result_mean_to_save[[j]] <- dataframe
      
      myarrow <- arrow(angle = 20, length = unit(0.25, "cm"), ends = "last", type = "closed")
      
      col <- c("#D73027", "#4575B4", "#74ADD1", "#FDAE61")  
      line_color <- "black"
      ribbon_color <- "#4575B4"
      
      mean_plt_fun <- ggplot(data = dataframe, aes(x = x_num, y = Mean)) +
        # Standard deviation ribbon
        geom_ribbon(aes(ymin = Mean - sd, ymax = Mean + sd), 
                    fill = ribbon_color, alpha = 0.15) +
        # Average line with arrow
        geom_line(aes(group = 1), color = line_color, linewidth = 1.2,
                  arrow = myarrow, lineend = "round") +
        # Points of each condition
        geom_point(color = col, size = 4, stroke = 0.6, shape = 21, fill = "white") +
        # Axis and labels
        scale_x_continuous(breaks = dataframe$x_num, labels = dataframe$Conditions) +
        ylim(-1.1,1.1)+
        labs(
          y = "Mean expression ± SD",
          x = "Condition",
          title = paste0("Cluster mean trend (k = ", j, ")") ) +
        theme_minimal(base_size = 16, base_family = "Helvetica") +
        theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 14, angle = 45, vjust = 0.9, hjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
          axis.line = element_line(color = "black", linewidth = 0.6),
          plot.margin = margin(10, 10, 10, 10)
        )
      
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
    
    pdf(file.path(trend_dir,"grid_plt.pdf"),width = 40,height = 38)
    # Use do.call with grid.arrange to display the plots with the appropriate number of rows and columns
    do.call(grid.arrange, c(plots_list, ncol = ncol_plots, nrow = nrow_plots))
    dev.off()
  }
  
  saveRDS(result_mean_to_save, file.path(cluster_input_dir,"002_df_Mean_sd_by_cluster.RDS"))
  saveRDS(result_mean_list,    file.path(cluster_input_dir,"002_Mean_sd_by_cluster.RDS"))
  
}

##############
##  Step#4  ##
##############
{
  
  #################################################
  ##  Plot Heatmaps between AIC and BIC results  ##
  #################################################
  
  hckmeans_selected <- hckmeans_model[cluster_selected_BIC_hckmeans:cluster_selected_AIC_hckmeans]
  cluster_out1 <- norm_betas
  
  # Loop to add selected columns from each data frame in the list to the cluster_out1 df
  for (i in seq_along(hckmeans_selected)) {
    # Select the column you want to add, e.g., cluster
    col_to_add <- hckmeans_selected[[i]][["cluster"]]
    
    # Add the column to cluster_out1 using the k clusters as a name
    cluster_out1[[paste0("k_", i+cluster_selected_BIC_hckmeans-1)]] <- col_to_add
  }
  
  vec <- 5:ncol(cluster_out1)
  
  for(i in vec){
    
    lfc_normalized <- cluster_out1[,c(1:4,i)]
    
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
                      cluster_row_slices = TRUE,
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
  
  write.table( cluster_out1, file.path(cluster_input_dir,"lfc_data_with_clusters_information.txt"))
  write.table( AIC_values, file.path(cluster_input_dir,"AIC_kmeans.txt"))
  write.table( BIC_values, file.path(cluster_input_dir,"BIC_kmeans.txt"))
  write.table( AIC_values_hckmeans,file.path(cluster_input_dir,"AIC_hckmeans.txt"))
  write.table( BIC_values_hckmeans,file.path(cluster_input_dir,"BIC_hckmeans.txt"))
  
  
}

##############
##  Step#5  ##
##############
{
  ###############################################################################
  ##  We selected k=16 from hckmeans. Now we plot AIC and BIC, and the trends  ##
  ###############################################################################
  
  cluster_out1 <- read.table(file.path(cluster_input_dir,"lfc_data_with_clusters_information.txt"))
  
  BIC_hckmeans_plt <- info_crit_plt(ic_df = ic_df[2:nrow(ic_df),],col_to_select = "BIC.Values.hckmeans",
                                    ylab = "BIC_hckmeans_wo_1st_point",color_optimal = "#f03b20",color_points = "#2c7fb8")
  cluster_selected_BIC_hckmeans <- BIC_hckmeans_plt[[2]]
  
  AIC_hckmeans_plt <- info_crit_plt(ic_df = ic_df[2:nrow(ic_df),],col_to_select = "AIC.Values.hckmeans",
                                    ylab = "AIC_hckmeans_wo_1st_point",color_optimal = "#f03b20",color_points = "#2c7fb8")
  cluster_selected_AIC_hckmeans <- AIC_hckmeans_plt[[2]]
  
  opt_AIC <- cluster_selected_AIC_hckmeans
  opt_BIC <- cluster_selected_BIC_hckmeans
  chosen_cluster <- 16
  
  combined_IC_plot <- ggplot(ic_df, aes(x = k.clusters)) +
    # AIC
    geom_point(aes(y = AIC.Values.hckmeans), color = "#2c7fb8", size = 3) +
    # BIC points
    geom_point(aes(y = BIC.Values.hckmeans), color = "#f03b20", size = 3) +
    # Optimal AIC point
    geom_point(data = subset(ic_df, k.clusters == opt_AIC),
               aes(y = AIC.Values.hckmeans),
               color = "#238b45", size = 6) +   # verde, triángulo
    # Optimal BIC point
    geom_point(data = subset(ic_df, k.clusters == opt_BIC),
               aes(y = BIC.Values.hckmeans),
               color = "#e6550d", size = 6) +   # naranja, rombo
    # Vertical line at chosen cluster
    geom_vline(xintercept = chosen_cluster, 
               linetype = "dashed", color = "black", linewidth = 1) +
    xlab("Number of Clusters") +
    ylab("Information Criterion") +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16, face = "bold"))
  
  # Guardar en PDF de alta resolución
  pdf(file.path(cluster_dir, "002_AIC_BIC_comparison.pdf"),width = 14, height = 8, useDingbats = FALSE)
  print(combined_IC_plot)
  dev.off()
  
  ###################
  ##  Trend plots  ##
  ###################
  
  {
    # Compute the mean and the sd by cluster
    result_mean_list <- list() # list to save
    result_sd        <- list() # list to save the sd
    result_mean      <- list() # list to save the mean
    result_mean_to_save <- list()
    
    df <- cluster_out1[,c(1:4,6)] # 6 is the cluster k=16
    folder_to_save <- paste0("hckmeans_",colnames(df)[5])
    df_por_clusters <- split(df, df[,5])
    names(df_por_clusters) <- paste0("k=",names(df_por_clusters))
    
    # list to save plt to create a grid
    plots_list <- list()
    
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
      result_mean_to_save[[j]] <- dataframe
      # dataframe <- rbind(data.frame(Conditions = "Start", Mean = 0, sd = 0, x_num = 0), dataframe)
      
      myarrow <- arrow(angle = 20, length = unit(0.25, "cm"), ends = "last", type = "closed")
      
      # Define the line color and SD shading
      col <- c("#D73027", "#4575B4", "#74ADD1", "#FDAE61")  
      line_color <- "black" # Change this value to the desired line color
      ribbon_color <- "#4575B4"# Change this value to the desired shading color
      
      data_ribbon <- data.frame(x = c(0, dataframe$x_num[1]),
                                y = c(0, dataframe$Mean[1]),
                                ymin = c(0, dataframe$Mean[1] - dataframe$sd[1]),
                                ymax = c(0, dataframe$Mean[1] + dataframe$sd[1]))
      
      mean_plt_fun <- ggplot(data = dataframe, aes(x = x_num, y = Mean)) +
        # Standard deviation ribbon
        geom_ribbon(aes(ymin = Mean - sd, ymax = Mean + sd), 
                    fill = ribbon_color, alpha = 0.15) +
        # Average line with arrow
        geom_line(aes(group = 1), color = line_color, linewidth = 1.2,
                  arrow = myarrow, lineend = "round") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray10") +
        # Points of each condition
        geom_point(color = col, size = 4, stroke = 0.6, shape = 21, fill = "white") +
        # Axis and labels
        scale_x_continuous(breaks = dataframe$x_num, labels = dataframe$Conditions) +
        ylim(-1.1,1.1)+
        labs(
          y = "Mean expression ± SD",
          x = "Condition",
          title = paste0("Cluster mean trend (k = ", j, ")") ) +
        theme_minimal(base_size = 16, base_family = "Helvetica") +
        theme(
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(size = 14, angle = 45, vjust = 0.9, hjust = 1),
          axis.text.y = element_text(size = 14),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray80", linewidth = 0.3),
          axis.line = element_line(color = "black", linewidth = 0.6),
          plot.margin = margin(10, 10, 10, 10)
        )
      
      # print(mean_plt_fun)
      
      # Save the plot in the list
      plots_list[[paste0("k=", j)]] <- mean_plt_fun
      
      trend_dir <- file.path(cluster_dir,"Trends_per_cluster",folder_to_save)
      dir.create(trend_dir,showWarnings = F,recursive = T)
      
      pdf(file.path(trend_dir,paste0("k=",j,"_plt.pdf")),width = 10,height = 10)
      print(mean_plt_fun)
      dev.off()
      
    }
    
    
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
  
  pdf(file.path(trend_dir,"grid_plt.pdf"),width = 32,height = 28)
  # Use do.call with grid.arrange to display the plots with the appropriate number of rows and columns
  # do.call(plot_grid, c(plotlist =plots_list, nrow = nrow_plots,ncol = ncol_plots,
  #                      rel_heights=rel_heights))
  plot_grid(plotlist = plots_list, align = "v", nrow = nrow_plots,
            ncol = ncol_plots,rel_heights=rel_heights)
  
  dev.off()
  
  saveRDS(result_mean_to_save, file.path(cluster_input_dir,"002_df_Mean_sd_by_cluster_k_16.RDS"))
  
}
