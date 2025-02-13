#############################
### 0. Load dependencies. ###
#############################

{
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
}

############################
### 1. Declare functions ### 
############################
dcols=function(x){data.frame(colnames(x))}

###################################################
### 2. Set working directory and create folders ###
###################################################
main_wd <- getwd()
output_dir <- "Analyses/Outputs"
input_dir  <- "Analyses/Inputs/2_Processed_data"
cluster_input_dir <- file.path("Analyses","Inputs","003_Clustering_tabs")

#################
##  Load Data  ##
#################

# Load the DE results
results <- read.table(file.path(output_dir,"003_DE_def","resultados.txt"))
cluster_out <- read.table(file.path(cluster_input_dir,"lfc_data_with_clusters_information.txt"))

#We selected the k=16 cluster results
cluster_out1 <- cluster_out[,c(1:4,7)]
genes.names <- rownames(cluster_out1)### Gene names with Ensemble IDs (Novel Proteins)
list_per_clusters <- split(rownames(cluster_out1), cluster_out1[,5])
names(list_per_clusters) <- paste0("k=",names(list_per_clusters))

file_name_clusters <- names(list_per_clusters)
### Save the clusters again
lapply(seq_along(list_per_clusters), function(i) {
  file_path <- file.path(cluster_input_dir, paste0(file_name_clusters[i],".txt"))  
  writeLines(list_per_clusters[[i]], file_path)  
})

#Inspect the data
grep("beta",colnames(results))
colnames(results)[c(1,7,13,19,25,31,37,43,49,55,61,67,73)]
grep("hit",colnames(results))
colnames(results)[grep("hit",colnames(results))]
grep("SE",colnames(results))
colnames(results)[grep("SE",colnames(results))]

df=results[,c("beta_inf_2","beta_inf_20","beta_inf_48","beta_inf_72","hit_inf_2","hit_inf_20","hit_inf_48","hit_inf_72")]

infected_2h  <- rownames(df[which(df$hit_inf_2==1),])
infected_20h <- rownames(df[which(df$hit_inf_20==1),])
infected_48h <- rownames(df[which(df$hit_inf_48==1),])
infected_72h <- rownames(df[which(df$hit_inf_72==1),])

list_DE_infected_overtime <- list(infected_2h=infected_2h,infected_20h=infected_20h,
                                  infected_48h=infected_48h,infected_72h=infected_72h)
# Save the list with the DE expressed gene names over time
saveRDS(list_DE_infected_overtime,file.path(cluster_input_dir,"list_DE_infected_overtime.RDS"))
file_name <- names(list_DE_infected_overtime)

# Save the gene names by conditions separately
writeLines(list_DE_infected_overtime[[1]], file.path(cluster_input_dir,file_name[1]))
writeLines(list_DE_infected_overtime[[2]], file.path(cluster_input_dir,file_name[2]))
writeLines(list_DE_infected_overtime[[3]], file.path(cluster_input_dir,file_name[3]))
writeLines(list_DE_infected_overtime[[4]], file.path(cluster_input_dir,file_name[4]))
  



