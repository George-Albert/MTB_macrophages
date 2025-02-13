#########################
## 0.Load Dependencies ##
#########################

{
  library(tidyverse)
  library(dplyr)
  library(ggthemes)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
  library(RColorBrewer)
  library(qvalue)
  library(igraph)
  library(xlsx)
  library(circlize)
  library(ggsci)
}

#########################
##  Declare functions  ##
#########################

{
  # See the column structure as a df
  dcols=function(x){data.frame(colnames(x))}
  
  # Measure the size of the Heatmap plot
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
  
  # Create the gexf file to create a network in Gephi
  write_gephi_gexf <- function(GO_object,label_indexes,output_file) {
    
    df_nodes_all=GO_object$result
    df_nodes_all$label=""
    set_base=which(!is.na(df_nodes_all$community))
    df_nodes_all$label[set_base[label_indexes]]=df_nodes_all$description[set_base[label_indexes]]
    df_nodes_all=df_nodes_all[,c(1:2,ncol(df_nodes_all),3:(ncol(df_nodes_all)-1))]
    GO_object_out=GO_object
    GO_object_out$result=df_nodes_all
    
    df_nodes=GO_object$result
    df_nodes=df_nodes[which(!is.na(df_nodes$community)),]
    df_nodes$label=""
    df_nodes$label[label_indexes]=df_nodes$description[label_indexes]
    df_nodes=df_nodes[,c(1:2,16,3:15)]
    df_edges=GO_object$network
    
    nodes <- data.frame(ID=df_nodes$GO_ID,Label=df_nodes$label)
    edges=df_edges[,1:2]
    edgesWeight=as.numeric(df_edges$Weight)
    
    #df_edges$edge_ID=c(1:nrow(df_edges))
    #edgesId=data.frame(df_edges$edge_ID)
    #edgesAtt=data.frame(df_edges[,4,3])
    df_nodes$significance <- -log10(df_nodes$fdr)
    nodesAtt=df_nodes
    
    # Escribir el archivo .gexf
    write.gexf(nodes=nodes,edges=edges,edgesWeight=edgesWeight,nodesAtt=nodesAtt,output=paste0(output_file,".gexf"))
    return(GO_object_out)
  }
  
  # Prepare the output
  output_preparer <- function(data_input) {
    data_input$a="A"
    data_input$b="B"
    
    for(i in 1:nrow(data_input)) {
      data_input$a[i] <- paste(data_input$Levels[[i]], collapse = ", ")
      data_input$b[i] <- paste(data_input$Members[[i]][1][[1]], collapse=", ")
    }
    data_input$Levels=data_input$a
    data_input$Members=data_input$b
    data_input=data_input[,1:(ncol(data_input)-2)]
    return(data_input)
  }
  
  # Build the network
  build_network <- function(tab, threshold) {
    links <- data.frame(node_source = character(0), node_target = character(0), P = numeric(0), stringsAsFactors = FALSE)
    
    for (i in 1:(nrow(tab))) {
      for (j in (i):nrow(tab)) {
        elements1 <- unlist(strsplit(tab$Members[i], ","))
        elements2 <- unlist(strsplit(tab$Members[j], ","))
        
        common_elements <- intersect(elements1, elements2)
        P <- length(common_elements) / min(length(elements1), length(elements2))
        
        if (P > threshold) {
          links <- rbind(links, c(tab[i, "Ontology"], tab[j, "Ontology"], P))
        }
      }
    }
    
    colnames(links) <- c("GO_source", "GO_target", "Weight")
    return(links)
  }
  
  # GO function 
  GO=function(query,background_matrix,metadata_matrix,bg=NA,min_GO_level=0,
              max_GO_level=15,test_mode="enrichment",pvalue_correction_method="BH",
              min_term_size = 1,max_term_size = Inf,th_hit=0.1,th_hit_OR=2,
              threshold_links =0.5,merge_duplicates=TRUE){
    
    ## Quality check 1: Are there cols, or rows, duplicated in the background matrix? 
    ## (In the row names, duplicates are not allowed). Are there genes in the query?.
    
    unicos=length(unique(colnames(background_matrix)))
    todos=ncol(background_matrix)
    if(unicos!=todos)
    {
      print(paste("Warning: repeated terms in the background matrix:",todos," submitted, ",unicos," are unique"))
      background_matrix=background_matrix[,which(!duplicated(colnames(background_matrix)))]
    }
    
    unicos=length(unique(rownames(background_matrix)))
    todos=nrow(background_matrix)
    if(unicos!=todos)
    {
      print(paste("Warning: repeated genes in the background matrix:",todos," submitted, ",unicos," are unique"))
      background_matrix=background_matrix[which(!duplicated(rownames(background_matrix))),]
    }
    
    unicos=length(unique(query))
    todos=length(query)
    if(unicos!=todos)
    {
      print(paste("Warning: repeated genes in query:",todos," submitted, ",unicos," are unique"))
      query=query[which(!duplicated(query))]
    }
    
    ## Si había background:
    ## 1. Completo QC: bg no contiene duplicados, y todo el query esta dentro del bg.
    ## 2. Subset the matrix with annotations row-wise: only genes within background 
    ## should be in. AND expand it to include genes in BG not annotated.
    
    ## Si NO había background:
    ## 1. Define the background como el set de genes con alguna anotación EN 
    ## CUALQUIER NIVEL. (Idea: si un gen tiene una anotación fuera de los 
    ## niveles target, aunque no la tenga en ninguno de los términos de los
    ## niveles target, podría haberla tenido, y por tanto puede y debe entrar en el bg)
    
    ## 2. Subset the query to live in that background.
    
    ## Finally in both cases:
    ## 3. Subset the annotation matrix to contain the terms that will be 
    ## tested, and only those. Declare also query_vector.
    ## 4. locate terms containing exactly the same set of genes in the
    ##  defined universe through a uniqueness_index attribute for terms
    
    if(length(bg)>1)
    {
      ## 1. completamos QC: dupes in BG??
      unicos=length(unique(bg))
      todos=length(bg)
      if(unicos!=todos)
      {
        print(paste("Warning: repeated genes in bg:",todos," submitted, ",unicos," are unique"))
        bg=bg[which(!duplicated(bg))]
      }
      ## 1. Completamos QC: Is query within bg?
      within=length(which(query %in% bg))
      todos=length(query)
      if(within!=todos)
      {
        print(paste("Warning: from ",todos," submitted in the query, only ",within," are in the bg, and therefore will be tested"))
        query=query[which(query %in% bg)]
      }
      #2. Subset the matrix with annotations row-wise:only genes within background 
      #should be in...
      background_matrix=background_matrix[which(rownames(background_matrix) %in% bg),]
      extra_genes=bg[which(!bg %in% rownames(background_matrix))]
      #2...AND expand it to include genes in BG not annotated.
      appendix <- data.frame(matrix(FALSE, nrow = length(extra_genes), 
                                    ncol = ncol(background_matrix)))
      colnames(appendix) <- colnames(background_matrix)
      rownames(appendix) <- extra_genes
      background_matrix=rbind(background_matrix,appendix)
      
      # Check that 2 worked well:
      background_matrix=background_matrix[order(rownames(background_matrix)),]
      bg=bg[order(bg)]
      #print("Vio bg")
    }else{
      # 1. define bg from annotated genes in all levels.
      genes_presence_vector=apply(background_matrix,1,sum)
      annotated_genes=rownames(background_matrix)[which(genes_presence_vector>0)]
      background_matrix=background_matrix[annotated_genes,]
      background_matrix=background_matrix[order(rownames(background_matrix)),]
      bg=rownames(background_matrix)
      ## 2. Subset the query to live in that background.
      query=query[which(query %in% bg)]
      #print("NO VIO bg")
    }
    
    ##3. Subset the annotation matrix to contain the terms that will  
    ##be tested, and only those.
    
    #query_vector=data.frame(background_matrix[,1])
    #query_vector[,1]=FALSE
    #colnames(query_vector)="query"
    #rownames(query_vector)=rownames(background_matrix)
    #query_vector[which(rownames(query_vector) %in% query),1]=TRUE
    
    
    ## filter terms: 1st according to levels:
    GO_levels <- min_GO_level:max_GO_level
    mode <- switch(test_mode, "enrichment" = "greater", "depletion"="less", "two.sided"="two.sided")
    terms_within_levels <- c()
    
    for (i in 1:nrow(metadata_matrix)) {
      if (sum(metadata_matrix$Levels[i][[1]] %in% GO_levels)>0) {
        terms_within_levels <- c(terms_within_levels, i)
      }
      # else{
      #   # print(paste0("GO_levels are less than or equal to 0"))
      # }
    }
    
    metadata_matrix <- metadata_matrix[terms_within_levels,]
    background_matrix=background_matrix[,terms_within_levels]
    ## filter terms: 2nd according to term size:
    
    number_genes_per_term=apply(background_matrix,2,sum)
    metadata_matrix$Size_in_universe=number_genes_per_term
    set_passing=which(number_genes_per_term>=min_term_size & number_genes_per_term<=max_term_size)
    
    metadata_matrix <- metadata_matrix[set_passing,]
    background_matrix=background_matrix[,set_passing]
    
    #4. Identify terms that have exactly the same annotations using a uniqueness_index column
    
    colapsos=apply(background_matrix,2,function(x){paste(x,collapse="")})
    set_first_instances=which(!duplicated(colapsos))
    metadata_matrix$uniqueness_index=0
    metadata_matrix$uniqueness_index[set_first_instances]=c(1:length(set_first_instances))
    
    dupe_instances=unique(colapsos[which(duplicated(colapsos))])
    
    for(i in 1:length(dupe_instances))
    {
      set=which(colapsos==dupe_instances[i])
      value=max(metadata_matrix$uniqueness_index[set])
      metadata_matrix$uniqueness_index[set]=value
    }
    
    ## Now, we load the number of genes in query and the percentage of the total 
    ## term that represent it
    query_matrix <- background_matrix[which(rownames(background_matrix) %in% query),]
    
    number_genes_vector <- c()
    percentage_genes_vector <- c()
    
    for (i in 1:nrow(metadata_matrix)) {
      
      number_of_genes <- sum(query_matrix[,metadata_matrix[i,2]])
      percentage <- (number_of_genes/metadata_matrix$Size_in_universe[i]) * 100
      
      number_genes_vector <- c(number_genes_vector, number_of_genes)
      percentage_genes_vector <- c(percentage_genes_vector, percentage)
    }
    
    metadata_matrix$Number_of_genes <- number_genes_vector
    metadata_matrix$Percentage <- percentage_genes_vector
    
    ## Omitiré filtrado por num genes en la intersección, percentage of 
    ## genes in the interseccion y tb por terminos con tamaño minimo.
    pval_vector <- c()
    OR_vector <- c()
    min_conf_int_vector <- c()
    max_conf_int_vector <- c()
    members_list <- list()
    
    for (i in  1:nrow(metadata_matrix)){
      
      in_both           <- sum(query_matrix[,metadata_matrix$Category[i]])
      in_term_out_query <- metadata_matrix$Size_in_universe[i] - in_both
      out_term_in_query <- length(query) - in_both
      out_both          <- nrow(background_matrix) - in_both - in_term_out_query - out_term_in_query
      
      contingency_table <- data.frame(
        c(in_both, in_term_out_query),
        c(out_term_in_query,out_both)
      )
      
      # od_ratio <- oddsratio(contingency_table)
      test=fisher.test(contingency_table,alternative=mode)
      pval_vector[i] <- test$p.value
      OR_vector[i] <- test$estimate
      min_conf_int_vector[i] <- test$conf.int[1]
      max_conf_int_vector[i] <- test$conf.int[2]
      
      members <- rownames(query_matrix)[which(query_matrix[,metadata_matrix$Category[i]])]
      members_list[[i]] <- list(members)
    }
    
    metadata_matrix$pvalue <- pval_vector
    metadata_matrix$Odds_Ratio <- OR_vector
    metadata_matrix$Members <- members_list
    metadata_matrix$Min_Confidence_Int <- min_conf_int_vector
    metadata_matrix$Max_Confidence_Int <- max_conf_int_vector
    metadata_matrix[order(metadata_matrix$pvalue),]
    #And now we apply the pvalue correction
    print("Nominal p values distribution summary:")
    print(summary(metadata_matrix$pvalue))
    #metadata_matrix <- metadata_matrix[which(is.na(metadata_matrix$pvalue) == FALSE & is.infinite(metadata_matrix$pvalue) == FALSE),]
    
    
    if(pvalue_correction_method=="qval"){
      set=which(!duplicated(metadata_matrix$uniqueness_index))
      metadata_matrix$P.adj=2 # why two
      # qs=qvalue(metadata_matrix$pvalue[set])
      qs=p.adjust(metadata_matrix$pvalue[set],method = "BH")
      # pi0=qs$pi0
      # metadata_matrix$P.adj[set]=qs$qvalues
      metadata_matrix$P.adj[set]=qs
    }else{
      set=which(!duplicated(metadata_matrix$uniqueness_index))
      metadata_matrix$P.adj=2
      metadata_matrix$P.adj[set]=p.adjust(metadata_matrix$pvalue[set], method=pvalue_correction_method)
      pi0=1
    }
    ## For the sets that are dupes, I copy the same FDR along the entire term: these are not many different terms, only one that is repeated.
    uniqueness_indexes_w_several_terms=unique(metadata_matrix$uniqueness_index[which(duplicated(metadata_matrix$uniqueness_index))])
    for(i in uniqueness_indexes_w_several_terms){
      set=which(metadata_matrix$uniqueness_index==i)
      metadata_matrix$P.adj[set]=min(metadata_matrix$P.adj[set])
    }
    
    ## Order things:
    background_matrix=background_matrix[,order(metadata_matrix$P.adj)]
    metadata_matrix=metadata_matrix[order(metadata_matrix$P.adj),]
    
    length(which(colnames(background_matrix)!=metadata_matrix$Category))
    rownames(metadata_matrix)=metadata_matrix$Ontology
    metadata_matrix <- output_preparer(metadata_matrix)
    metadata_matrix=metadata_matrix[,c(1:11,15,12:14)]
    
    # And finally we output the results: the background matrix containing only the
    # genes and terms considered in the analyses, the metadata table with the results and pi0
    print("P.adj values distribution summary:")
    print(summary(metadata_matrix$P.adj))
    
    hits_table=metadata_matrix[which(metadata_matrix$P.adj<th_hit & metadata_matrix$Odds_Ratio>th_hit_OR),]
    if(merge_duplicates) hits_table=hits_table[which(!duplicated(hits_table$uniqueness_index)),]
    print(paste0("The rows of the hit table that fulfill P.adj < ",th_hit," and Odds_Ratio > ",th_hit_OR,
                 " are:",dim(hits_table)[1], " GO_terms")) # Muestra el número de filas
    if (nrow(hits_table)>0 ) {
      
      network <- build_network(hits_table, threshold=threshold_links)
      #print(network)
      graph <- graph_from_data_frame(network, directed = FALSE)
      community <- cluster_louvain(graph)
      
      #length(which(community$names==hits_table$Ontology))
      hits_table$community=community$membership
      hits_table=hits_table[order(hits_table$community,hits_table$P.adj),]
      hits_table$index=c(1:nrow(hits_table))
      hits_table=hits_table[,c(1:12,16:17,14:15,13)]
      
      ## This is selected manually: these indexes are actually not good for this case
      
      results_output=merge(metadata_matrix,hits_table,by=colnames(metadata_matrix),all=TRUE)
      results_output=results_output[,c(1:12,16:17,14:15,13)]
      
      colnames(results_output)=c("GO_ID","description","levels","size_in_annotation","ontology_group",
                                 "size_in_universe","uniqueness_index","number_of_genes_in_query",
                                 "percentage_of_genes_in_query","p", "OR","fdr","community","node_numeric_id",
                                 "Min_value_CI","Max_value_CI","genes_in_query")
      results_output=results_output[,c(1,2,3,5,14,13,7,4,6,8,9,11,10,12,15:17)]
      results_output=results_output[order(results_output$community,results_output$fdr),]
      
      factor_original <- factor(results_output$uniqueness_index,levels=unique(results_output$uniqueness_index))
      niveles <- levels(factor_original)
      niveles_numeros <- seq_along(niveles)
      factor_numeros <- factor(factor_original, levels = niveles, labels = niveles_numeros)
      
      results_output$uniqueness_index=factor_numeros
      results_output=results_output[order(results_output$uniqueness_index,results_output$fdr),]
      
      output <- list(background_matrix=background_matrix,result=results_output,network=network)
      return(output)
      
    }else{
      print("There are no hits (rows) in `hits_table` because the `P.adj` and odds ratio values do not meet the specified thresholds.")
    }
    
  }
  
  # Heatmap plot
  heatmaplotter <- function (namelist,data,feature_data,columns=NULL,name){
    if (is.null(columns)) {
      data_mod<-data[which(rownames(data) %in% namelist),]
    }else { data_mod<-data[which(rownames(data) %in% namelist),columns]
    }
    data_mod<-data_mod[match(namelist,rownames(data_mod)),]
    cols=colnames(data_mod)
    
    data_mod=t(apply(data_mod,1,scale))
    colnames(data_mod)=cols
    data_mod<-data_mod[match(namelist,rownames(data_mod)),]
    
    # meanvec<-apply(data_mod, 1, mean,na.rm=TRUE)
    #for(i in 1:nrow(data_mod)){
    #     data_mod[i,]<-data_mod[i,]-meanvec[i]
    #}
    data_mod<-as.matrix(data_mod)
    chunk=feature_data[which(feature_data$Locus.Tag %in% rownames(data_mod)),]
    chunk<-chunk[match(rownames(data_mod),chunk$Locus.Tag),]
    rownames(data_mod)=chunk$Gene_name
    
    data_expanded=melt(data_mod)
    colnames(data_expanded)=c("Gene","Sample","Z")
    
    finalplot<-ggplot(data_expanded, aes(x=Sample,y=Gene, fill= Z)) +
      geom_tile(color = "white",
                lwd = 0.5,
                linetype = 1) +
      scale_fill_gradient2(low="blue", high="red")+coord_fixed()+ggtitle(name)
    return(finalplot)
  }
  
  # Names to select
  names_select=function(tab,attribute,values,features)
  {
    genes=tab$genes_in_query[which(tab[,attribute] %in% values)]
    genes <- paste(genes, collapse = ", ")
    genes <- unique(unlist(strsplit(genes, ", ")))
    print(length(genes))
    features=features[which(features$Locus.Tag %in% genes),]
    if(length(genes)!=nrow(features))
    {
      print("Se pierden alguno de los ",length(genes)," genes")
    }
    return(features)
  }
  
}

##################################
##  Load GO annotation objects  ##
##################################

{
  
  output_dir <- "Analyses/Outputs"
  input_dir  <- "Analyses/Inputs"
  
  enrichment_dir <- file.path(input_dir,"005_Enrichment")
  dir.create(enrichment_dir,showWarnings = F)
  input_folder <- file.path(enrichment_dir,"GO_terms")
  dir.create(file.path(input_folder),recursive = T,showWarnings = F)
  
  GO_gene_term_matrix_bp <- readRDS(file.path(input_folder,"boolean_matrix_bp.RDS"))
  BP_metadata_matrix     <- readRDS(file.path(input_folder,"BP_metadata_matrix.RDS"))
  GO_gene_term_matrix_cc <- readRDS(file.path(input_folder,"boolean_matrix_cc.RDS"))
  CC_metadata_matrix     <- readRDS(file.path(input_folder,"CC_metadata_matrix.RDS"))
  GO_gene_term_matrix_mf <- readRDS(file.path(input_folder,"boolean_matrix_mf.RDS"))
  MF_metadata_matrix     <- readRDS(file.path(input_folder,"MF_metadata_matrix.RDS"))
  
  ##------------------------------
  ##  Load GO annotation objects  
  ##------------------------------
  
  print(paste("Dimensions of the matrix with information about terms",
              "related to Cellular Component (CC) and their genes are", 
              dim(GO_gene_term_matrix_cc)[1],dim(GO_gene_term_matrix_cc)[2]))
  ###This matrix contains information about terms related to cellular component 
  ###and their genes genes are 19678 1989
  
  dim(CC_metadata_matrix) #This matrix contains metadata about the cellular 
  # component terms
  # 1989    5
  
  dim(GO_gene_term_matrix_bp) #This matrix contains information about terms 
  # related to biological process and their genes
  # 19678 16372
  
  dim(GO_gene_term_matrix_mf) #This matrix contains information about terms 
  # related to mollecular functions  and their genes
  # 20067  5114
  
  dim(BP_metadata_matrix)     #This matrix contains metadata about the 
  # biological process terms
  # 16372     5
  
  dim(MF_metadata_matrix) #This matrix contains metadata about the 
  # molecular functions terms
  # 5114    5
  
  # Check if the genes are sorted
  length(which(rownames(GO_gene_term_matrix_cc)!=rownames(GO_gene_term_matrix_bp)))
  # length(which(rownames(GO_gene_term_matrix_cc)!=rownames(GO_gene_term_matrix_mf)))
  
  # Check if the number of columns (Category) of the matrix are the same as 
  # the number of rows of the respective metadata
  length(which(colnames(GO_gene_term_matrix_cc)!=CC_metadata_matrix$Category))
  length(which(colnames(GO_gene_term_matrix_bp)!=BP_metadata_matrix$Category))
  
  GO_gene_term_matrix_all=cbind(GO_gene_term_matrix_bp,GO_gene_term_matrix_cc)
  
  ## Add the Ontology_group column:
  BP_metadata_matrix$Ontology_group="Biological_process"
  CC_metadata_matrix$Ontology_group="Cell_Components"
  ALL_metadata_matrix=rbind(BP_metadata_matrix,CC_metadata_matrix)
  #This way we have merged both terms
}


################
##  Pipeline  ##
################

# output_dir <- "Analyses/Outputs"
cluster_input_dir <- file.path("Analyses","Inputs","003_Clustering_tabs")

# Load lfc with clustering information
cluster_out <- read.table(file.path(cluster_input_dir,"lfc_data_with_clusters_information.txt"))
list_DE_infected_overtime <- readRDS(file.path(cluster_input_dir,"list_DE_infected_overtime.RDS"))
names_to_save <- names(list_DE_infected_overtime)

#We selected the k=16 cluster results
cluster_out1 <- cluster_out[,c(1:4,7)]
genes.names <- rownames(cluster_out1)### Gene names with Ensemble IDs (Novel Proteins)
list_per_clusters <- split(rownames(cluster_out1), cluster_out1[,5])
names(list_per_clusters) <- paste0("k=",names(list_per_clusters))

### Lets create clusters attending to the trend they have
### k= 1,4,7 Graduall Down
### k= 3,8,9 Graduall Up
### k= 10 Persistent Down
### k= 14 Persistent Down

# Define how to group the data frames
group_indices <- list(c(1,4,7), c(3,8,9),10,14,c(3,7),c(1,8), c(4,9),3,7,1,8,4,9)  # Each sublist represents a group of indices
group_names <- c("Gradual.Down","Gradual.Up","Persistent.Down", "Persistent.UP",
                 "Early","Mid","Late","k=3","k=7","k=1","k=8","k=4","k=9")
# Create a new list where each element is the combination of the selected groups
combined_list_by_trend <- lapply(group_indices, function(indices) {
  # Combine the selected data frames using rbind
  unique(unlist(list_per_clusters[indices], use.names = FALSE))
  
})

# Assign names to the elements of the combined list
names(combined_list_by_trend) <- group_names

# Save each element of the list as a separate file
# vec_to_loop <- seq_along(list_per_clusters)
vec_to_loop <- seq_along(combined_list_by_trend)
# vec_to_loop <- seq_along(combined_list_by_stage)
enrichment_dir <- "Analyses/Inputs/005_Enrichment"

for (i in vec_to_loop) {
  # Create a file name for each list element
  # file_name <- file.path(enrichment_dir,paste0("gene_names_k_",i, ".txt"))
  file_name <- file.path(enrichment_dir,paste0(group_names[i],".txt"))
  
  # Save the character vector to a file
  # write.table(list_per_clusters[[i]], file_name,col.names = F,row.names = F,
  #             quote = F)
  write.table(combined_list_by_trend[[i]], file_name,col.names = F,row.names = F,
  quote = F)
}

background_df <-read.table(file.path(input_dir,"002_Processed/ready_for_DE/reads.txt"))
background <- rownames(background_df)
length(background)

write.table(background,file.path(input_folder,"background_list.txt"))
### name that is going to be used on the bellow lists
# vector_of_names <- names(list_per_clusters)
vector_of_names <- names(combined_list_by_trend)

### Zero Configuration: To get a view of the stats
# config_dir <- "zero_Configuration_by_clusters"
config_dir <- "zero_Configuration_by_trend_and_stage"

{
  min_level <- 3
  max_level <- 8
  min_size  <- 1
  fdr       <- 0.05
  OR        <- 1.5
}


### 1st Configuration more permissive values
config_dir <- "1st_Configuration_by_trend_and_stage"

{
  min_level <- 8
  max_level <- 12
  min_size  <- 8
  fdr       <- 0.05
  OR        <- 2
}

### 2nd Configuration
config_dir <- "2nd_Configuration_by_trend_and_stage"

{
  min_level <- 7
  max_level <- 15
  min_size  <- 7
  fdr       <- 0.05
  OR        <- 2
}

### 3rd Configuration
config_dir <- "3rd_Configuration_by_trend_and_stage"
{
  min_level <- 1
  max_level <- 4
  min_size  <- 1
  fdr       <- 0.05
  OR        <- 2
}
### 4th Configuration
config_dir <- "4st_Configuration_by_trend_and_stage"
{
  min_level <- 1
  max_level <- 8
  min_size  <- 1
  fdr       <- 0.05
  OR        <- 2
}
# vec_to_loop <- seq_along(list_DE_infected_overtime)
# # vec_to_loop <- c(combined_list_by_trend)
# vec_to_loop <- c(12)

# config_dir <- "test_Configuration_by_clusters"
config_dir <- "6th_Configuration_by_trend_and_stage"

{
  min_level <- 3
  max_level <- 8
  min_size  <- 3
  fdr       <- 0.05
  OR        <- 2
}


for (i in vec_to_loop) {
  
  # type_of_enhanced <- readline("introduce the name of the list of genes your are using:")
  # min_level <- readline("introduce the minimum GO level you have to use:")
  # min_level <- as.integer(min_level)
  # max_level <- readline("introduce the maximum GO level you have to use:")
  # max_level <- as.integer(max_level)
  # min_size <- readline("introduce the minimum size a term needs to have to be taken into account:")
  # min_size <- as.integer(min_size)
  # fdr <- readline("introduce the threshold for fdr:")
  # fdr <- as.integer(fdr)
  
  # query_genes <- unlist(list_per_clusters[i])
  query_genes <- unlist(combined_list_by_trend[i])
  
  # print(paste(" I´m working with the condition of:",names(list_per_clusters[i])))
  
  print(paste(" I´m working with the condition of:",names(combined_list_by_trend[i])))
  print(paste("the length of the gene query is:",length(query_genes)))
  
  tryCatch(enrichment <- GO(query = query_genes,                           #Here we enter the query
                            background_matrix = GO_gene_term_matrix_all, #Here we load the matrix     
                            metadata_matrix = ALL_metadata_matrix,       #Here we load the metadata matrix
                            bg=background,                               #Here we load the background, THIS IS OPTIONAL, YOU CAN RUN IT WITHOUT BACKGROUND
                            min_GO_level = min_level,                            #Here you set the minimum GO level you have to use
                            max_GO_level = max_level,                            #Here you set the maximum GO level you have to use
                            test_mode = "two.sided",                    #Here you select the mode of the test
                            pvalue_correction_method = "qval",  #Here you set the pvalue correction method
                            min_term_size= min_size,                   #Here you set the minimum size a term needs to have to be taken into account
                            max_term_size = Inf,                #Here you set the maximum allowed term size
                            th_hit=fdr,                         #Here you select the threshold for fdr
                            th_hit_OR=OR,                        #Here you select the threshold for Odds Ratio
                            threshold_links = 0.5,              #And here you select  the threshold for connectivity, to create the networks
                            merge_duplicates=TRUE ),              #There are terms that have exactly the same genes. This option removes them and leaves only one
           error=function(e) e)
  if(inherits("Error in if (P > threshold) { : 
      valor ausente donde TRUE/FALSE es necesario", "error")){
    print("There was an error")
    next
  } 
  
  # Check if 'enrichment' exists and is not empty
  if (exists("enrichment") && !is.null(enrichment) && 
      is.list(enrichment) && nrow(enrichment$result) > 0) {
    res_enrichment <- enrichment$result
    
    # Create the directory for enrichment output
    enrichment_dir <- file.path(output_dir, "005_Enrichment_GO")
    dir.create(file.path(enrichment_dir,config_dir), recursive = TRUE, showWarnings = FALSE)
    
    # Create a subdirectory for the current element in the loop
    file_name_enrichment <- file.path(config_dir,vector_of_names[i])
    dir.create(file.path(enrichment_dir, file_name_enrichment), recursive = TRUE, showWarnings = FALSE)
    
    # Write the enrichment results to text and Excel files
    write.table(res_enrichment, file = file.path(enrichment_dir, file_name_enrichment, paste0("Enrichment_GO_",vector_of_names[i],".txt")))
    write.xlsx(res_enrichment, file = file.path(enrichment_dir, file_name_enrichment, paste0("Enrichment_GO_",vector_of_names[i],".xlsx")))
    
    # Filter the results for entries with non-NA 'community' values
    enrichment <- res_enrichment[which(!is.na(res_enrichment$community)), ]
    
    # Write the filtered results to text and Excel files
    write.table(enrichment, file = file.path(enrichment_dir, file_name_enrichment,paste0("Enrichment_GO_community_",vector_of_names[i],".txt")))
    write.xlsx(enrichment, file = file.path(enrichment_dir, file_name_enrichment, paste0("Enrichment_GO_community_",vector_of_names[i],".xlsx")))
   
    config_vector <- c(min_level=min_level,
                       max_level=max_level, 
                       min_size = min_size,
                       fdr=fdr,
                       OR=2)
    write.table(config_vector, file = file.path(enrichment_dir, file_name_enrichment,paste0("Configuration_used_",vector_of_names[i],".txt")))
    } 
  else {
    # Optional message to indicate that this element is being skipped
    message("Skipping ", vector_of_names[i], ": enrichment list is not defined or empty.")
    next # Skip to the next iteration of the loop
  }
}





# #You can set the terms with labels this way: label_indexes_enrichment=c(1,5,22,29,34,39,44,49,54:59)
# RNA_all_up_out=write_gephi_gexf(GO_object=RNA_all_up, #output GO
#                                 label_indexes=label_indexes_RNA_up,#numeric vector name i want in Gephi
#                                 output_file="Outputs/GO_enrichments/gephi_networks/RNA_up")
#And this way you can write your gephi network









