#########################
## 0.Load Dependencies ##
#########################

{library(tidyverse)
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
    library(reshape2)}
    
### Set dir
main_dir <- getwd()
setwd(main_dir)
input_dir  <- "Analyses/Inputs"
script_id  <- "002_Differential_expression"
output_dir <- file.path("Analyses/Outputs",script_id)
data_to_plot_dir <- file.path(output_dir,"002_Data_to_plot")

dir.create(data_to_plot_dir,recursive = T,showWarnings = F)
dir.create(file.path(input_dir),recursive = T,showWarnings = F)

### load data
reads=read.table(file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
cols=read.table(file.path(input_dir,"002_Processed","ready_for_DE","metadata.txt"))
cols_whole=read.table(file.path(input_dir,"002_Processed","whole","metadata_whole.txt"))

length(which(colnames(reads)!=rownames(cols)))
length(which(paste0(cols$Individual,"_",cols$Setup)!=rownames(cols)))

##############
##  Step #1 ##
##############

{
    #############################
    ##  1. Run linear models.  ##
    #############################
    design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols)

    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=T)
    
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    
    betas=fit$coefficients[,103:109]
    ps=fit$p.value[,103:109]
    SE <- as.data.frame(fit$stdev.unscaled*fit$sigma)
    SE <- SE[,103:109]
    
    fdrs=ps
    
    for(i in 1:ncol(ps)){
        fdrs[,i]=p.adjust(ps[,i],method="BH")
        print(paste(colnames(ps)[i],length(which(fdrs[,i]<0.05 & abs(betas[,i])>0.2)),"hits."))
    }
    
    # [1] "TimePoint20h 5425 hits."
    # [1] "TimePoint48h 6337 hits."
    # [1] "TimePoint72h 6720 hits."
    # [1] "TimePoint02h:TreatmentTB 2270 hits."
    # [1] "TimePoint20h:TreatmentTB 5550 hits."
    # [1] "TimePoint48h:TreatmentTB 6359 hits."
    # [1] "TimePoint72h:TreatmentTB 6875 hits."
    
    dir.create(file.path(output_dir,"003_DE_def","tables"),recursive=T,showWarnings = F)
    
    save(file=file.path(output_dir,"003_DE_def","tables","ajuste.Rdata"),fit)
    
    ### Define function to retrieve contrasts.
    
    get_contrast=function(fit,design,contrast,name,table,th=0.05,th_size=0.2){
        vec=rep(0,length(colnames(design)))
        vec[abs(contrast)]=contrast/abs(contrast)
        fit2 <- contrasts.fit(fit, vec)
        fit2 <- eBayes(fit2)
        
        # SE <- as.data.frame(fit$stdev.unscaled*fit$sigma)
        tab=data.frame(beta=fit2$coefficients,p=fit2$p.value,SE=fit2$stdev.unscaled*fit2$sigma,t=fit2$t)
        #tab$fdr=qvalue(tab$p)$qvalues
        tab$fdr=p.adjust(tab$p,method="BH")
        tab$hit=0
        tab$hit[which(tab$fdr<th & abs(tab$beta)>th_size)]=1
        print(paste0(name,": ",length(which(tab$hit==1))," hits"))
        colnames(tab)=paste0(colnames(tab),"_",name)
        
        if(missing(table)){
            return(tab)
        }else{
            tab=cbind(table,tab)
            return(tab)
        }
    }
    
    inf_2=106
    inf_20=107
    inf_48=108
    inf_72=109
    t_2_20_NI=103
    t_20_48_NI=c(104,-103)
    t_48_72_NI=c(105,-104)
    t_2_20_TB=c(103,107,-106)
    t_20_48_TB=c(104,-103,108,-107)
    t_48_72_TB=c(105,-104,109,-108)
    dinf_2_20_NI=c(107,-106)
    dinf_20_48_NI=c(108,-107)
    dinf_48_72_NI=c(109,-108)
    
    results=get_contrast(fit,design,inf_2,"inf_2")
    results=get_contrast(fit,design,inf_20,"inf_20",results)
    results=get_contrast(fit,design,inf_48,"inf_48",results)
    results=get_contrast(fit,design,inf_72,"inf_72",results)
    results=get_contrast(fit,design,t_2_20_NI,"t_2_20_NI",results)
    results=get_contrast(fit,design,t_20_48_NI,"t_20_48_NI",results)
    results=get_contrast(fit,design,t_48_72_NI,"t_48_72_NI",results)
    results=get_contrast(fit,design,t_2_20_TB,"t_2_20_TB",results)
    results=get_contrast(fit,design,t_20_48_TB,"t_20_48_TB",results)
    results=get_contrast(fit,design,t_48_72_TB,"t_48_72_TB",results)
    results=get_contrast(fit,design,dinf_2_20_NI,"dinf_2_20_NI",results)
    results=get_contrast(fit,design,dinf_20_48_NI,"dinf_20_48_NI",results)
    results=get_contrast(fit,design,dinf_48_72_NI,"dinf_48_72_NI",results)
    results=data.frame(results)

    write.table(results,file.path(output_dir,"003_DE_def","resultados.txt"))
    
    x <- results
    
    grep("hit",colnames(x))
    colnames(x)[grep("hit",colnames(x))]
    
    x <- x[,colnames(x)[grep("hit",colnames(x))]]
    
    DE.genes <- data.frame(colSums(x == 1, na.rm = TRUE))
    colnames(DE.genes) <- "DE.Genes"
    
    DE_tab <- data.frame(Condition=rownames(DE.genes),DE.genes=DE.genes$DE.Genes)
    write.table(DE_tab,file.path(output_dir,"003_DE_def","DE_table.txt"))
    write.xlsx(DE_tab,file.path(output_dir,"003_DE_def","DE_table.xlsx"))
    
    ### Intercept is TimePoint2h
    # Condition DE.genes
    # 1          hit_inf_2     2270
    # 2         hit_inf_20     5550
    # 3         hit_inf_48     6359
    # 4         hit_inf_72     6875
    # 5      hit_t_2_20_NI     5425
    # 6     hit_t_20_48_NI     2456
    # 7     hit_t_48_72_NI      687
    # 8      hit_t_2_20_TB     5944
    # 9     hit_t_20_48_TB     1902
    # 10    hit_t_48_72_TB     1010
    # 11  hit_dinf_2_20_NI     4344
    # 12 hit_dinf_20_48_NI     1759
    # 13 hit_dinf_48_72_NI      407
   
     ### Joaquin previous results 
    #[1] "inf_2: 2270 hits"
    #[1] "inf_20: 5567 hits"
    #[1] "inf_48: 6378 hits"
    #[1] "inf_72: 6896 hits"
    #[1] "t_2_20_NI: 5440 hits"
    #[1] "t_20_48_NI: 2469 hits"
    #[1] "t_48_72_NI: 693 hits"
    #[1] "t_2_20_TB: 5980 hits"
    #[1] "t_20_48_TB: 1909 hits"
    #[1] "t_48_72_TB: 1012 hits"
    #[1] "dinf_2_20_NI: 4363 hits"
    #[1] "dinf_20_48_NI: 1758 hits"
    #[1] "dinf_48_72_NI: 411 hits"
}

