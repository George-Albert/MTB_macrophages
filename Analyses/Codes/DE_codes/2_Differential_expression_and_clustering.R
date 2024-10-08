###
### 0. Load dependencies & files
###

{
    options(width=10000)
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

    reads=read.table("inputs/processed/ready_for_DE/reads.txt")
    cols=read.table("inputs/processed/ready_for_DE/metadata.txt")
    cols_whole=read.table("inputs/processed/whole/metadata_whole.txt")
    length(which(colnames(reads)!=rownames(cols)))
    length(which(paste0(cols$Individual,"_",cols$Setup)!=rownames(cols)))
}

###
### 1. Run linear models.
###

{
    design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols)
    
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=FALSE)
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    
    betas=fit$coefficients[,103:109]
    ps=fit$p.value[,103:109]
    
    fdrs=ps
    
    for(i in 1:ncol(ps))
    {
        fdrs[,i]=p.adjust(ps[,i],method="BH")
        print(paste(colnames(ps)[i],length(which(fdrs[,i]<0.05 & abs(betas[,i])>0.2)),"hits."))
    }
    
    #[1] "TimePoint20h 5440 hits."
    #[1] "TimePoint48h 6354 hits."
    #[1] "TimePoint72h 6734 hits."
    #[1] "TimePoint02h:TreatmentTB 2270 hits."
    #[1] "TimePoint20h:TreatmentTB 5567 hits."
    #[1] "TimePoint48h:TreatmentTB 6378 hits."
    #[1] "TimePoint72h:TreatmentTB 6896 hits."
    
    system("mkdir -p Outputs/DE_def/tables")
    save(file="Outputs/DE_def/tables/ajuste.Rdata",fit)
    
    system("mkdir -p Outputs/DE_def/tables")
    system("mkdir -p Outputs/DE_def/clustering_plots")
    
    ## Define function to retrieve contrasts.
    
    get_contrast=function(fit,design,contrast,name,table,th=0.05,th_size=0.2){
        vec=rep(0,length(colnames(design)))
        vec[abs(contrast)]=contrast/abs(contrast)
        fit2 <- contrasts.fit(fit, vec)
        fit2 <- eBayes(fit2)
        
        tab=data.frame(beta=fit2$coefficients,p=fit2$p.value)
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
    
    write.table(results,"Outputs/DE_def/resultados.txt")

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


###
### 2. Do clustering. Version 1: Hierarchical clustering.
###

{
## Build table with the infection logFCs of genes that respond to infection at at least one timepoint.

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

pdf("Outputs/DE_def/clustering_plots/arbol1_mod.pdf",height=30,width=10)
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

pdf("Outputs/DE_def/clustering_plots/heatmap_exp_ok.pdf",width=5,height=30)
print(pl_tot)
dev.off()

### This is problematic: if we run HCs cutting the tree at arbitrarily defined cutoffs (or sleecting a certain number of clusters), we produce too many small clusters at some sections of the tree, or too large clusters at others... alternatively (which I did provisionally, but don´t quite like) one needs to define arbitrary cutoffs at different parts of the dendrogram to retrieve clusters that look ok, ít is too arbitrary/hard to justify...

}

###
### 3. Do clustering. Version 2: k-means
###

{
number_clusters=12
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

pdf("Outputs/DE_def/clustering_plots/kmneans_plot.pdf",width=7,height=30)
print(kmneans_plot)
dev.off()


}
