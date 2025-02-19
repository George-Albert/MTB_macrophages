
#########################
##  Load Dependencies  ##
#########################

{
    
    library(ggplot2)
    library(ggrepel)
    library(limma)
    library(edgeR)
    library(qvalue)
    library(rtracklayer)
    library(RColorBrewer)
    library(reshape2)
    library(cowplot)
    library(readr)
    library(Rtsne)
    library(HGNChelper)
     
}


###############################
##  Set labeling of samples  ##
###############################



main_dir <- setwd(getwd())
input_dir <- file.path("Analyses/Inputs")
output_dir <- file.path("Analyses/Outputs")



{
    meta_dir <- file.path(input_dir,"001_Raw_data","metadata.txt")
    cols=read.table(meta_dir,header=TRUE)
    # Create Setup columns with the union of Condition and Time Point
    cols$Setup=paste0(cols$Condition,"_",cols$TimePoint)
    setups=unique(cols$Setup)
    
    reads_dir <- file.path(input_dir,"001_Raw_data","reads.txt")
    reads=read.table(reads_dir,header=TRUE,check.names=FALSE)
    rownames(reads)=reads[,1]
    reads=reads[,2:ncol(reads)]
    
    ##########################################################################
    ##                 Hay algunas samples con el SampleID                  ##
    ##                   duplicado en la tabla de metadatos                 ##
    ##  ,en las colnames de reads se diferencian con un .1. Asumo que esto  ##
    ##                 no es problema para saber cuál de las                ##
    ##    columnas se corresponde a un duplicado porque están en el mismo   ##
    ##               orden en data (columnas) y metadata (rows)             ##
    ##########################################################################
    
    set_differ=which(cols$SampleID!=colnames(reads))
    cols$SampleID=as.character(cols$SampleID)
    cols$Duplicated=0
    dupes=cols$SampleID[which(duplicated(cols$SampleID))]
    cols$Duplicated[which(cols$SampleID %in% dupes)]=1
    cols$Duplicated[which(cols$SampleID!=colnames(reads))]=2
    
    #######################################################################
    ##       Así que añado a mano los .1 también en los metadatos,       ##
    ##   en una nueva columna (que almacenará los sampleIDs corregidos.  ##
    #######################################################################
    
    cols$Sample=cols$SampleID
    cols$Sample[set_differ]=paste0(cols$SampleID[set_differ],".1")
    
    ### This is to see that all these dups appeared only two times 
    length(unique(cols$Sample))
    ## 883
    nrow(cols)
    ## 883
    
    ## This is provisional, rownames will be re declared later
    rownames(cols)=cols$Sample
    
    length(which(rownames(cols)!=colnames(reads)))
    length(which(cols$Sample!=colnames(reads)))
    
    ##############################################################
    ##       AHORA CORRIJO LOS .1 Y LOS PONGO EN LAS MUESTRAS   ##
    ##       CON DUPLICATED==1, QUE TERMINARAN POR NO USARSE    ##
    ##############################################################
    
    cols$Sample[which(cols$Duplicated==1)]=paste0(cols$SampleID[which(cols$Duplicated==1)],".1")
    cols$Sample[which(cols$Duplicated==2)]=cols$SampleID[which(cols$Duplicated==2)]
    length(unique(cols$Sample))
    
    rownames(cols)=cols$Sample
    colnames(reads)=cols$Sample

    
    cols=cols[order(rownames(cols)),]
    reads=reads[,order(colnames(reads))]
    
    length(which(colnames(reads)!=rownames(cols)))

    ### Conclusion: labeling OK.
}

###################################################################################
##            2. Check sample depths (specially in Duplicated samples)           ##
##           & Check for flowcell merging (i.e. are the repeated samples         ##
##   new sequencing instances, or the merging between the old and the new one?)  ##
###################################################################################

{
    sumas=apply(reads,2,sum)
    
    cols$Original_depth=sumas
    cols$QC_labels=""
    
    ### check the distribution
    tabla=cols[order(-cols$Duplicated),]
    pl=ggplot()+geom_jitter(aes(x=Flowcell,y=Original_depth,
                                color=as.character(Duplicated)),
                            data=tabla)+coord_flip()
    
    ##################################################################
    ##    From visual inspection, these are samples that will be    ##
    ##           good to track because of different reasons         ##
    ##            (To do that I define the variable QC_label        ##
    ##################################################################
    
    ## Deeper outliers than their Flowcell averages
    set_1=which((cols$Flowcell %in% c("171211_K00243_0104_BHLVWHBBXX",
                                      "171124_K00243_0098_AHLNGVBBXX") & 
                     cols$Original_depth>3E7) | cols$Original_depth>3.5E7)
    set_2=which(cols$Flowcell!="170928_D00416_0387_ACBL36ANXX" & cols$Original_depth<1.5E7)
    set_3=which(cols$Flowcell!="171211_K00243_0104_BHLVWHBBXX" & cols$Duplicated==2)

    samples_to_label=c(set_1,set_2,set_3)
    cols$QC_labels[samples_to_label]=cols$SampleID[samples_to_label]
    
    cols$Notes=""
    cols$Notes[set_1]="Depth_outlier_above"
    cols$Notes[set_2]="Depth_outlier_below"
    cols$Notes[set_3]="Resequenced_in_minoritary_batch"

    pl1=ggplot(cols)+
    geom_jitter(aes(x=Flowcell,y=Original_depth,color=as.character(Duplicated)),
                position = position_jitter(seed = 1))+
    geom_text_repel(aes(x=Flowcell,y=Original_depth,label=QC_labels),
                    position = position_jitter(seed = 1),size=5)+
    coord_flip()+
        theme_linedraw()+
        theme(
            axis.text = element_text(size = 12),  # Aumenta el tamaño del texto de los ejes
            axis.title = element_text(size = 14),  # Aumenta el tamaño de los títulos de los ejes
            legend.text = element_text(size = 12),  # Aumenta el tamaño del texto de la leyenda
            legend.title = element_text(size = 14),  # Aumenta el tamaño del título de la leyenda
            plot.title = element_text(size = 16)  # Aumenta el tamaño del título del gráfico
        )
    
    pl1
    
    dir.create(file.path(output_dir,"001_Figures","001_QC"),recursive = T)
    QC_dir <- file.path(output_dir,"001_Figures","001_QC")
    
    pdf(file.path(QC_dir,"001_Depths_per_FC.pdf"),width=12,height=6)
    print(pl1)
    dev.off()
    
    ###################################################################
    ##                    Two more sanity checks:                    ##
    ##           Check that all samples that are not unique          ##
    ##   form pairs (i.e. no sample was sequenced more than twice=)  ##
    ###################################################################
    
    # Remember that, for each sample ID that appear more than once, 
    # Duplicated=1 is the first instance and =2 all the others).
    length(which(cols$Duplicated==1))
    length(which(cols$Duplicated==2))
    
    ## So, it is ok.
     
    ## Another sanity check: in all the samples that are repeated, 
    ## the repetitions are resequenced instances, and one is not a merging of 
    ## two sequencing instances instead.
    
    unos=reads[,which(cols$Duplicated==1)]
    doses=reads[,which(cols$Duplicated==2)]
    test=data.frame(colnames(unos),colnames(doses))
    
    # Check that the pairing is correct
    test$check=paste0(test[,1],".1")
    length(which(test[,2]!=test[,3]))
     
    difs=doses-unos
    mins=apply(difs,2,min)
    maxs=apply(difs,2,max)
    
    summary(mins)
    summary(maxs)
       
    # All samples have some entries that are higher than the resequenced pair, 
    # and some that are lower: none is a merging.
    
    #####################################################################################
    ##                                   Conclusion:                                   ##
    ##       There was a funky flowcell, and a couple of other bad samples that,       ##
    ##             once repeated, yielded good samples with two exceptions             ##
    ##       (C50_0_72 & C60_TB_72). A third sample was also bad, .-but was not        ##
    ##              resequenced (C12_0_2h). As for the rest, some of them              ##
    ##   could be eventually removed if we were to get picky regarding library depths  ##
    ##                distributions, but I don´t think that is necesary.               ##
    #####################################################################################
    
    ### For now, let's take the 1E7 QC cutoff here.
    cols$Depth_criterium=0
    cols$Depth_criterium[which(cols$Original_depth>1E7)]=1
    
}

###############################################################################################
##  3. First sample subsetting (according to Depth_criterium, declare _safe objects first).  ##
###############################################################################################

{
    cols_whole=cols
    reads_whole=reads
    
    cols=cols[which(cols$Depth_criterium==1),]
    reads=reads[,which(colnames(reads) %in% rownames(cols))]
    
    tab=cols
    tab$x=1
    tab$x_jit=1+rnorm(nrow(cols),sd=0.1)
    
    pl_Depths_distribution=ggplot(tab)+
        geom_jitter(aes(x=x_jit,y=Original_depth,color=as.factor(Duplicated)))+
        geom_violin(aes(x=x,y=Original_depth),fill=NA,linewidth=1.5)+
    geom_hline(yintercept=5000000,color="grey90")+
    geom_hline(yintercept=10000000,color="grey80")+
    geom_hline(yintercept=40000000,color="grey60")+
    geom_hline(yintercept=15000000,color="grey40")+
    geom_hline(yintercept=35000000,color="black")+
        theme_classic()+
        theme(
            axis.text = element_text(size = 12),  # Aumenta el tamaño del texto de los ejes
            axis.title = element_text(size = 14),  # Aumenta el tamaño de los títulos de los ejes
            legend.text = element_text(size = 12),  # Aumenta el tamaño del texto de la leyenda
            legend.title = element_text(size = 14),  # Aumenta el tamaño del título de la leyenda
            plot.title = element_text(size = 16)  # Aumenta el tamaño del título del gráfico
        )
    pl_Depths_distribution
    
    pdf(file.path(QC_dir,"002_Global_depth_distribution.pdf"),width=7.5,height=6)
    print(pl_Depths_distribution)
    dev.off()
}

#########################################################################################
##      4. Basic exploration on how are the individuals and flowcells distributed      ##
##   across conditions and times after that first subsetting to identify pity samples  ##
#########################################################################################

{
    tab=data.frame(ind=unique(cols$IndividualID))
    tab$mult=0
    
    for(i in c(1:nrow(tab)))
    {
        tab$mult[i]=length(which(cols$IndividualID %in% tab$ind[i]))
    }
    
    tab=tab[order(tab$mult),]
    tab$ah=0
    tab$b=0
    tab$c=0
    tab$d=0
    tab$TB_20h=0
    tab$TB_2h=0
    tab$TB_48h=0
    tab$TB_72h=0
    tab$e=0
    tab$TB_18h=0
    colnames(tab)=c("ind","mult",as.character(setups))
    
    for(i in c(1:length(setups)))
    {
        given_setup=setups[i]
        
        ## These are the individuals having at least one sample in the given_setup
        individuals_present=unique(cols$IndividualID[which(cols$Setup %in% given_setup)])
        ## I put 1s in the rows corresponding to them
        rows=which(tab$ind %in% individuals_present)
        tab[rows,i+2]=1
    }
    
    ## This is to collapse the Setup columns populated with 0s and 1s into a 
    ## binary string of length(setusps) bits
    
    colapsar=function(x){return(paste(x,collapse=""))}
    
    obj=tab[,c(3:12)]
    tab$config=apply(obj,1,colapsar)
    tab=tab[order(tab$mult,tab$config),]
    
    tab$fcs=0
    
    for(i in 1:nrow(tab))
    {
        given_individual=tab$ind[i]
        tabla=cols[which(cols$IndividualID==given_individual),]
        valor=length(unique(tabla$Flowcell))
        print(paste(i,valor,"\n"))
        if(valor>0)
        tab$fcs[i]=valor
    }
    
    length(which(tab$fcs>1))
    
    ###############################################################################
    ##  5. These ´pity´ guys are the individuals whose samples got sequenced in  ##
    ##       more than one flowcell. This is a pity, bc if Individuals were      ##
    ##     perfectly nested within flowcells that would have simplified a lot    ##
    ##                             the model designs).                           ##
    ###############################################################################
    
    pity=tab[(which(tab$fcs>1)),]
    pity_cols=cols[which(cols$IndividualID %in% pity$ind),]
    
    ### Exploring this pity_cols table, I identify visually the 14 samples that 
    ### one should remove (that's the minimum amount) to recover a perfect 
    ### nesting of individuals within flowcells. Alternatively, we could use 
    ### the variable intercept to model fixed-effect intercepts (it's almost
    ### always the individual, except for these 6 pity guys for which an 
    ### individual is splitted in 2 flowcells (we would fit 2 intercepts for 
    ### each of them).
    
    pity_samples=c("10_0_2h","10_TB_2h","C053_0_72h","C13_0_48h","C13_0_72h",
                   "C16_0_72h","C23_0_18h","C23_0_2h","C23_TB_18h","C23_TB_2h",
                   "C9_0_18h","C9_0_2h","C9_TB_18h","C9_TB_2h")
    
    cols$Ind_FC_nesting_criterium=1
    cols$Ind_FC_nesting_criterium[which(rownames(cols) %in% pity_samples)]=0
    cols$Notes[which(rownames(cols) %in% pity_samples)]=
        paste0(cols$Notes[which(rownames(cols) %in% pity_samples)],
               "_Breaks_Ind_in_FC_nesting")
    
    cols_whole$Ind_FC_nesting_criterium=1
    cols_whole$Ind_FC_nesting_criterium[which(rownames(cols_whole) %in% pity_samples)]=0
    cols_whole$Notes[which(rownames(cols_whole) %in% pity_samples)]=
        paste0(cols_whole$Notes[which(rownames(cols_whole) %in% pity_samples)],
               "_Breaks_Ind_in_FC_nesting")
    
    ## From C23 and C9 I label as pity the early ones because I lack late 
    ## samples for more individuals...
    
    tab_ind=tab
    
    tab=data.frame(fc=unique(cols$Flowcell))
    tab$mult=0
    
    for(i in c(1:nrow(tab)))
    {
        tab$mult[i]=length(which(cols$Flowcell %in% tab$fc[i]))
    }
    
    tab=tab[order(tab$mult),]
    tab$ah=0
    tab$b=0
    tab$c=0
    tab$d=0
    tab$TB_20h=0
    tab$TB_2h=0
    tab$TB_48h=0
    tab$TB_72h=0
    tab$e=0
    tab$TB_18h=0
    
    colnames(tab)=c("fc","mult",as.character(setups))
    
    for(i in c(1:length(setups)))
    {
        given_setup=setups[i]
        
        for(j in c(1:nrow(tab)))
        {
            given_fc=tab$fc[j]
            tab[j,i+2]=length(which(cols$Flowcell==given_fc & cols$Setup==given_setup))
        }
    }
    
    dir.create(file.path(QC_dir,"Tables"),recursive=T,showWarnings = F)
    write.table(tab, file.path(QC_dir,"Tables","Distribution_of_flowcells_across_conditions_upon_filtering_one.txt"))
    write.table(tab_ind,file.path(QC_dir,"Tables","Distribution_of_individuals_across_conditions_upon_filtering_one.txt"))
    
    tabm=tab[,c(1,3:ncol(tab))]
    tabm=melt(tabm,ID="fc")
    
    pl_tiles=ggplot(tabm, aes(variable,fc)) +
        geom_tile(aes(fill = value),color="black") +
        geom_text(aes(label = round(value, 1))) +
        scale_fill_gradient(low = "white", high = "red")+
        theme_bw()+
        theme(panel.border = element_blank(),  # Quita el borde del panel
              panel.grid = element_blank(),
              axis.text = element_text(size = 12),  # Aumenta el tamaño del texto de los ejes
              axis.title = element_text(size = 14),  # Aumenta el tamaño de los títulos de los ejes
              legend.text = element_text(size = 12),  # Aumenta el tamaño del texto de la leyenda
              legend.title = element_text(size = 14),  # Aumenta el tamaño del título de la leyenda
              plot.title = element_text(size = 16) ) +  # Quita las líneas de la cuadrícula
        scale_x_discrete(expand = c(0, 0)) +  # Quita el espacio en blanco en el eje x
        scale_y_discrete(expand = c(0, 0))   # Quita el espacio en blanco en el eje y
        
    
    pl_tiles
    
    tabm=tab_ind[,c(1,3:12)]
    tabm=melt(tabm,ID="ind")
       
   pl_tiles2=ggplot(tabm, aes(variable,ind)) +
       geom_tile(aes(fill = value),color="black") +
       scale_fill_gradient(low = "white", high = "red")+
       ylab("Individual")+xlab("Setup")+
       theme(axis.text.y   = element_blank(),
             legend.position="left",
             panel.border = element_blank(),  # Quita el borde del panel
             panel.grid = element_blank(),   # Quita las líneas de la cuadrícula
             axis.text = element_text(size = 12),  # Aumenta el tamaño del texto de los ejes
             axis.title = element_text(size = 14),  # Aumenta el tamaño de los títulos de los ejes
             legend.text = element_text(size = 12),  # Aumenta el tamaño del texto de la leyenda
             legend.title = element_text(size = 14),  # Aumenta el tamaño del título de la leyenda
             plot.title = element_text(size = 16))+ 
       scale_x_discrete(expand = c(0, 0)) +  # Quita el espacio en blanco en el eje x
       scale_y_discrete(expand = c(0, 0))    # Quita el espacio en blanco en el eje y
   
   pl_tiles2  
   
       
   pl_bars=ggplot(tab_ind)+
       geom_point(aes(x=ind,y=as.character(fcs)))+
       theme_classic()+
       theme(axis.title.y= element_blank(),
             legend.position="left")+
       coord_flip()+ylab("# of FCs")
   
   pl_tot=plot_grid(pl_tiles2,pl_bars,ncol=2,rel_widths=c(5,1))
   
   output_dir
   
   pdf(file.path(QC_dir,"003_Setup_distribution_across_flowcells.pdf"),width=10,height=3)
   print(pl_tiles)
   dev.off()
     
   pdf(file.path(QC_dir,"004_Setup_distribution_across_individuals.pdf"),width=8,height=12)
   print(pl_tot)
   dev.off()
    
}

###################################################################################
##  6. Second filtering: pity samples out: now individual is nested in flowcell  ##
###################################################################################

{
    set=which(rownames(cols) %in% pity_samples)
    cols=cols[-set,]
    reads=reads[,which(colnames(reads) %in% rownames(cols))]
}

################################################################################
##   Low expression filtering first on those 2 filtering (depth outliers and  ##
##                          FC_ind_nesting breakers).                         ##
################################################################################

{
    design=model.matrix(~1,data=cols)
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=TRUE)
    
    exp=data.frame(v$E)
    medians=exp[,1:10]
    
    
    colnames(medians)=setups
    
    for(i in 1:ncol(medians))
    {
        case=colnames(medians)[i]
        set=which(cols$Setup %in% case)
        print(paste(i,colnames(medians)[i],length(set)))
        tab=exp[,set]
        medians[,i]=apply(tab,1,median)
    }
    
    medians$max_all=apply(medians,1,max)
    medians$max_used=apply(medians[,c(2,5,6,7,8,10)],1,max)
    
    th=1
    length(which(medians$max_all>th))
    ## 12989
    th=2
    length(which(medians$max_all>th))
    ## 11548
    length(which(medians$max_all>th & medians$max_used<th))
    ## 262
    length(which(medians$max_all>th & medians$max_used<1))
    ## 13
    
    ### Two options to do this:
    ### Option A: I select genes that are expressed at a threshold of
    ### log2(cpm)=1 in at least one of the conditions of interest 
    ### (the first control and all the infected conditions) (I got 12989 genes)
    ### Option B: for considering a gene, I'll ask it to have at least a
    ### median(log2(cpm))>2 in any of the ten combinations condition-time. 
    ### At that threshold, two things happen:
    ### 1-The # of genes that pass is still nice (11548) and doesn't drop much 
    ### from the maybe more usual threshold of log(cpm)=1 (12989).
    ### One can disregard the dychotomy of whether using all the 10 setups to 
    ### do this filtering, only the relevant (or even fancier, model all genes,
    ### remove the effects of time in the control as we'll do later, and do
    ### the filtering on the relevant conditions, once the time has been 
    ### regressed out, model based). This is because, of the 262 genes that
    ### pass the threshold by considering all conditions which would have
    ### failed if considering only the relevant ones, only 13 would have 
    ### failed also, considering the relevant ones, if using a median 
    ### log(cpm)=1 as a threshold (so most of them are still fairly expressed 
    ### in at least one of the interesting conditions too.
    
    ## So I simplify and do this:
    
    reads=reads[which(medians$max_all>th),]

}

#################################################################################
##      7. Let's groom a bit the files to make things easier from here on.     ##
##     Creating new columns (Batch, Individual, Treatment, Time and, Sample)   ##
##   with the same info of what it is in there already, but clean, for common  ##
##                                  use later.                                 ##
#################################################################################

{
    whole=cols_whole
    setted=cols
    ###
    cols=cols_whole
    
    cols$Flowcell=factor(cols$Flowcell)
    cols$Batch=""
    lets=letters[1:length(levels(cols$Flowcell))]
    
    for(i in 1:length(levels(cols$Flowcell)))
    {
        fl=levels(cols$Flowcell)[i]
        set=which(cols$Flowcell %in% fl)
        cols$Batch[set]=lets[i]
    }
    cols$Batch=factor(cols$Batch,levels=lets)
    
    cols$IndividualID=factor(cols$IndividualID)
    cols$Individual=paste0("I_",cols$IndividualID)
    cols$Individual=factor(cols$Individual)
    
    cols$Treatment=as.character(cols$Condition)
    cols$Treatment[which(cols$Condition==0)]="NI"
    cols$Treatment=factor(cols$Treatment)
    cols$Condition=factor(cols$Condition)
    
    cols$TimePoint=as.character(cols$TimePoint)
    cols$TimePoint[which(cols$TimePoint=="2h")]="02h"
    cols$Time=substr(cols$TimePoint,1,nchar(cols$TimePoint)-1)
    cols$Time=as.numeric(cols$Time)
    cols$TimePoint=factor(cols$TimePoint,levels=c("02h","18h","20h","48h","72h"))
    
    cols$Setup=paste0(cols$Treatment,"_",cols$TimePoint)
    length(which(cols$Sample!=colnames(reads_whole)))

    cols$Sample=paste0(cols$Individual,"_",cols$Setup)
    cols$Sample[which(cols$Duplicated==1)]=paste0(cols$Sample[which(cols$Duplicated==1)],".1")
    rownames(cols)=cols$Sample
    colnames(reads_whole)=rownames(cols)
    cols$Setup=factor(cols$Setup)
    cols$Intercept=paste0(cols$Individual,"_",cols$Batch)
    cols$Intercept=factor(cols$Intercept)


    whole_groomed=cols
    setted_groomed=whole_groomed[which(whole_groomed$Depth_criterium==1 & whole_groomed$Ind_FC_nesting_criterium==1),]
    cols=setted_groomed
    
    
    cols$Flowcell=factor(cols$Flowcell)
    cols$Batch=factor(cols$Batch)
    cols$IndividualID=factor(cols$IndividualID)
    cols$Individual=factor(cols$Individual)
    cols$Setup=factor(cols$Setup)
    cols$Intercept=factor(cols$Intercept)

    cols$Sample=cols$SampleID
    length(which(cols$Sample!=colnames(reads)))

    cols$Sample=paste0(cols$Individual,"_",cols$Setup)
    rownames(cols)=cols$Sample
    colnames(reads)=rownames(cols)
    
    length(which(rownames(cols) %in% rownames(whole_groomed)))
}

#######################################################################################
##  8. PCA, first attempt useful to identify misslabels, stim. failuers & outliers.  ##
#######################################################################################

{
    dir.create(file.path(output_dir,"001_Figures","002_PCA"),showWarnings = F)
    
    PCA_dir <- file.path(output_dir,"001_Figures","002_PCA")
    
    get_PCA=function(setups,name,type="infection"){
        
        colores_ref=brewer.pal(10,"RdYlBu")[c(6:10,5:1)]
        colores=colores_ref[which(levels(factor(cols$Setup)) %in% setups)]
        
        dir.create(file.path(PCA_dir,"clean","plots"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"clean","summaries"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"clean","data"),recursive = T,showWarnings = F)

        dir.create(file.path(PCA_dir,"dirty","plots"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"dirty","summaries"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"dirty","data"),recursive = T,showWarnings = F)
        
        set=which(cols$Setup %in% setups)
        cols_chunk=cols[set,]
        reads_chunk=reads[,set]
        cols_chunk$Individual=factor(cols_chunk$Individual)
        cols_chunk$Treatment=factor(cols_chunk$Treatment)
        cols_chunk$TimePoint=factor(cols_chunk$TimePoint)
        
        if(type=="infection"){design=model.matrix(~Individual+Treatment,data=cols_chunk)}
        if(type=="time"){design=model.matrix(~Individual+TimePoint,data=cols_chunk)}
        if(type=="all"){design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols_chunk)}
        
        dge <- DGEList(counts=reads_chunk)
        dge <- calcNormFactors(dge)
        v <- voom(dge,design,plot=FALSE)
        fit <-lmFit(v,design)
        fit <- eBayes(fit)
        model=fit$coefficients
        
        exp_dirty=v$E
        exp_clean=exp_dirty
        
        for(i in 2:length(unique(cols_chunk$Individual)))
        exp_clean=exp_clean-model[,i]%*%t(design[,i])
        
        ### First process dirty
        
        corr_mat <- cor(exp_dirty,method="pearson")
        pca <-prcomp(corr_mat,scale=T)
        sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
        
        datos <- data.frame(pca$x)
        colnames(datos)=paste0("PC",c(1:ncol(datos)))
        rownames(datos)=colnames(exp_dirty)
        datos=merge(datos,cols[colnames(exp_dirty),],by=0)
        
        pl=ggplot(datos) +
        geom_point(data=datos,aes(x=PC1, y=PC2,colour=Setup))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1,linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines'))+
            xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
            ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
        
        pdf(paste0(PCA_dir,"/dirty/plots/",name,".pdf"),width=6,height=4)
        print(pl)
        dev.off()
        
        sink(paste0(PCA_dir,"/dirty/summaries/",name,".txt"))
        print(sum_pca)
        sink()
        
        write.table(datos,paste0(PCA_dir,"/dirty/data/",name,".txt"))
        
        datos_dirty=datos
        ### Then process clean
        
        corr_mat <- cor(exp_clean,method="pearson")
        pca <-prcomp(corr_mat,scale=T)
        sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
        
        datos <- data.frame(pca$x)
        colnames(datos)=paste0("PC",c(1:ncol(datos)))
        rownames(datos)=colnames(exp_clean)
        datos=merge(datos,cols[colnames(exp_clean),],by=0)
        
        pl=ggplot(datos) +
        geom_point(data=datos,aes(x=PC1, y=PC2,colour=Setup))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1,linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines'))+
            xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
            ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
        
        pdf(paste0(PCA_dir,"/clean/plots/",name,".pdf"),width=6,height=4)
        print(pl)
        dev.off()
        
        sink(paste0(PCA_dir,"/clean/summaries/",name,".txt"))
        print(sum_pca)
        sink()
        
        write.table(datos,paste0(PCA_dir,"/clean/data/",name,".txt"))
        
        datos_clean=datos
        
        return(list(datos_dirty=datos_dirty,datos_clean=datos_clean,sum_pca=sum_pca))
    }
    
    pca_Infection_2h=get_PCA(setups=c("NI_02h","TB_02h"),name="Infection_2h",type="infection")
    pca_Infection_18h=get_PCA(setups=c("NI_18h","TB_18h"),name="Infection_18h",type="infection")
    pca_Infection_20h=get_PCA(setups=c("NI_20h","TB_20h"),name="Infection_20h",type="infection")
    pca_Infection_48h=get_PCA(setups=c("NI_48h","TB_48h"),name="Infection_48h",type="infection")
    pca_Infection_72h=get_PCA(setups=c("NI_72h","TB_72h"),name="Infection_72h",type="infection")
    
    pca_Time_2_18_20_NI=get_PCA(setups=c("NI_02h","NI_18h","NI_20h"),name="Time_2_18_20_NI",type="time")
    pca_Time_18_20_48_NI=get_PCA(setups=c("NI_18h","NI_20h","NI_48h"),name="Time_18_20_48_NI",type="time")
    pca_Time_48_72_NI=get_PCA(setups=c("NI_48h","NI_72h"),name="Time_48_72_NI",type="time")
    
    pca_Time_2_18_20_TB=get_PCA(setups=c("TB_02h","TB_18h","TB_20h"),name="Time_2_18_20_TB",type="time")
    pca_Time_18_20_48_TB=get_PCA(setups=c("TB_18h","TB_20h","TB_48h"),name="Time_18_20_48_TB",type="time")
    pca_Time_48_72_TB=get_PCA(setups=c("TB_48h","TB_72h"),name="Time_48_72_TB",type="time")
    
    pca_Time_all_NI=get_PCA(setups=c("NI_02h","NI_18h","NI_20h","NI_48h","NI_72h"),name="Time_all_NI",type="time")
    pca_Time_all_TB=get_PCA(setups=c("TB_02h","TB_18h","TB_20h","TB_48h","TB_72h"),name="Time_all_TB",type="time")
    pca_Time_all=get_PCA(setups=levels(factor(cols$Setup)),name="Time_all",type="all")
    
}

##########################################################################
##  9. Here I look the plots and spot the samples that I need to label  ##
##########################################################################

{
    put_labels=function(tab,ranges){
        tab$label=""
        rownames(tab)=tab$Row.names
        for(i in 1:nrow(ranges))
        {
            set=which(
            (tab$Setup %in% ranges$setups[i])
            &
            (
            (tab$PC1>ranges$minima_PC1[i] & tab$PC1<ranges$maxima_PC1[i])
            |
            (tab$PC2>ranges$minima_PC2[i] & tab$PC2<ranges$maxima_PC2[i])
            )
            )
            tab$label[set]=rownames(tab)[set]
        }
        return(tab)
    }
    
    label_PCA=function(data,ranges,name){
        
        dir.create(file.path(PCA_dir,"clean","plots_labelled"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"clean","data_labelled"),recursive = T,showWarnings = F)
        
        colores_ref=brewer.pal(10,"RdYlBu")[c(6:10,5:1)]
        colores=colores_ref[which(levels(factor(cols$Setup)) %in% ranges$setups)]
        
        datos_dirty=data[[1]]
        datos_clean=data[[2]]
        sum_pca=data[[3]]
        
        datos_clean=put_labels(datos_clean,ranges)
        
        pl=ggplot(datos_clean) +
        geom_point(aes(x=PC1, y=PC2,colour=Setup))+
        geom_text_repel(aes(x=PC1, y=PC2,label=label))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,
                                    linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines'))+
            xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
            ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
        
        pdf(paste0(PCA_dir,"/clean/plots_labelled/",name,".pdf"),width=8,height=6)
        print(pl)
        dev.off()
        
        write.table(datos_clean,paste0(PCA_dir,"/clean/data_labelled/",name,".txt"))
        return(list(datos_dirty=datos_dirty,datos_clean=datos_clean))
    }
    
    ranges_df=data.frame(
    setups=c("NI_02h","TB_02h"),
    minima_PC1=c(0,-Inf),
    maxima_PC1=c(Inf,4),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Infection_2h=label_PCA(data=pca_Infection_2h,ranges=ranges_df,name="Infection_2h")
    
    ranges_df=data.frame(
    setups=c("NI_18h","TB_18h"),
    minima_PC1=c(2,-Inf),
    maxima_PC1=c(Inf,3),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Infection_18h=label_PCA(data=pca_Infection_18h,ranges=ranges_df,name="Infection_18h")
    
    ranges_df=data.frame(
    setups=c("NI_20h","TB_20h"),
    minima_PC1=c(0,-Inf),
    maxima_PC1=c(Inf,0),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Infection_20h=label_PCA(data=pca_Infection_20h,ranges=ranges_df,name="Infection_20h")
    
    ranges_df=data.frame(
    setups=c("NI_48h","TB_48h"),
    minima_PC1=c(0,-Inf),
    maxima_PC1=c(Inf,5),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Infection_48h=label_PCA(data=pca_Infection_48h,ranges=ranges_df,name="Infection_48h")
    
    ranges_df=data.frame(
    setups=c("NI_72h","TB_72h"),
    minima_PC1=c(0,-Inf),
    maxima_PC1=c(Inf,5),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Infection_72h=label_PCA(data=pca_Infection_72h,ranges=ranges_df,name="Infection_72h")
    
    ranges_df=data.frame(
    setups=c("NI_02h","NI_18h","NI_20h"),
    minima_PC1=c(Inf,Inf,Inf),
    maxima_PC1=c(Inf,Inf,Inf),
    minima_PC2=c(Inf,Inf,-Inf),
    maxima_PC2=c(Inf,Inf,5))
    
    pca_Time_2_18_20_NI=label_PCA(data=pca_Time_2_18_20_NI,ranges=ranges_df,name="Time_2_18_20_NI")
    
    ranges_df=data.frame(
    setups=c("NI_18h","NI_20h","NI_48h"),
    minima_PC1=c(20,20,20),
    maxima_PC1=c(Inf,Inf,Inf),
    minima_PC2=c(Inf,Inf,Inf),
    maxima_PC2=c(Inf,Inf,Inf))
    
    pca_Time_18_20_48_NI=label_PCA(data=pca_Time_18_20_48_NI,ranges=ranges_df,name="Time_18_20_48_NI")
    
    ranges_df=data.frame(
    setups=c("NI_48h","NI_72h"),
    minima_PC1=c(25,25),
    maxima_PC1=c(Inf,Inf),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Time_48_72_NI=label_PCA(data=pca_Time_48_72_NI,ranges=ranges_df,name="Time_48_72_NI")
    
    #pca_Time_2_18_20_TB=label_PCA(data=pca_Time_2_18_20_TB,ranges=ranges_df,name="Time_2_18_20_TB")
    
    ranges_df=data.frame(
    setups=c("TB_18h","TB_20h","TB_48h"),
    minima_PC1=c(20,20,20),
    maxima_PC1=c(Inf,Inf,Inf),
    minima_PC2=c(0,0,-Inf),
    maxima_PC2=c(Inf,Inf,0))
    
    pca_Time_18_20_48_TB=label_PCA(data=pca_Time_18_20_48_TB,ranges=ranges_df,name="Time_18_20_48_TB")
    
    ranges_df=data.frame(
    setups=c("TB_48h","TB_72h"),
    minima_PC1=c(15,15),
    maxima_PC1=c(Inf,Inf),
    minima_PC2=c(Inf,Inf),
    maxima_PC2=c(Inf,Inf))
    
    pca_Time_48_72_TB=label_PCA(data=pca_Time_48_72_TB,ranges=ranges_df,name="Time_48_72_TB")
    
    ranges_df=data.frame(
    setups=c("NI_02h","NI_18h","NI_20h","NI_48h","NI_72h"),
    minima_PC1=c(50,50,50,50,50),
    maxima_PC1=c(Inf,Inf,Inf,Inf,Inf),
    minima_PC2=c(-Inf,-Inf,-Inf,-Inf,-Inf),
    maxima_PC2=c(-30,-30,-30,-30,-30))
    
    pca_Time_all_NI=label_PCA(data=pca_Time_all_NI,ranges=ranges_df,name="Time_all_NI")
    
    ranges_df=data.frame(
    setups=c("TB_02h","TB_18h","TB_20h","TB_48h","TB_72h"),
    minima_PC1=c(50,50,50,50,50),
    maxima_PC1=c(Inf,Inf,Inf,Inf,Inf),
    minima_PC2=c(-Inf,-Inf,-Inf,-Inf,-Inf),
    maxima_PC2=c(-25,-25,-25,-25,-25))
    
    pca_Time_all_TB=label_PCA(data=pca_Time_all_TB,ranges=ranges_df,name="Time_all_TB")
    
    ranges_df=data.frame(
    setups=c("NI_02h","NI_18h","NI_20h","NI_48h","NI_72h","TB_02h","TB_18h","TB_20h","TB_48h","TB_72h"),
    minima_PC1=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
    maxima_PC1=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
    minima_PC2=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
    maxima_PC2=rep(-50,10))
    
    pca_Time_all=label_PCA(data=pca_Time_all,ranges=ranges_df,name="Time_all")
}

##############################################################################
##  10. Upon visual identification of misslabels, I correct their metadata  ##
##############################################################################

{
    correct_misslabeled_condition=function(cols_loc,guy,time){
        
        cols_loc$Original_Treatment=cols_loc$Treatment
        cols_loc$Treatment=as.character(cols_loc$Treatment)
        cols_loc$Original_Setup=cols_loc$Setup
        cols_loc$Setup=as.character(cols_loc$Setup)
        cols_loc$Original_Sample=rownames(cols_loc)

        
        str_NI=paste0("I_",guy,"_NI_",time)
        str_TB=paste0("I_",guy,"_TB_",time)
        
        index_NI=which(cols_loc$Original_Sample==str_NI)
        index_TB=which(cols_loc$Original_Sample==str_TB)
        
        print(paste(index_NI,index_TB))
        cols_loc$Treatment[index_NI]="TB"
        cols_loc$Treatment[index_TB]="NI"
        
        cols_loc$Setup[index_NI]=as.character(cols_loc$Original_Setup)[index_TB]
        cols_loc$Setup[index_TB]=as.character(cols_loc$Original_Setup)[index_NI]
        
        rownames(cols_loc)[index_NI]="a"
        rownames(cols_loc)[index_TB]="b"
        rownames(cols_loc)[index_NI]=cols_loc$Original_Sample[index_TB]
        rownames(cols_loc)[index_TB]=cols_loc$Original_Sample[index_NI]
        
        cols_loc$Sample=rownames(cols_loc)
        
        cols_loc$Individual=factor(cols_loc$Individual)
        cols_loc$Treatment=factor(cols_loc$Treatment)
        cols_loc$TimePoint=factor(cols_loc$TimePoint,levels=c("02h","18h","20h","48h","72h"))
        cols_loc$Setup=factor(cols_loc$Setup)
    
        return(cols_loc[,1:18])
    }
    
    cols=correct_misslabeled_condition(cols_loc=cols,guy="48",time="02h")
    cols=correct_misslabeled_condition(cols_loc=cols,guy="C16",time="18h")
    cols=correct_misslabeled_condition(cols_loc=cols,guy="C052",time="20h")
    cols=correct_misslabeled_condition(cols_loc=cols,guy="29",time="72h")
    cols=correct_misslabeled_condition(cols_loc=cols,guy="C31",time="72h")

    colnames(reads)=rownames(cols)
    
    seguro=whole_groomed
    
    whole_groomed=correct_misslabeled_condition(cols_loc=whole_groomed,guy="48",time="02h")
    whole_groomed=correct_misslabeled_condition(cols_loc=whole_groomed,guy="C16",time="18h")
    whole_groomed=correct_misslabeled_condition(cols_loc=whole_groomed,guy="C052",time="20h")
    whole_groomed=correct_misslabeled_condition(cols_loc=whole_groomed,guy="29",time="72h")
    whole_groomed=correct_misslabeled_condition(cols_loc=whole_groomed,guy="C31",time="72h")
    
    colnames(reads_whole)=rownames(whole_groomed)
}

############################################################################
##    11. Repeat PCAs after correcting misslabels, and label outliers a   ##
##  bit less aggressively (or I would be thrashing too many 20h samples,  ##
##         probably.(Omit all-times-together,non-interesting here)        ##
############################################################################

{
    get_PCA_bis=function(setups,name,type="infection"){
        
        colores_ref=brewer.pal(10,"RdYlBu")[c(6:10,5:1)]
        colores=colores_ref[which(levels(cols$Setup) %in% setups)]
        
        dir.create(file.path(PCA_dir,"clean_bis","plots"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"clean_bis","summaries"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"clean_bis","data"),recursive = T,showWarnings = F)
        
        dir.create(file.path(PCA_dir,"dirty_bis","plots"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"dirty_bis","summaries"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"dirty_bis","data"),recursive = T,showWarnings = F)
        
        set=which(cols$Setup %in% setups)
        cols_chunk=cols[set,]
        reads_chunk=reads[,set]
        cols_chunk$Individual=factor(cols_chunk$Individual)
        cols_chunk$Treatment=factor(cols_chunk$Treatment)
        cols_chunk$TimePoint=factor(cols_chunk$TimePoint)
        
        if(type=="infection"){design=model.matrix(~Individual+Treatment,data=cols_chunk)}
        if(type=="time"){design=model.matrix(~Individual+TimePoint,data=cols_chunk)}
        if(type=="all"){design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols_chunk)}
        
        dge <- DGEList(counts=reads_chunk)
        dge <- calcNormFactors(dge)
        v <- voom(dge,design,plot=FALSE)
        fit <-lmFit(v,design)
        fit <- eBayes(fit)
        model=fit$coefficients
        
        exp_dirty=v$E
        exp_clean=exp_dirty
        
        for(i in 2:length(unique(cols_chunk$Individual)))
        exp_clean=exp_clean-model[,i]%*%t(design[,i])
        
        ### First process dirty
        
        corr_mat <- cor(exp_dirty,method="pearson")
        pca <-prcomp(corr_mat,scale=T)
        sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
        
        datos <- data.frame(pca$x)
        colnames(datos)=paste0("PC",c(1:ncol(datos)))
        rownames(datos)=colnames(exp_dirty)
        datos=merge(datos,cols[colnames(exp_dirty),],by=0)
        
        pl=ggplot(datos) +
        geom_point(data=datos,aes(x=PC1, y=PC2,colour=Setup))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines'))+
            xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
            ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
        
        pdf(paste0(PCA_dir,"/dirty_bis/plots/",name,".pdf"),width=6,height=4)
        print(pl)
        dev.off()
        
        sink(paste0(PCA_dir,"/dirty_bis/summaries/",name,".txt"))
        print(sum_pca)
        sink()
        
        write.table(datos,paste0(PCA_dir,"/dirty_bis/data/",name,".txt"))
        
        datos_dirty=datos
        ### Then process clean
        
        corr_mat <- cor(exp_clean,method="pearson")
        pca <-prcomp(corr_mat,scale=T)
        sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
        
        datos <- data.frame(pca$x)
        colnames(datos)=paste0("PC",c(1:ncol(datos)))
        rownames(datos)=colnames(exp_clean)
        datos=merge(datos,cols[colnames(exp_clean),],by=0)
        
        pl=ggplot(datos) +
        geom_point(data=datos,aes(x=PC1, y=PC2,colour=Setup))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines'))+
            xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
            ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
        
        pdf(paste0(PCA_dir,"/clean_bis/plots/",name,".pdf"),width=6,height=4)
        print(pl)
        dev.off()
        
        sink(paste0(PCA_dir,"/clean_bis/summaries/",name,".txt"))
        print(sum_pca)
        sink()
        
        datos_clean=datos
        write.table(datos_clean,paste0(PCA_dir,"/clean_bis/data/",name,".txt"))
        
        return(list(datos_dirty=datos_dirty,datos_clean=datos_clean,sum_pca=sum_pca))
    }
    
    bis_pca_Infection_2h=get_PCA_bis(setups=c("NI_02h","TB_02h"),name="Infection_2h",type="infection")
    bis_pca_Infection_18h=get_PCA_bis(setups=c("NI_18h","TB_18h"),name="Infection_18h",type="infection")
    bis_pca_Infection_20h=get_PCA_bis(setups=c("NI_20h","TB_20h"),name="Infection_20h",type="infection")
    bis_pca_Infection_48h=get_PCA_bis(setups=c("NI_48h","TB_48h"),name="Infection_48h",type="infection")
    bis_pca_Infection_72h=get_PCA_bis(setups=c("NI_72h","TB_72h"),name="Infection_72h",type="infection")
    
    bis_pca_Time_2_18_20_NI=get_PCA_bis(setups=c("NI_02h","NI_18h","NI_20h"),name="Time_2_18_20_NI",type="time")
    bis_pca_Time_18_20_48_NI=get_PCA_bis(setups=c("NI_18h","NI_20h","NI_48h"),name="Time_18_20_48_NI",type="time")
    bis_pca_Time_48_72_NI=get_PCA_bis(setups=c("NI_48h","NI_72h"),name="Time_48_72_NI",type="time")
    
    bis_pca_Time_2_18_20_TB=get_PCA_bis(setups=c("TB_02h","TB_18h","TB_20h"),name="Time_2_18_20_TB",type="time")
    bis_pca_Time_18_20_48_TB=get_PCA_bis(setups=c("TB_18h","TB_20h","TB_48h"),name="Time_18_20_48_TB",type="time")
    bis_pca_Time_48_72_TB=get_PCA_bis(setups=c("TB_48h","TB_72h"),name="Time_48_72_TB",type="time")
    
    #bis_pca_Time_all_NI=get_PCA_bis(setups=c("NI_02h","NI_18h","NI_20h","NI_48h","NI_72h"),name="Time_all_NI",type="time")
    #bis_pca_Time_all_TB=get_PCA_bis(setups=c("TB_02h","TB_18h","TB_20h","TB_48h","TB_72h"),name="Time_all_TB",type="time")
    #bis_pca_Time_all=get_PCA_bis(setups=levels(cols$Setup),name="Time_all",type="all")
    
    label_PCA_bis=function(data,ranges,name){
        
        dir.create(file.path(PCA_dir,"clean_bis","plots_labelled"),recursive = T,showWarnings = F)
        dir.create(file.path(PCA_dir,"clean_bis","data_labelled"),recursive = T,showWarnings = F)
        
        colores_ref=brewer.pal(10,"RdYlBu")[c(6:10,5:1)]
        colores=colores_ref[which(levels(cols$Setup) %in% ranges$setups)]
        
        datos_dirty=data[[1]]
        datos_clean=data[[2]]
        sum_pca=data[[3]]
        
        datos_clean=put_labels(datos_clean,ranges)
        
        pl=ggplot(datos_clean) +
        geom_point(aes(x=PC1, y=PC2,colour=Setup))+
        geom_text_repel(aes(x=PC1, y=PC2,label=label))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines'))+
            xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
            ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
        
        pdf(paste0(PCA_dir,"/clean_bis/plots_labelled/",name,".pdf"),width=8,height=6)
        print(pl)
        dev.off()
        
        write.table(datos_clean,paste0(PCA_dir,"/clean_bis/data_labelled/",name,".txt"))
        return(list(datos_dirty=datos_dirty,datos_clean=datos_clean,sum_pca=sum_pca))
    }
    {
        ranges_df=data.frame(
        setups=c("NI_02h","TB_02h"),
        minima_PC1=c(0,-Inf),
        maxima_PC1=c(Inf,4),
        minima_PC2=c(Inf,Inf),
        maxima_PC2=c(Inf,Inf))
        
        bis_pca_Infection_2h=label_PCA_bis(data=bis_pca_Infection_2h,ranges=ranges_df,name="Infection_2h")
        
        ranges_df=data.frame(
        setups=c("NI_18h","TB_18h"),
        minima_PC1=c(2,-Inf),
        maxima_PC1=c(Inf,3),
        minima_PC2=c(Inf,Inf),
        maxima_PC2=c(Inf,Inf))
        
        bis_pca_Infection_18h=label_PCA_bis(data=bis_pca_Infection_18h,ranges=ranges_df,name="Infection_18h")
        
        ranges_df=data.frame(
        setups=c("NI_48h","TB_48h"),
        minima_PC1=c(0,-Inf),
        maxima_PC1=c(Inf,5),
        minima_PC2=c(Inf,Inf),
        maxima_PC2=c(Inf,Inf))
        
        bis_pca_Infection_48h=label_PCA_bis(data=bis_pca_Infection_48h,ranges=ranges_df,name="Infection_48h")
        
        ranges_df=data.frame(
        setups=c("NI_72h","TB_72h"),
        minima_PC1=c(0,-Inf),
        maxima_PC1=c(Inf,5),
        minima_PC2=c(Inf,Inf),
        maxima_PC2=c(Inf,Inf))
        
        bis_pca_Infection_72h=label_PCA_bis(data=bis_pca_Infection_72h,ranges=ranges_df,name="Infection_72h")
        
        ranges_df=data.frame(
        setups=c("NI_02h","NI_18h","NI_20h"),
        minima_PC1=c(Inf,Inf,Inf),
        maxima_PC1=c(Inf,Inf,Inf),
        minima_PC2=c(Inf,Inf,-Inf),
        maxima_PC2=c(Inf,Inf,0))
        
        bis_pca_Time_2_18_20_NI=label_PCA_bis(data=bis_pca_Time_2_18_20_NI,ranges=ranges_df,name="Time_2_18_20_NI")
        
        ranges_df=data.frame(
        setups=c("NI_18h","NI_20h","NI_48h"),
        #minima_PC1=c(30,30,30),
        minima_PC1=c(50,50,50),
        maxima_PC1=c(Inf,Inf,Inf),
        minima_PC2=c(Inf,Inf,Inf),
        maxima_PC2=c(Inf,Inf,Inf))
        
        bis_pca_Time_18_20_48_NI=label_PCA_bis(data=bis_pca_Time_18_20_48_NI,ranges=ranges_df,name="Time_18_20_48_NI")
        
        ranges_df=data.frame(
        setups=c("NI_48h","NI_72h"),
        minima_PC1=c(40,40),
        maxima_PC1=c(Inf,Inf),
        minima_PC2=c(Inf,Inf),
        maxima_PC2=c(Inf,Inf))
        
        bis_pca_Time_48_72_NI=label_PCA_bis(data=bis_pca_Time_48_72_NI,ranges=ranges_df,name="Time_48_72_NI")
        
        
        ranges_df=data.frame(
        setups=c("TB_02h","TB_18h","TB_20h"),
        minima_PC1=c(Inf,Inf,Inf),
        maxima_PC1=c(Inf,Inf,Inf),
        minima_PC2=c(Inf,Inf,-Inf),
        maxima_PC2=c(Inf,Inf,-30))
        
        
        bis_pca_Time_2_18_20_TB=label_PCA_bis(data=bis_pca_Time_2_18_20_TB,ranges=ranges_df,name="Time_2_18_20_TB")
        
        ranges_df=data.frame(
        setups=c("TB_18h","TB_20h","TB_48h"),
        minima_PC1=c(25,25,25),
        maxima_PC1=c(Inf,Inf,Inf),
        minima_PC2=c(Inf,Inf,Inf),
        maxima_PC2=c(Inf,Inf,Inf))
        
        
        bis_pca_Time_18_20_48_TB=label_PCA_bis(data=bis_pca_Time_18_20_48_TB,ranges=ranges_df,name="Time_18_20_48_TB")
        
        ranges_df=data.frame(
        setups=c("TB_48h","TB_72h"),
        minima_PC1=c(20,20),
        maxima_PC1=c(Inf,Inf),
        minima_PC2=c(Inf,Inf),
        maxima_PC2=c(Inf,Inf))
        
        bis_pca_Time_48_72_TB=label_PCA_bis(data=bis_pca_Time_48_72_TB,ranges=ranges_df,name="Time_48_72_TB")
    }
    if(FALSE){
        
        ranges_df=data.frame(
        setups=c("NI_02h","NI_18h","NI_20h","NI_48h","NI_72h"),
        minima_PC1=c(50,50,50,50,50),
        maxima_PC1=c(Inf,Inf,Inf,Inf,Inf),
        minima_PC2=c(-Inf,-Inf,-Inf,-Inf,-Inf),
        maxima_PC2=c(-30,-30,-30,-30,-30))
        
        bis_pca_Time_all_NI=label_PCA_bis(data=bis_pca_Time_all_NI,ranges=ranges_df,name="Time_all_NI")
        
        ranges_df=data.frame(
        setups=c("TB_02h","TB_18h","TB_20h","TB_48h","TB_72h"),
        minima_PC1=c(50,50,50,50,30),
        maxima_PC1=c(Inf,Inf,Inf,Inf,Inf),
        minima_PC2=c(-Inf,-Inf,-Inf,-Inf,-Inf),
        maxima_PC2=c(-25,-25,-25,-25,-25))
        
        bis_pca_Time_all_TB=label_PCA_bis(data=bis_pca_Time_all_TB,ranges=ranges_df,name="Time_all_TB")
        
        ranges_df=data.frame(
        setups=c("NI_02h","NI_18h","NI_20h","NI_48h","NI_72h","TB_02h","TB_18h","TB_20h","TB_48h","TB_72h"),
        minima_PC1=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
        maxima_PC1=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf),
        minima_PC2=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),
        maxima_PC2=c(rep(-40,5),rep(-60,5)))
        
        bis_pca_Time_all=label_PCA_bis(data=bis_pca_Time_all,ranges=ranges_df,name="Time_all")
    }
    
    samples_to_remove=unique(c(
    bis_pca_Time_18_20_48_NI[[2]]$label,
    bis_pca_Time_48_72_NI[[2]]$label,
    bis_pca_Time_18_20_48_TB[[2]]$label,
    bis_pca_Time_48_72_TB[[2]]$label))
    
    
    samples_to_remove=samples_to_remove[-which(samples_to_remove=="")]
    
    #####################################################################################
    ##          I will add to samples_to_remove I_C14_NI_72h and I_C13_TB_20h:         ##
    ##               these belong to the only pairs of TB/NI samples that,             ##
    ##   once the PCA outliers are removed, remains clumping together  (I_C14_NI_72h)  ##
    ##        , or still go with|n the other cluster in the tsne (I_C13_TB_20h).       ##
    #####################################################################################
    
    samples_to_remove=c(samples_to_remove,"I_C14_NI_72h","I_C13_TB_18h")
    
    #####################################################################
    ##  samples_to_remove=c("I_C59_NI_48h","I_4_NI_48h","I_4_NI_72h",  ##
    ##     "I_C59_NI_72h","I_C4_TB_18h","I_C5_TB_18h","I_C5_TB_48h",   ##
    ##                  "I_C60_TB_20h","I_C60_TB_48h",                 ##
    ##                    "I_3_TB_48h","I_3_TB_72h",                   ##
    ##    "I_42_TB_48h","I_42_TB_72h","I_C32_TB_48h","I_C32_TB_72h",   ##
    ##                   "I_C5_TB_72h","I_C14_NI_72h")                 ##
    #####################################################################
    
    cols$PCA_criterion=1
    cols$PCA_criterion[which(rownames(cols) %in% samples_to_remove)]=0
    
    whole_groomed$PCA_criterion=1
    whole_groomed$PCA_criterion[which(rownames(whole_groomed) %in% samples_to_remove)]=0

}

##################################################################
##        12. See that removing these last samples do not       ##
##       reduce the time-wise signal (before r´setting the      ##
##        geneset upon final low expression thresholding)       ##
##################################################################

{
    
    ###Check signal before removing them.
    design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols)
    
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=FALSE)
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    betas=fit$coefficients[,c(103:111)]
    ps=fit$p.value[,c(103:111)]
    
    fdrs=ps
    
    for(i in 1:ncol(fdrs)){
        #fdrs[,i]=pmin(rep(1,nrow(ps)),nrow(ps)*(ps[,i]))
        fdrs[,i]=p.adjust(ps[,i])
        #FDR<0.01
        print(paste(colnames(fdrs)[i],length(which(fdrs[,i]<0.01))))
    }
    
    #[1] "TimePoint18h 6291"
    #[1] "TimePoint20h 3570"
    #[1] "TimePoint48h 7541"
    #[1] "TimePoint72h 7736"
    #[1] "TimePoint02h:TreatmentTB 2696"
    #[1] "TimePoint18h:TreatmentTB 6771"
    #[1] "TimePoint20h:TreatmentTB 1646"
    #[1] "TimePoint48h:TreatmentTB 7895"
    #[1] "TimePoint72h:TreatmentTB 8202"
    
    contrast_fitter=function(contrast){
        contrast_last=abs(contrast)+102
        mults=contrast/abs(contrast)
        contrast=contrast_last*mults
        vec=rep(0,length(colnames(design)))
        vec[abs(contrast)]=contrast/abs(contrast)
        fit2 <- contrasts.fit(fit, vec)
        fit2 <- eBayes(fit2)
        res <- topTable(fit2,adjust="BH",n=nrow(reads))
        res$adj.P.Val=p.adjust(res$P.Value)
        print(length(which(res$adj.P.Val<0.01)))
        return(res)
    }
    
    contrast_2_to_18_NI=1
    contrast_2_to_20_NI=2
    contrast_2_to_48_NI=3
    contrast_2_to_72_NI=4
    contrast_18_to_20_NI=c(2,-1)
    contrast_18_to_48_NI=c(3,-1)
    contrast_18_to_72_NI=c(4,-1)
    contrast_20_to_48_NI=c(3,-2)
    contrast_20_to_72_NI=c(4,-2)
    contrast_48_to_72_NI=c(4,-3)
    
    contrast_2_to_18_TB=c(1,6,-5)
    contrast_2_to_20_TB=c(2,7,-5)
    contrast_2_to_48_TB=c(3,8,-5)
    contrast_2_to_72_TB=c(4,9,-5)
    contrast_18_to_20_TB=c(2,7,-1,-6)
    contrast_18_to_48_TB=c(3,8,-1,-6)
    contrast_18_to_72_TB=c(4,9,-1,-6)
    contrast_20_to_48_TB=c(3,8,-2,-7)
    contrast_20_to_72_TB=c(4,9,-2,-7)
    contrast_48_to_72_TB=c(4,9,-3,-8)
    
    contrast_2_to_18_TB_int=c(6,-5)
    contrast_2_to_20_TB_int=c(7,-5)
    contrast_2_to_48_TB_int=c(8,-5)
    contrast_2_to_72_TB_int=c(9,-5)
    contrast_18_to_20_TB_int=c(7,-6)
    contrast_18_to_48_TB_int=c(8,-6)
    contrast_18_to_72_TB_int=c(9,-6)
    contrast_20_to_48_TB_int=c(8,-7)
    contrast_20_to_72_TB_int=c(9,-7)
    contrast_48_to_72_TB_int=c(9,-8)
    
    res_contrast_2_to_18_NI=contrast_fitter(contrast_2_to_18_NI)
    #6291
    res_contrast_2_to_20_NI=contrast_fitter(contrast_2_to_20_NI)
    #3570
    res_contrast_2_to_48_NI=contrast_fitter(contrast_2_to_48_NI)
    #7541
    res_contrast_2_to_72_NI=contrast_fitter(contrast_2_to_72_NI)
    #7736
    res_contrast_18_to_20_NI=contrast_fitter(contrast_18_to_20_NI)
    #53
    res_contrast_18_to_48_NI=contrast_fitter(contrast_18_to_48_NI)
    #2719
    res_contrast_18_to_72_NI=contrast_fitter(contrast_18_to_72_NI)
    #4179
    res_contrast_20_to_48_NI=contrast_fitter(contrast_20_to_48_NI)
    #1022
    res_contrast_20_to_72_NI=contrast_fitter(contrast_20_to_72_NI)
    #1942
    res_contrast_48_to_72_NI=contrast_fitter(contrast_48_to_72_NI)
    #463
    res_contrast_2_to_18_TB=contrast_fitter(contrast_2_to_18_TB)
    #6770
    res_contrast_2_to_20_TB=contrast_fitter(contrast_2_to_20_TB)
    #3943
    res_contrast_2_to_48_TB=contrast_fitter(contrast_2_to_48_TB)
    #7629
    res_contrast_2_to_72_TB=contrast_fitter(contrast_2_to_72_TB)
    #8120
    res_contrast_18_to_20_TB=contrast_fitter(contrast_18_to_20_TB)
    #844
    res_contrast_18_to_48_TB=contrast_fitter(contrast_18_to_48_TB)
    #2151
    res_contrast_18_to_72_TB=contrast_fitter(contrast_18_to_72_TB)
    #4073
    res_contrast_20_to_48_TB=contrast_fitter(contrast_20_to_48_TB)
    #1567
    res_contrast_20_to_72_TB=contrast_fitter(contrast_20_to_72_TB)
    #3518
    res_contrast_48_to_72_TB=contrast_fitter(contrast_48_to_72_TB)
    #781
    
    res_contrast_2_to_18_TB_int=contrast_fitter(contrast_2_to_18_TB_int)
    #[1] 3472 3810
    res_contrast_2_to_20_TB_int=contrast_fitter(contrast_2_to_20_TB_int)
    #[1] 621 590
    res_contrast_2_to_48_TB_int=contrast_fitter(contrast_2_to_48_TB_int)
    #[1] 5062 5134
    res_contrast_2_to_72_TB_int=contrast_fitter(contrast_2_to_72_TB_int)
    #[1] 5809 5982
    res_contrast_18_to_20_TB_int=contrast_fitter(contrast_18_to_20_TB_int)
    #[1] 452 389
    res_contrast_18_to_48_TB_int=contrast_fitter(contrast_18_to_48_TB_int)
    #[1] 825 820
    res_contrast_18_to_72_TB_int=contrast_fitter(contrast_18_to_72_TB_int)
    #[1] 1658 1784
    res_contrast_20_to_48_TB_int=contrast_fitter(contrast_20_to_48_TB_int)
    #[1] 373 319
    res_contrast_20_to_72_TB_int=contrast_fitter(contrast_20_to_72_TB_int)
    #[1] 1184 1202
    res_contrast_48_to_72_TB_int=contrast_fitter(contrast_48_to_72_TB_int)
    #[1] 34 49
    
    summary(cols$Setup)
    #NI_02h NI_18h NI_20h NI_48h NI_72h TB_02h TB_18h TB_20h TB_48h TB_72h
    #98     74     24     95     90     99     74     25     98     91
    
    coef_18_NI=74/(74+24)
    coef_20_NI=24/(74+24)
    
    coef_18_TB=74/(74+25)
    coef_20_TB=25/(74+25)
    
    ### Now remove the PCA outliers and repeat signal estimation
    
    cols_last=cols
    reads_last=reads
    
    cols=cols[-which(rownames(cols) %in% samples_to_remove),]
    reads=reads[,which(colnames(reads) %in% rownames(cols))]
    
    cols$Individual=factor(cols$Individual)
    cols$Treatment=factor(cols$Treatment)
    cols$TimePoint=factor(cols$TimePoint,levels=c("02h","18h","20h","48h","72h"))
    cols$Setup=factor(cols$Setup)
    
    #############################################################################
    ##  Get a signal estimate of contrasts associated with this current setup  ##
    ##                    (i.e. pity samples out, outliers out)                ##
    #############################################################################
    
    design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols)
    
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=FALSE)
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    betas=fit$coefficients[,c(103:111)]
    ps=fit$p.value[,c(103:111)]
    
    fdrs=ps
    
    for(i in 1:ncol(fdrs)){
        #fdrs[,i]=pmin(rep(1,nrow(ps)),nrow(ps)*(ps[,i]))
        fdrs[,i]=p.adjust(ps[,i])
        
        print(paste(colnames(fdrs)[i],length(which(fdrs[,i]<0.01))))
    }
    
    # "TimePoint18h 6443"
    #[1] "TimePoint20h 3684"
    #[1] "TimePoint48h 7594"
    #[1] "TimePoint72h 7863"
    #[1] "TimePoint02h:TreatmentTB 2803"
    #[1] "TimePoint18h:TreatmentTB 7021"
    #[1] "TimePoint20h:TreatmentTB 1635"
    #[1] "TimePoint48h:TreatmentTB 7872"
    #[1] "TimePoint72h:TreatmentTB 8260"

    res_contrast_2_to_18_NI=contrast_fitter(contrast_2_to_18_NI)
    #6291 6443
    res_contrast_2_to_20_NI=contrast_fitter(contrast_2_to_20_NI)
    #3570 3684
    res_contrast_2_to_48_NI=contrast_fitter(contrast_2_to_48_NI)
    #7541 7594
    res_contrast_2_to_72_NI=contrast_fitter(contrast_2_to_72_NI)
    #7736 7863
    res_contrast_18_to_20_NI=contrast_fitter(contrast_18_to_20_NI)
    #53 62
    res_contrast_18_to_48_NI=contrast_fitter(contrast_18_to_48_NI)
    #2719 2775
    res_contrast_18_to_72_NI=contrast_fitter(contrast_18_to_72_NI)
    #4179 4312
    res_contrast_20_to_48_NI=contrast_fitter(contrast_20_to_48_NI)
    #1022 1057
    res_contrast_20_to_72_NI=contrast_fitter(contrast_20_to_72_NI)
    #1942 2124
    res_contrast_48_to_72_NI=contrast_fitter(contrast_48_to_72_NI)
    #463 499
    res_contrast_2_to_18_TB=contrast_fitter(contrast_2_to_18_TB)
    #6770 6930
    res_contrast_2_to_20_TB=contrast_fitter(contrast_2_to_20_TB)
    #3943 3888
    res_contrast_2_to_48_TB=contrast_fitter(contrast_2_to_48_TB)
    #7629 7705
    res_contrast_2_to_72_TB=contrast_fitter(contrast_2_to_72_TB)
    #8120 8159
    res_contrast_18_to_20_TB=contrast_fitter(contrast_18_to_20_TB)
    #844 635
    res_contrast_18_to_48_TB=contrast_fitter(contrast_18_to_48_TB)
    #2151 2100
    res_contrast_18_to_72_TB=contrast_fitter(contrast_18_to_72_TB)
    #4073 4081
    res_contrast_20_to_48_TB=contrast_fitter(contrast_20_to_48_TB)
    #1567 1369
    res_contrast_20_to_72_TB=contrast_fitter(contrast_20_to_72_TB)
    #3518 3329
    res_contrast_48_to_72_TB=contrast_fitter(contrast_48_to_72_TB)
    #781 806
    
    res_contrast_2_to_18_TB_int=contrast_fitter(contrast_2_to_18_TB_int)
    #[1] 3472 3913
    res_contrast_2_to_20_TB_int=contrast_fitter(contrast_2_to_20_TB_int)
    #[1] 621 599
    res_contrast_2_to_48_TB_int=contrast_fitter(contrast_2_to_48_TB_int)
    #[1] 5062 5129
    res_contrast_2_to_72_TB_int=contrast_fitter(contrast_2_to_72_TB_int)
    #[1] 5809 6014
    res_contrast_18_to_20_TB_int=contrast_fitter(contrast_18_to_20_TB_int)
    #[1] 452 418
    res_contrast_18_to_48_TB_int=contrast_fitter(contrast_18_to_48_TB_int)
    #[1] 825 788
    res_contrast_18_to_72_TB_int=contrast_fitter(contrast_18_to_72_TB_int)
    #[1] 1658 1763
    res_contrast_20_to_48_TB_int=contrast_fitter(contrast_20_to_48_TB_int)
    #[1] 373 323
    res_contrast_20_to_72_TB_int=contrast_fitter(contrast_20_to_72_TB_int)
    #[1] 1184 1260
    res_contrast_48_to_72_TB_int=contrast_fitter(contrast_48_to_72_TB_int)
    #[1] 34 63
    
}

###########################################
##  13. Conclussion: Identification of:  ##
###########################################

{
### 1 global threshold for library depth (1E7, before low expressed genes get filtered out.)
### 2. samples with Duplicated=1 should be out.
### 3. samples to remove to recover nesting of individual effects within flowcell
### "10_0_2h"     "10_TB_2h"    "C053_0_72h"  "C13_0_48h.1" "C13_0_72h.1" "C16_0_72h"   "C23_0_18h"   "C23_0_2h"    "C23_TB_18h"  "C23_TB_2h"   "C9_0_18h"    "C9_0_2h"     "C9_TB_18h"   "C9_TB_2h"

### 4. Samples that look as condition misslabels:

## guy="48",time="02h"
## guy="C16",time="18h")
## guy="C052",time="20h")
## guy="29",time="72h")
## guy="C31",time="72h")

### 5. Samples that are outliers in the PCAs: removing these time until 48 h
###  separates in PC1 in adjacent time points wout lowering the signal, as 
###  tested in the part of the coded commented below

## "I_C59_NI_48h" "I_4_NI_48h"   "I_4_NI_72h"   "I_C59_NI_72h" "I_C4_TB_18h"  
## "I_C5_TB_18h"  "I_C5_TB_48h"  "I_C60_TB_20h" "I_C60_TB_48h" "I_3_TB_48h"  
##  "I_3_TB_72h"   "I_42_TB_48h"  "I_42_TB_72h"  "I_C32_TB_48h" "I_C32_TB_72h"
##   "I_C5_TB_72h"

## One more sample that looks like a contamination after doing tsne & PCA,
##  once all other samples have been corrected/removed:

## 6. I_C14_NI_72h I_C13_TB_18h

## On top of that, nothing seems to indicate a technical or biological 
## difference between time=20h and time=18h.

}

###################################################################################
##    14. Build final tables to use: the reads are filtered for low expression   ##
##    using the dataset that verifies Depth, nesting and PCA; and the metadata   ##
##   only has the important, groomed columns. Also, merge timepoints 18 ahd 20.  ##
###################################################################################

{
    ## First, I add one last covariate: the useful depths at the end (using the good sampleset)
    cols$Depth=apply(reads,2,sum)

    ## And add that column to the whole_groomed table
    test=cols[order(rownames(cols)),]
    whole_groomed$index=c(1:nrow(whole_groomed))
    whole_groomed=whole_groomed[order(rownames(whole_groomed)),]
    set=which(rownames(whole_groomed) %in% rownames(test))
    length(which(rownames(test)!=rownames(whole_groomed)[set]))

    whole_groomed$Depth=NA
    whole_groomed$Depth[set]=test$Depth
    test_whole=whole_groomed[which(rownames(whole_groomed) %in% rownames(test)),]
    length(which(test$Depth!=test_whole$Depth))
    whole_groomed=whole_groomed[order(whole_groomed$index),c(1:19,21)]

    ## Leave only the necessary column names to work with in the metadata. (The meaning of each column is further explained in the readme_metadata.txt file)
    cols=cols[,c(8,15,16,17,5,6,14,18, 9,20,10,11)]

    cols$Original_TimePoint=cols$TimePoint
    cols$TimePoint[which(cols$TimePoint=="18h")]="20h"
    cols$TimePoint=factor(cols$TimePoint)
    cols$Time[which(cols$Time==18)]=20
    cols$Setup=paste0(cols$Treatment,"_",cols$TimePoint)
    cols$Sample=paste0(cols$Individual,"_",cols$Setup)
    length(which(colnames(reads)!=rownames(cols)))
    rownames(cols)=cols$Sample
    colnames(reads)=cols$Sample
    cols=cols[,1:12]

    whole_groomed$Original_TimePoint=whole_groomed$TimePoint
    whole_groomed$TimePoint[which(whole_groomed$TimePoint=="18h")]="20h"
    whole_groomed$TimePoint=factor(whole_groomed$TimePoint)
    whole_groomed$Time[which(whole_groomed$Time==18)]=20
    whole_groomed$Setup=paste0(whole_groomed$Treatment,"_",whole_groomed$TimePoint)
    whole_groomed$Sample=paste0(whole_groomed$Individual,"_",whole_groomed$Setup)
    whole_groomed$Sample[which(whole_groomed$Duplicated==1)]=paste0(whole_groomed$Sample[which(whole_groomed$Duplicated==1)],".1")

    whole_groomed=whole_groomed[,c(8,15,16,17,5,6,14,18,9,20,10,11,2,3,4,21,1,7,12,13,19)]

    length(which(colnames(reads_whole)!=rownames(whole_groomed)))
    rownames(whole_groomed)=whole_groomed$Sample
    colnames(reads_whole)=whole_groomed$Sample
}

##################################################################
##   15. Do low expression filtering from the final sampleset   ##
##################################################################

{
    reads=reads_whole[,which(colnames(reads_whole) %in% rownames(cols))]
    
    design=model.matrix(~1,data=cols)
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=TRUE)
    
    exp=data.frame(v$E)
    medians=exp[,1:8]
    
    
    colnames(medians)=unique(cols$Setup)
    
    for(i in 1:ncol(medians))
    {
        case=colnames(medians)[i]
        set=which(cols$Setup %in% case)
        print(paste(i,colnames(medians)[i],length(set)))
        tab=exp[,set]
        medians[,i]=apply(tab,1,median)
    }
    
    medians$max_all=apply(medians,1,max)
    medians$max_used=apply(medians[,c(3,4,5,6,8)],1,max)
    
    th=1
    length(which(medians$max_all>th))
    ## 12920
    th=2
    length(which(medians$max_all>th))
    ## 11494
    length(which(medians$max_all>th & medians$max_used<th))
    ## 298
    length(which(medians$max_all>th & medians$max_used<1))
    ## 20
        
    reads=reads[which(medians$max_all>th),]
}

########################################
##  16. Get dictionary ENSEMBL/HUGOS  ##
########################################

{   #host="jul2018.archive.ensembl.org" Joaquin used this one

    library(biomaRt)
    # host="https://may2024.archive.ensembl.org",
    ensembl <- useMart( 
                       biomart="ENSEMBL_MART_ENSEMBL", 
                       dataset="hsapiens_gene_ensembl")
    # NCBI gene (formerly Entrezgene) accession(s) [e.g. A1BG]
    # uniprot_gn_id
    # entrezgene_accession
    # hgnc_symbol
    # "uniprot_gn_symbol"
    # entrezgene_id
    
    filter_list <- data.frame(listFilters(ensembl))
    
    hugos <- getBM(attributes=c("hgnc_symbol","hgnc_id", "ensembl_gene_id","gene_biotype"),
                   filters = "ensembl_gene_id", values = rownames(reads), 
                   uniqueRows = TRUE, 
                   mart = ensembl)
    
    hugos=hugos[which(hugos$gene_biotype=="protein_coding"),]
    ### Check duplicates by ensemble id
    length(which(duplicated(hugos$ensembl_gene_id)))
    # duplicated_ensembl <- hugos[which(duplicated(hugos$ensembl_gene_id)),]
    # duplicated_entrez_id <- hugos[which(duplicated(hugos$entrezgene_id)),]
    # duplicated_hugo <- hugos[which(duplicated(hugos$hgnc_symbol)),]
    duplicated_hugo_id <- hugos[which(duplicated(hugos$hgnc_symbol)),]
    
    ### SAve duplicated data
    # write.xlsx(duplicated_ensembl,file.path(input_dir,"002_Processed","Duplicated_genes_ensembl.xlsx"))
    # write.table(duplicated_ensembl,file.path(input_dir,"002_Processed","Duplicated_genes_ensembl.txt"))
    # 
    hugos <- hugos[!duplicated(hugos$ensembl_gene_id), ]
    # duplicated_ensembl_gene_id <- hugos[duplicated(hugos$ensembl_gene_id),]
    
    length(which(hugos$hgnc_symbol==""))
    ## 50 non-available gene names, I give them the ensembl ID also in the 
    ## gene name column of the dictionary
    entrez_wo_symbol <- hugos[which(hugos$hgnc_symbol==""),]
    ## duplicated_symbol <- hugos[duplicated(hugos$hgnc_symbol),]
    safe <- hugos
    hugos$hgnc_symbol[which(hugos$hgnc_symbol=="")]=hugos$ensembl_gene_id[which(hugos$hgnc_symbol=="")]
    
    # Check the gene names
    check_genes <- checkGeneSymbols(hugos$hgnc_symbol, unmapped.as.na = F, 
                                    map = NULL, species = "human")
    
    length(which(check_genes$Approved==FALSE & !grepl("ENS",check_genes$x)))
    genes_to_substitute <- check_genes[which(check_genes$Approved==FALSE & !grepl("ENS",check_genes$x)),]
    
    length(which(hugos$hgnc_symbol!=check_genes$x))
    hugos$hgnc_symbol <- check_genes$Suggested.Symbol
    
    ### There are no HUGO duplicated names.
    length(which(duplicated(hugos$hgnc_symbol)))
    # 0

    write.xlsx(hugos,file.path(input_dir,"002_Processed","hugo_sheme_annotation.xlsx"))
    write.table(hugos,file.path(input_dir,"002_Processed","hugo_sheme_annotation.txt"))
    # 
    dictionary=hugos[,1:2]
    dictionary=dictionary[order(dictionary$ensembl_gene_id),]

    reads=reads[which(rownames(reads) %in% dictionary$ensembl_gene_id),]
    reads=reads[order(rownames(reads)),]
    length(which(rownames(reads)!=dictionary$ensembl_gene_id))
    rownames(reads)=dictionary$hgnc_symbol

    dictionary=dictionary[order(dictionary$hgnc_symbol),]
    reads=reads[order(rownames(reads)),]

    cols=cols[order(rownames(cols)),]
    reads=reads[,order(colnames(reads))]

    ### But this is a shitty order for samples
    reads=reads[,order(cols$Setup,cols$Individual)]
    cols=cols[order(cols$Setup,cols$Individual),]

    length(which(colnames(reads)!=rownames(cols)))
    length(which(rownames(reads)!=dictionary$hgnc_symbol))

    ## 0 & 0, so, ok.
}

#################################################################################
##    17. (homogenize genotype samples IDs in the DNA data with Individuals´   ##
##   notation here), and annotate in the RNA metadata those samples having vs  ##
##                            lacking genotype data.                           ##
#################################################################################

{

    dir.create(file.path(input_dir,"002_Processed","genotype_matrix_building_files"),
               showWarnings = F,recursive=T)
    sample_IDs=read.table(file.path(input_dir,"001_Raw_data","genotype_IDs.txt"),header=TRUE)
    sample_IDs_down=read.table(file.path(input_dir,"001_Raw_data","genotype_IDs_tuned_preliminar.txt"),header=TRUE)
    sample_IDs_down$Individual=paste0("I_",sample_IDs_down$Individual)

    metadata=whole_groomed

    ### Aqui ya veo que hay individual IDs en distinto formato en el RNA que en el DNA

    RNA=unique(as.character(metadata$Individual))
    RNA=RNA[order(RNA)]

    DNA=unique(as.character(sample_IDs_down$Individual))
    DNA=DNA[order(DNA)]

    absent_in_DNA=RNA[which(!RNA %in% DNA)]
    absent_in_RNA=DNA[which(!DNA %in% RNA)]

    ### Así que salvo y corrijo a mano.

    length(which(sample_IDs_down$Gencove_Sample_ID !=sample_IDs$Gencove_Sample_ID))
    sample_IDs_down$Gencove_external_IDs=sample_IDs$External_ID
    colnames(sample_IDs_down)[2:4]=paste0("Raw_",colnames(sample_IDs_down)[2:4])
    
    ## Ya que las columnas de ordinal, individual y type van a ser raw, 
    ## dejo vacio el que estaba vacio.
    sample_IDs_down$Raw_Type[which(sample_IDs_down$Raw_Type=="DNA_tuned")]=NA

    sample_IDs_down$Individual=sample_IDs_down$Raw_Individual

    ### La escribo en un csv file para editar a mano los errores
    library(openxlsx)
    write.xlsx(sample_IDs_down,file.path(input_dir,"001_Raw_data","genotype_IDs_preliminar2.xlsx"))

    ### Escribo la columna de individuos tuneada y releo:
    tuned=read.xlsx(file.path(input_dir,"002_Processed","genotype_matrix_building_files","genotype_IDs_tuned.xlsx"))

    RNA=unique(as.character(metadata$Individual))
    RNA=RNA[order(RNA)]

    DNA=unique(as.character(tuned$Individual))
    DNA=DNA[order(DNA)]

    absent_in_DNA=RNA[which(!RNA %in% DNA)]
    absent_in_RNA=DNA[which(!DNA %in% RNA)]

    absent_in_DNA
    # "I_15"     "I_2"      "I_29"     "I_34"     "I_EFS889"
    absent_in_RNA
    ## "I_36_bis"  "I_8_bis"   "I_C29_bis"
    
    ## Likely, one of the two DNA samples named I_C29 correspond to I_29 in 
    ## the RNAseq: in there, we do have a guy named I_29, and another one named I_C29.
    
    whole_groomed$Genotype_criterion=1
    whole_groomed$Genotype_criterion[which(whole_groomed$Individual %in% absent_in_DNA)]=0

}

######################################
##  18. Write ready-to-use tables.  ##
######################################

{
    ## We have two versions of reads and metadata, right now:
    ## 1. reads_whole & cols_whole: they contain all samples and reads
    
    ## 2. reads & cols: they only contain the samples useful for DE analyses, which verify:
    ##
    
    ## reads
    ## cols
    ## dictionary.
    dir.create(file.path(input_dir,"002_Processed","whole"),showWarnings = F)
    dir.create(file.path(input_dir,"002_Processed","ready_for_DE"),showWarnings = F)
    dir.create(file.path(input_dir,"002_Processed","ready_for_EQTL"),showWarnings = F)

    write.table(reads,file.path(input_dir,"002_Processed","ready_for_DE","reads.txt"))
    write.table(cols,file.path(input_dir,"002_Processed","ready_for_DE","metadata.txt"))
    write.table(dictionary,file.path(input_dir,"002_Processed","ready_for_DE","dictionary.txt"))

    ## In the whole tables I do not use HUGOs (didn´t check how many issues are
    ##  in there, but I don´t quite need that)
    write.table(reads_whole,file.path(input_dir,"002_Processed","whole","reads_whole.txt"))
    write.table(reads_whole,file.path(input_dir,"002_Processed","whole","metadata_whole.txt"))
    
    cols_EQTL=whole_groomed[which(whole_groomed$Depth_criterium==1 & whole_groomed$PCA_criterion==1 & whole_groomed$Genotype_criterion==1),]
    reads_EQTL=reads_whole[,which(colnames(reads_whole)%in%rownames(cols_EQTL))]
    reads_EQTL=reads_EQTL[which(rownames(reads_EQTL) %in% dictionary$ensembl_gene_id),]
    
    ## In the files ready for EQTL mapping I will use the same geneset to simplify the situation.
    
    local_dict=dictionary[order(dictionary$ensembl_gene_id),]
    length(which(local_dict$ensembl_gene_id != rownames(reads_EQTL)))
    rownames(reads_EQTL)=local_dict$hgnc_symbol
    
    write.table(reads_EQTL,file.path(input_dir,"002_Processed","ready_for_EQTL","reads.txt"))
    write.table(cols_EQTL,file.path(input_dir,"002_Processed","ready_for_EQTL","metadata.txt"))
    write.table(dictionary,file.path(input_dir,"002_Processed","ready_for_EQTL","dictionary.txt"))

}

########################################################
##  19. Do dim reduction with the final dataset: PCA  ##
########################################################

{
    
get_PCA_tris=function(setups,name,type="infection",shift1=1,shift2=1){
    
    colores_ref=brewer.pal(10,"RdYlBu")[c(7:10,4:1)]
    colores=colores_ref[which(levels(factor(cols$Setup)) %in% setups)]
    
    dir.create(file.path(PCA_dir,"clean_tris","plots"),recursive = T,showWarnings = F)
    dir.create(file.path(PCA_dir,"clean_tris","summaries"),recursive = T,showWarnings = F)
    dir.create(file.path(PCA_dir,"clean_tris","data"),recursive = T,showWarnings = F)
    
    dir.create(file.path(PCA_dir,"dirty_tris","plots"),recursive = T,showWarnings = F)
    dir.create(file.path(PCA_dir,"dirty_tris","summaries"),recursive = T,showWarnings = F)
    dir.create(file.path(PCA_dir,"dirty_tris","data"),recursive = T,showWarnings = F)
    
    set=which(cols$Setup %in% setups)
    cols_chunk=cols[set,]
    reads_chunk=reads[,set]
    cols_chunk$Individual=factor(cols_chunk$Individual)
    cols_chunk$Treatment=factor(cols_chunk$Treatment)
    cols_chunk$TimePoint=factor(cols_chunk$TimePoint)
    
    if(type=="infection"){design=model.matrix(~Individual+Treatment,data=cols_chunk)}
    if(type=="time"){design=model.matrix(~Individual+TimePoint,data=cols_chunk)}
    if(type=="all"){design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols_chunk)}
    
    dge <- DGEList(counts=reads_chunk)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=FALSE)
    fit <-lmFit(v,design)
    fit <- eBayes(fit)
    model=fit$coefficients
    
    exp_dirty=v$E
    exp_clean=exp_dirty
    
    for(i in 2:length(unique(cols_chunk$Individual)))
    exp_clean=exp_clean-model[,i]%*%t(design[,i])
    
    ### First process dirty
    
    corr_mat <- cor(exp_dirty,method="pearson")
    pca <-prcomp(corr_mat,scale=T)
    sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
    
    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    rownames(datos)=colnames(exp_dirty)
    datos=merge(datos,cols[colnames(exp_dirty),],by=0)
    
    
    datos$PC1=shift1*datos$PC1
    datos$PC2=shift2*datos$PC2

    pl=ggplot(datos) +
    geom_point(data=datos,aes(x=PC1, y=PC2,colour=Setup))+
    guides(color = guide_legend(override.aes = list(size=3)))+
    scale_colour_manual(values=colores)+
    #xlim(-6,3)+ylim(-2,2)+
    theme_classic()+
    theme(legend.position="top",
    axis.text.y   = element_text(size=18),
    axis.text.x   = element_text(size=18),
    axis.title.y  = element_text(size=18),
    axis.title.x  = element_text(size=18),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
    legend.title=element_blank(),
    legend.text=element_text(size=16),
    legend.key.size = unit(1, 'lines'))+
        xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
        ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
    
    pdf(paste0(PCA_dir,"/dirty_tris/plots/",name,".pdf"),width=6,height=4)
    print(pl)
    dev.off()
    
    sink(paste0(PCA_dir,"/dirty_tris/summaries/",name,".txt"))
    print(sum_pca)
    sink()
    
    write.table(datos,paste0(PCA_dir,"/dirty_tris/data/",name,".txt"))
    
    datos_dirty=datos
    ### Then process clean
    
    corr_mat <- cor(exp_clean,method="pearson")
    pca <-prcomp(corr_mat,scale=T)
    sum_pca=data.frame(summary(pca)$importance[,c(1:5)])
    
    datos <- data.frame(pca$x)
    colnames(datos)=paste0("PC",c(1:ncol(datos)))
    rownames(datos)=colnames(exp_clean)
    datos=merge(datos,cols[colnames(exp_clean),],by=0)
    
    datos$PC1=shift1*datos$PC1
    datos$PC2=shift2*datos$PC2
    
    pl=ggplot(datos) +
    geom_point(data=datos,aes(x=PC1, y=PC2,colour=Setup))+
    guides(color = guide_legend(override.aes = list(size=3)))+
    scale_colour_manual(values=colores)+
    #xlim(-6,3)+ylim(-2,2)+
    theme_classic()+
    theme(legend.position="top",
    axis.text.y   = element_text(size=18),
    axis.text.x   = element_text(size=18),
    axis.title.y  = element_text(size=18),
    axis.title.x  = element_text(size=18),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
    legend.title=element_blank(),
    legend.text=element_text(size=16),
    legend.key.size = unit(1, 'lines'))+
        xlab(paste0("PC1 (",round(sum_pca$PC1[2],3)*100,"% variance)"))+
        ylab(paste0("PC2 (",round(sum_pca$PC2[2],3)*100,"% variance)"))
    
    pdf(paste0(PCA_dir,"/clean_tris/plots/",name,".pdf"),width=6,height=4)
    print(pl)
    dev.off()
    
    sink(paste0(PCA_dir,"/clean_tris/summaries/",name,".txt"))
    print(sum_pca)
    sink()
    
    write.table(datos,paste0(PCA_dir,"/clean_tris/data/",name,".txt"))
    
    datos_clean=datos
    
    return(list(datos_dirty=datos_dirty,datos_clean=datos_clean,sum_pca=sum_pca))
}

tris_pca_Infection_2h=get_PCA_tris(setups=c("NI_02h","TB_02h"),name="Infection_2h",type="infection")
tris_pca_Infection_20h=get_PCA_tris(setups=c("NI_20h","TB_20h"),name="Infection_20h",type="infection",shift2=-1)
tris_pca_Infection_48h=get_PCA_tris(setups=c("NI_48h","TB_48h"),name="Infection_48h",type="infection")
tris_pca_Infection_72h=get_PCA_tris(setups=c("NI_72h","TB_72h"),name="Infection_72h",type="infection")

tris_pca_Time_2_20_NI=get_PCA_tris(setups=c("NI_02h","NI_20h"),name="Time_2_20_NI",type="time")
tris_pca_Time_20_48_NI=get_PCA_tris(setups=c("NI_20h","NI_48h"),name="Time_20_48_NI",type="time",shift2=-1)
tris_pca_Time_48_72_NI=get_PCA_tris(setups=c("NI_48h","NI_72h"),name="Time_48_72_NI",type="time")

tris_pca_Time_2_20_TB=get_PCA_tris(setups=c("TB_02h","TB_20h"),name="Time_2_20_TB",type="time")
tris_pca_Time_20_48_TB=get_PCA_tris(setups=c("TB_20h","TB_48h"),name="Time_20_48_TB",type="time",shift2=-1)
tris_pca_Time_48_72_TB=get_PCA_tris(setups=c("TB_48h","TB_72h"),name="Time_48_72_TB",type="time")

tris_pca_Time_all_NI=get_PCA_tris(setups=c("NI_02h","NI_20h","NI_48h","NI_72h"),name="Time_all_NI",type="time")
tris_pca_Time_all_TB=get_PCA_tris(setups=c("TB_02h","TB_20h","TB_48h","TB_72h"),name="Time_all_TB",type="time")
tris_pca_Time_all=get_PCA_tris(setups=levels(factor(cols$Setup)),name="Time_all",type="all",shift1=-1,shift2=-1)

}

##########################################################
##  20. Do dim reduction with the final dataset: t-SNE  ##
##########################################################

{
    
    design=model.matrix(~Individual+TimePoint+Treatment:TimePoint,data=cols)
    dge <- DGEList(counts=reads)
    dge <- calcNormFactors(dge)
    v <- voom(dge,design,plot=TRUE)
    fit <-lmFit(v,design)
    fit <- eBayes(fit)

    model=fit$coefficients
    exp_dirty=v$E
    exp_clean=exp_dirty
    
    for(i in 2:length(unique(cols$Individual)))
    exp_clean=exp_clean-model[,i]%*%t(design[,i])
    
    exp_time_out=exp_clean
    for(i in 103:105)
    exp_time_out=exp_time_out-model[,i]%*%t(design[,i])
    
    sds_dirty=apply(exp_dirty,1,sd)
    sds_clean=apply(exp_clean,1,sd)
    sds_time_out=apply(exp_time_out,1,sd)
    
    tab_variances=data.frame(gene=rownames(reads),sd_dirty=sds_dirty,sd_clean=sds_clean,sd_time_out=sds_time_out,index=c(1:nrow(reads)))
    
    do_tsne=function(genes_set=0,input="clean",perp=30,iters=1000,genes=1000,seed){
        
        set.seed(seed)
        tab_variances=tab_variances[order(-tab_variances$sd_dirty),]
        top_n_dirty=rownames(tab_variances)[1:genes]
        tab_variances=tab_variances[order(-tab_variances$sd_clean),]
        top_n_clean=rownames(tab_variances)[1:genes]
        
        tab_variances=tab_variances[order(-tab_variances$sd_time_out),]
        top_n_time_out=rownames(tab_variances)[1:genes]
        
        if(genes_set[1]==0)
        {
            exp_tsne_dirty=exp_dirty[which(rownames(reads) %in% top_n_dirty),]
            exp_tsne_clean=exp_clean[which(rownames(reads) %in% top_n_clean),]
            exp_tsne_time_out=exp_time_out[which(rownames(reads) %in% top_n_time_out),]
        }else{
            exp_tsne_dirty=exp_dirty[which(rownames(reads) %in% genes_set),]
            exp_tsne_clean=exp_clean[which(rownames(reads) %in% genes_set),]
            exp_tsne_time_out=exp_time_out[which(rownames(reads) %in% genes_set),]
        }
        
        if(input=="dirty"){exp=exp_tsne_dirty}
        if(input=="clean"){exp=exp_tsne_clean}
        if(input=="time_out"){exp=exp_tsne_time_out}
        
        tsne_obj <- Rtsne(t(exp), dims = 2, perplexity=perp, verbose=TRUE, max_iter = iters)
        
        tab=data.frame(tsne_obj$Y)
        rownames(tab)=rownames(cols)
        colnames(tab)=c("tsne_1","tsne_2")
        tab=cbind(tab,cols)
        
        colores=brewer.pal(10,"RdYlBu")[c(7:10,4:1)]
        pl=ggplot(tab)+
        geom_point(aes(x=tsne_1, y=tsne_2,colour=Setup))+
        #geom_text_repel(aes(x=PC1, y=PC2,label=label))+
        guides(color = guide_legend(override.aes = list(size=3)))+
        scale_colour_manual(values=colores)+
        #xlim(-6,3)+ylim(-2,2)+
        theme_classic()+
        theme(legend.position="top",
        axis.text.y   = element_text(size=18),
        axis.text.x   = element_text(size=18),
        axis.title.y  = element_text(size=18),
        axis.title.x  = element_text(size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth = 1,linetype="solid"),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.key.size = unit(1, 'lines')
        )
        
        return(list(pl,tab))
    }
    tsne_clean=do_tsne(perp=60,genes=nrow(reads),seed=123)
    tabb=tsne_clean[[2]]
    tabb[(which(tabb$Setup=="NI_72h" & tabb$tsne_1>0)),]
    
    tsne_dir <- file.path(output_dir,"001_Figures","003_TSNE")
    dir.create(tsne_dir,showWarnings = F)
    pdf(file.path(tsne_dir,"tsne.pdf"),width=6,height=6)
    print(tsne_clean[[1]])
    dev.off()
    
    write.table(tsne_clean[[2]],file.path(tsne_dir,"tsne_data.txt"))
}




