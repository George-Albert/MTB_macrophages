##############################################################
#### 1. Load dependencies ####################################
##############################################################
#!/usr/bin/env Rscript
args=(commandArgs(TRUE))
for(i in 1:length(args))
{
    eval(parse(text=args[[i]]))
    print(args[i])
}


#debug=TRUE
### Note: declare condition equal to "CTL, "LPS", or "GARD" to obtain cis-EQTLs for each condition.
#condition="TB"
## Number of random permutations (for empiric fdr corrections).
#setup=94
#time=2
#n_pcs=0
#iterations=2

iters=iterations

file_name=paste0("inputs/processed/EQTL_mapping/samples_sets_tables/tab_",condition,"_",time,"_",setup,".txt")

library(MatrixEQTL)
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)
library(qvalue)
library(stats)

if(iterations>0){
current=getwd()
setwd("codes/common_functions")
source("permFDR.R")
setwd(current)
}
##########################################################
#### 2. Load input files & clean/declare output files ####
##########################################################

genepos = read.table("inputs/processed/EQTL_mapping/gene_positions.txt",header = TRUE, stringsAsFactors = FALSE)

if(debug){
    snpspos = read.table("inputs/processed/EQTL_mapping/further_stuff_to_keep/chunk_SNP_positions.txt",header = TRUE, stringsAsFactors = FALSE)
    gtypes = read.table(paste0("inputs/processed/EQTL_mapping/further_stuff_to_keep/chunk_genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
}else{
    snpspos = read.table(paste0("inputs/processed/EQTL_mapping/SNP_positions.txt"),header = TRUE, stringsAsFactors = FALSE)
    gtypes = read.table(paste0("inputs/processed/EQTL_mapping/genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)
}

system(paste0("mkdir -p Outputs/6_EQTL_mapping_fv/temp_files/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs))
system(paste0("mkdir -p Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/raw_results"))

## Erase and create again the temp file. Not necessary for now.
## system(paste0("rm -rf Outputs/6_EQTL_mapping_fv/temp_files/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs))
## system(paste0("mkdir -p Outputs/6_EQTL_mapping_fv/temp_files/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs))

############################################################################################################
#### 3. Obtain clean input files:                                                                       ####
####   - Voomed Expression tables per condition, where first n principal components are regressed out.  ####
####   - Covariate tables for matrixEQTL, including:                                                    ####
####     2 first genotypes PCs, Flowcell,Sex,fraction of reads assigned and tissue composition          ####
############################################################################################################

#####
##### Clean global tables of expression and metadata
#####


if(setup==94){
    metadata_whole = read.table(paste0("inputs/processed/EQTL_mapping/metadata_w_PCs_94.txt"))
}else{metadata_whole = read.table(paste0("inputs/processed/EQTL_mapping/metadata_w_PCs_71.txt"))
}
#Instead of metadata_whole = read.table(paste0("inputs/processed/EQTL_mapping/global_metadata.txt"))
reads_whole = read.table("inputs/processed/EQTL_mapping/global_expression.txt")
tab_elements=read.table(file_name)

item=load("inputs/processed/EQTL_mapping/samples_sets_tables/individuals_sets.Rdata")
if(setup==94)
{
    inds_set=individuals_set[[1]]
}else{
    inds_set=individuals_set[[2]]
}



## Select only samples for which genotype data is available in the global metadata table, and order samples:
## Esto corría antes, ahora not need metadata_whole=metadata_whole[which(metadata_whole$Individual %in% inds_set),]
rownames(metadata_whole)=paste0(metadata_whole$Individual,"_",metadata_whole$Setup)
metadata_whole=metadata_whole[order(rownames(metadata_whole)),]

## Subset the corresponding columns in the reads matrix, and order samples and genes.
reads_whole=reads_whole[,which(colnames(reads_whole) %in% rownames(metadata_whole))]
reads_whole=reads_whole[,order(colnames(reads_whole))]


## Check that order of elements in the metadata_whole and reads_whole now are congruent, samples-wise
length(which(rownames(metadata_whole)!=colnames(reads_whole)))

## Obtain log(cpm), voomed expression data: the design is irrelevant here
design = model.matrix(~1, data=metadata_whole)
dge <- DGEList(counts=reads_whole)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=FALSE)
voomed_reads_whole = as.data.frame(v$E)

#####
##### Filter expression and metadata per condition.
#####

## Filter metadata per condition
metadata=metadata_whole[which(rownames(metadata_whole) %in% tab_elements$samples),]

## Clean factor variables, and mean center numeric ones.
metadata$Treatment=factor(metadata$Treatment)
metadata$Individual=factor(metadata$Individual)
metadata$Batch=factor(metadata$Batch)
metadata$Original_depth=metadata$Original_depth-mean(metadata$Original_depth)

## Filter voomed expression per condition
voomed_reads=voomed_reads_whole[,which(colnames(voomed_reads_whole) %in% rownames(metadata))]

## Check again the con=herence of samples order
length(which(colnames(voomed_reads)!=rownames(metadata)))


individuals=read.table("inputs/processed/EQTL_mapping/samples_sets_tables/individuals_order.txt")
individuals=individuals[which(individuals$Individual %in% metadata$Individual),]
individuals=individuals[order(individuals$order),2:3]
metadata=merge(metadata,individuals,by="Individual")
rownames(metadata)=paste0(metadata$Individual,"_",metadata$Setup)
metadata=metadata[order(rownames(metadata)),]
voomed_reads=voomed_reads[,order(colnames(voomed_reads))]
length(which(rownames(metadata)!=colnames(voomed_reads)))

voomed_reads=voomed_reads[,order(metadata$order)]
metadata=metadata[order(metadata$order),]
length(which(rownames(metadata)!=colnames(voomed_reads)))


### Shift from sampleIDs to Genotyping_IDs (there only will be one sample per genotype in every analysis).
colnames(voomed_reads)=metadata$Individual
rownames(metadata)=metadata$Individual



### Check that the orders are (and will be, with respect to gtypes, when I'll filter it) ok
length(which(colnames(voomed_reads)!=rownames(metadata)))
samples=colnames(gtypes)[which(colnames(gtypes) %in% colnames(voomed_reads))]
length(which(colnames(voomed_reads)!=samples))


###
### Build matrixEQTL input: expression tables: regressing out first n PCs from voomed_reads:
###



### Regress those out.
pca_rm <- function(input_data, n_pcs) {
    if(n_pcs>0){pc_set=c(1:n_pcs)}else{return(input_data)}
    pca = prcomp(t(input_data), na.action = na.omit)
    new = input_data
    new = apply(new, 1, FUN = function(x){return(lm(as.numeric(x) ~ -1 + pca$x[, as.numeric(pc_set)])$resid)})
    new = t(new)
    colnames(new) = colnames(input_data)
    rownames(new) = rownames(input_data)
    return(new)
}
expression = pca_rm(voomed_reads, n_pcs)

###
### Build matrixEQTL input: covariates tables: transpose metadata & add genotype 1st and 2nd PCs:
###

covariates=t(model.matrix(~PC1+PC2+PC3+PC4+PC5+Batch,data=metadata))
covariates=covariates[2:nrow(covariates),]


##################################################################################################
####  4. Before calling matrixEQTL subset individuals present in the condition to analyze and ####
####     remove genes and SNPs from genotype and expression tables                            ####
####     for which there is no available position (These wouldn't be tested for cis-EQTL,     ####
####     but had useful info to include in PC analyses performed above)                       ####
##################################################################################################

### Subset the genotypes file with the individuals present in each condition
### and the SNPs present in the SNP positions file.
genotypes=gtypes[,which(colnames(gtypes) %in% colnames(covariates))]

###############################################
#### 5. Check input data files congruence. ####
###############################################

## Samples-wise
length(which(rownames(metadata)!=colnames(covariates)))
length(which(rownames(metadata)!=colnames(expression)))
length(which(rownames(metadata)!=colnames(genotypes)))

## SNPs-wise
length(which(rownames(genotypes)!=snpspos$snp))

## Gene-wise
length(which(rownames(expression)!=genepos$Gene_ID))
## All 0, everything is coherent.

##############################################
#### 6. save matrixEQTL temp input files. ####
##############################################


#snps_positions_file_name="Inputs/6_EQTL_mapping/SNP_positions.txt"
#gene_positions_file_name="Inputs/6_EQTL_mapping/gene_positions.txt"
expression_file_name=paste0("Outputs/6_EQTL_mapping_fv/temp_files/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/expression.txt")
covariates_file_name=paste0("Outputs/6_EQTL_mapping_fv/temp_files/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/covariates.txt")
SNP_file_name=paste0("Outputs/6_EQTL_mapping_fv/temp_files/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/genotypes.txt")

write.table(genotypes,SNP_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(expression,expression_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(covariates,covariates_file_name, quote=F, sep="\t", row.names=TRUE)

### In this loop iter=0 runs the actual analyses, iters 1 to iterations (10), run permutation instances for
### FDR correction

permuted_pvalues_folder=paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/raw_results/")
for(iteration in 0:iterations){
    
    print(iteration)
    ##############################################################
    #### 7. Permute genotype data (only for iterations>0 #########
    ##############################################################
    
    if(iteration>0){
        
        cols<-colnames(genotypes)
        cols.perm<-sample(cols)
        if(iteration==1){
        random_individuals_df=data.frame(cols.perm)
        }else{
            random_individuals_df=cbind(random_individuals_df,cols.perm)
        }
        genotypes<-genotypes[,cols.perm]
        colnames(genotypes)<-cols
        write.table(genotypes,SNP_file_name, sep="\t", quote = FALSE)
    }
    
    ##############################################################
    #### 8. Prepare & Run Matrix EQTL ############################
    ##############################################################
    
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";     # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## Load covariates
    
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values;
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
    }
    
    ## Load genotype data
    
    snps = SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA";
    snps$fileOmitCharacters = "-9" # denote missing values;
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name)
    
    useModel = modelLINEAR
    output_file_name_cis = tempfile()
    pvOutputThreshold_cis = 1
    pvOutputThreshold = 0;
    errorCovariance = numeric()
    cisDist = 1e5
    output_file_name = tempfile()
    output_file_name_cis = tempfile()
    
    me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold = pvOutputThreshold,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = FALSE);
    
    ##############################################################
    #### 9. Write temporal output files ##########################
    ##############################################################
    
    unlink(output_file_name_cis);
    
    if(iteration==0){
        write.table(me$cis$eqtls, file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/raw_results/result_original.txt"))}else{
            write.table(me$cis$eqtls, file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/raw_results/result_permuted_",iteration,".txt"))}
        
}


#debug=FALSE
### Note: declare condition equal to "CTL, "LPS", or "GARD" to obtain cis-EQTLs for each condition.
#condition="TB"
## Number of random permutations (for empiric fdr corrections).
#setup=94
#time=2
#n_pcs=0
#iterations=10

iters=iterations

permuted_pvalues_folder=paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/raw_results/")
library(MatrixEQTL)
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)
library(qvalue)
library(stats)

if(iterations>0){
    current=getwd()
    setwd("codes/common_functions")
    source("permFDR.R")
    setwd(current)
}


for(iteration in 0:iterations){
    
    if(iteration==0){file_name="result_original.txt"}else{file_name=paste0("result_permuted_",iteration,".txt")}
    print(paste(iteration,file_name))
    
    event=read.table(paste0(permuted_pvalues_folder,file_name),header=TRUE)
    event.sort<-event[order(event[,2],event[,5]),]
    event.bestQTL<-event.sort[!duplicated(event.sort$gene),]
    event.bestQTL<-event.bestQTL[order(event.bestQTL[,5]),]
    
    if(iteration==0){
        original_best_EQTL=event.bestQTL
        original_all_EQTL=event
    }else{
        if(iteration==1)
        {
            Permutation_Input_t_best = abs(event.bestQTL[3])
            Permutation_Input_p_best = (event.bestQTL[4])
            Permutation_Input_t_all = abs(event[3])
            Permutation_Input_p_all = (event[4])
            
        }else{
            Permutation_Input_t_best=cbind(Permutation_Input_t_best,abs(event.bestQTL[3]))
            Permutation_Input_p_best=cbind(Permutation_Input_p_best,(event.bestQTL[4]))
            Permutation_Input_t_all=cbind(Permutation_Input_t_all,abs(event[3]))
            Permutation_Input_p_all=cbind(Permutation_Input_p_all,(event[4]))
        }
    }
}


#original_all_EQTL=read.table(file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/all_EQTLs.txt"))
#original_best_EQTL=read.table(file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/best_EQTLs.txt"))
#Permutation_Input_t_best=read.table(file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_abs_ts_best.txt"))
#Permutation_Input_p_best=read.table(file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_ps_best.txt"))
#Permutation_Input_t_all=read.table(file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_abs_ts_all.txt"))
#Permutation_Input_p_all=read.table(file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_ps_all.txt"))

write.table(original_all_EQTL,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/all_EQTLs.txt"), quote=FALSE)
write.table(original_best_EQTL,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/best_EQTLs.txt"), quote=FALSE)
write.table(Permutation_Input_t_best,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_abs_ts_best.txt"), quote=FALSE)
write.table(Permutation_Input_p_best,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_ps_best.txt"), quote=FALSE)
write.table(Permutation_Input_t_all,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_abs_ts_all.txt"), quote=FALSE)
write.table(Permutation_Input_p_all,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/permuted_ps_all.txt"), quote=FALSE)

true_ts_all=abs(original_all_EQTL$statistic)
true_ts_best=abs(original_best_EQTL$statistic)

null_ts_all=unlist(Permutation_Input_t_all,use.names=F)
null_ts_best=unlist(Permutation_Input_t_best,use.names=F)

original_all_EQTL$emp_pvalue <- empPvals(stat = true_ts_all, stat0 = null_ts_all)
original_best_EQTL$emp_pvalue <- empPvals(stat = true_ts_best, stat0 = null_ts_best)

original_all_EQTL$emp_BH=p.adjust(original_all_EQTL$emp_pvalue, method = "BH")
original_best_EQTL$emp_BH=p.adjust(original_best_EQTL$emp_pvalue, method = "BH")

original_all_EQTL_sort<-original_all_EQTL[order(original_all_EQTL[,2],original_all_EQTL[,5]),]
original_all_EQTL_sort<-original_all_EQTL_sort[!duplicated(original_all_EQTL_sort$gene),]
original_best_EQTL_sort<-original_all_EQTL_sort[order(original_all_EQTL_sort[,5]),]

## Save before calling qvalue just in case

write.table(original_best_EQTL,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/results_best_SNPs.txt"), quote=FALSE)
write.table(original_all_EQTL,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/results_all_SNPs.txt"), quote=FALSE)
write.table(original_best_EQTL_sort,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/results_best_EQTL_sort.txt"), quote=FALSE)

sink(paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/summary.txt"))
cat(paste0("Statistics for all SNPs:\n"))
#cat(paste0("pi0_=",pi0_all,"\n"))
cat(paste0("Emp_BH Hits @ 1,5 and 10%=",length(which(original_all_EQTL$emp_BH<0.01))," ",length(which(original_all_EQTL$emp_BH<0.05))," ",length(which(original_all_EQTL$emp_BH<0.1)),"\n"))
#cat(paste0("Emp_q Hits @ 1,5 and 10%=",length(which(original_all_EQTL$emp_qvalue<0.01))," ",length(which(original_all_EQTL$emp_qvalue<0.05))," ",length(which(original_all_EQTL$emp_qvalue<0.1)),"\n"))
cat(paste0("Statistics for best SNPs corrected for FDR against best SNPs:\n"))
#cat(paste0("pi0_=",pi0_best,"\n"))
cat(paste0("Emp_BH Hits @ 1,5 and 10%=",length(which(original_best_EQTL$emp_BH<0.01))," ",length(which(original_best_EQTL$emp_BH<0.05))," ",length(which(original_best_EQTL$emp_BH<0.1)),"\n"))
#cat(paste0("Emp_q Hits @ 1,5 and 10%=",length(which(original_best_EQTL$emp_qvalue<0.01))," ",length(which(original_best_EQTL$emp_qvalue<0.05))," ",length(which(original_best_EQTL$emp_qvalue<0.1)),"\n"))
cat(paste0("Statistics for best SNPs subsetted from all SNPs:\n"))
cat(paste0("Emp_BH Hits @ 1,5 and 10% (GENES HAVING EQTL)=",length(which(original_best_EQTL_sort$emp_BH<0.01))," ",length(which(original_best_EQTL_sort$emp_BH<0.05))," ",length(which(original_best_EQTL_sort$emp_BH<0.1)),"\n"))
#cat(paste0("Emp_q Hits @ 1,5 and 10%=",length(which(original_best_EQTL_sort$emp_qvalue<0.01))," ",length(which(original_best_EQTL_sort$emp_qvalue<0.05))," ",length(which(original_best_EQTL_sort$emp_qvalue<0.1)),"\n"))
sink()


qobj_best=qvalue(original_best_EQTL$emp_pvalue)
qobj_all=qvalue(original_all_EQTL$emp_pvalue)

pi0_best=qobj_best$pi0
pi0_all=qobj_all$pi0

original_best_EQTL$emp_qvalue=qobj_best$qvalues
original_all_EQTL$emp_qvalue=qobj_all$qvalues

original_all_EQTL_sort<-original_all_EQTL[order(original_all_EQTL[,2],original_all_EQTL[,5]),]
original_all_EQTL_sort<-original_all_EQTL_sort[!duplicated(original_all_EQTL_sort$gene),]
original_best_EQTL_sort<-original_all_EQTL_sort[order(original_all_EQTL_sort[,5]),]

## Rewrite, if everything worked well, qith the qvalues in there.

write.table(original_best_EQTL,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/results_best_SNPs.txt"), quote=FALSE)
write.table(original_all_EQTL,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/results_all_SNPs.txt"), quote=FALSE)
write.table(original_best_EQTL_sort,file = paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/results_best_EQTL_sort.txt"), quote=FALSE)

sink(paste0("Outputs/6_EQTL_mapping_fv/",condition,"/time_",time,"/setup_",setup,"/iters_",iters,"/n_pcs_",n_pcs,"/summary.txt"))
cat(paste0("Statistics for all SNPs:\n"))
cat(paste0("pi0_=",pi0_all,"\n"))
cat(paste0("Emp_BH Hits @ 1,5 and 10%=",length(which(original_all_EQTL$emp_BH<0.01))," ",length(which(original_all_EQTL$emp_BH<0.05))," ",length(which(original_all_EQTL$emp_BH<0.1)),"\n"))
cat(paste0("Emp_q Hits @ 1,5 and 10%=",length(which(original_all_EQTL$emp_qvalue<0.01))," ",length(which(original_all_EQTL$emp_qvalue<0.05))," ",length(which(original_all_EQTL$emp_qvalue<0.1)),"\n"))
cat(paste0("Statistics for best SNPs corrected for FDR against best SNPs:\n"))
cat(paste0("pi0_=",pi0_best,"\n"))
cat(paste0("Emp_BH Hits @ 1,5 and 10%=",length(which(original_best_EQTL$emp_BH<0.01))," ",length(which(original_best_EQTL$emp_BH<0.05))," ",length(which(original_best_EQTL$emp_BH<0.1)),"\n"))
cat(paste0("Emp_q Hits @ 1,5 and 10%=",length(which(original_best_EQTL$emp_qvalue<0.01))," ",length(which(original_best_EQTL$emp_qvalue<0.05))," ",length(which(original_best_EQTL$emp_qvalue<0.1)),"\n"))
cat(paste0("Statistics for best SNPs subsetted from all SNPs:\n"))
cat(paste0("Emp_BH Hits @ 1,5 and 10% (GENES HAVING EQTL)=",length(which(original_best_EQTL_sort$emp_BH<0.01))," ",length(which(original_best_EQTL_sort$emp_BH<0.05))," ",length(which(original_best_EQTL_sort$emp_BH<0.1)),"\n"))
cat(paste0("Emp_q Hits @ 1,5 and 10%=",length(which(original_best_EQTL_sort$emp_qvalue<0.01))," ",length(which(original_best_EQTL_sort$emp_qvalue<0.05))," ",length(which(original_best_EQTL_sort$emp_qvalue<0.1)),"\n"))
sink()

print("Done!")



