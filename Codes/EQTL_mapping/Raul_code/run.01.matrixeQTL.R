# the following Rscript is intended to be submitted to a cluster system
# the accompanying .sbatch file and submission scripts (.sh) files are included 
# and should be modified as necessary (currently, they are written for the slurm scheduler)
# RUN WITH 90G RAM

##########################################################
#### 1. load args & dependencies #########################
args <- commandArgs(TRUE)

## declare condition equal to "NI or "YP" to obtain cis-eQTLs for each condition
condition <- args[1]
## declare pseudobulk type
pseudobulk <- args[2]
## declare which cell type to use
celltype <- args[3]

if(length(args)==0)
  {
    print("WARNING: No arguments supplied.")
  }
print(args)

## number of iteractions for permutations that will be used for FDR correction 
iterations = 10

## load required libraries
library(MatrixEQTL)
library(edgeR)
library(limma)
library(gdsfmt)
library(SNPRelate)
library(qvalue)
library(dplyr)
library(data.table)
library(PCAForQTL)
# devtools::install_github("heatherjzhou/PCAForQTL")

##########################################################
#### 2. load input files & clean/declare output files ####

## load in the expression and counts data
load("./DATA/01.05.quantile_normalized_data.RData")
row.names(meta_data) <- meta_data$data_name
rm(QN_all_data)

genepos = read.table("./DATA/gt_inputs/GRCh38.92_gene_positions.txt",header = TRUE, stringsAsFactors = FALSE)
colnames(genepos)[1:2] <- c("Gene_ID", "chromosome")
snpspos = read.table(paste0("./DATA/gt_inputs/SNP_positions.txt"),header = TRUE, stringsAsFactors = FALSE)
gtypes = read.table(paste0("./DATA/gt_inputs/genotypes.txt"),header = TRUE, stringsAsFactors = FALSE)

## create directory structure to save outputs
if (! dir.exists("./outputs")) {system(paste0("mkdir -p outputs"))}

if (! dir.exists("./outputs/round1_real/")) {system(paste0("mkdir -p outputs/round1_real/"))}
if (! dir.exists("./outputs/round1_perms/")) {system(paste0("mkdir -p outputs/round1_perms/"))}
if (! dir.exists("./outputs/tmpdir/")) {system(paste0("mkdir -p outputs/tmpdir/"))}


############################################################################################################
#### 3. obtain clean input files:                                                                       ####
####   - expression tables output from pseudobulk, where first n principal components are regressed out ####
####   - covariate tables for matrixEQTL, including:                                                    ####
####     2 first genotypes PCs, age, bimodal proportion                                                 ####

# add cell counts to meta_data
counts <- get(paste("counts", pseudobulk, sep="_"))
counts <- as.data.frame(counts[[celltype]])
colnames(counts)[1] <- "ncells"
counts$ncells <- as.numeric(counts$ncells)
counts$data_name <- row.names(counts)
meta_data <- merge(x=meta_data, y=counts, by="data_name")

# scale variables for age and African Admixture
meta_data$age_scale <- scale(meta_data$age)
meta_data$YRI_scale <- scale(meta_data$Afr_admixture)

# get meta data name to match reads name
metadata_whole <- meta_data[order(meta_data$data_name),]
rm(meta_data)

reads_whole = QN_by_condition[[paste(celltype, pseudobulk, sep="_")]]
### subset genes position file to only genes that are present in the RNA-seq data
genepos <- genepos[genepos$Gene_ID %in% rownames(reads_whole),]

## remove duplicate gene positions
genepos <- genepos %>% distinct(Gene_ID, .keep_all = TRUE)

## remove any genes from the RNAseq data which don't have a window
reads_whole <- reads_whole[rownames(reads_whole) %in% genepos$Gene_ID,]

## select only samples for which genotype data is available in the global metadata table and order samples
metadata_whole <- metadata_whole[complete.cases(metadata_whole$Afr_admixture), ]
## remove samples from meta data that are not in the RNA-seq data
metadata_whole <- metadata_whole[which(metadata_whole$data_name %in% colnames(reads_whole)),]
metadata_whole = metadata_whole[order(metadata_whole$data_name),]

## subset the corresponding columns in the reads matrix and order samples and genes
# this removes samples that are in the reads_whole matrix that are not in the metadata_whole matrix
reads_whole = reads_whole[,which(colnames(reads_whole) %in% metadata_whole$data_name)]
# make sure that the samples are ordered
reads_whole = reads_whole[order(rownames(reads_whole)),order(colnames(reads_whole))]

## this is a final check to make sure the order of elements in the metadata_whole and reads_whole
if (sum(rownames(metadata_whole)!=colnames(reads_whole)) != nrow(metadata_whole)) {message("Different numbers of samples in data and metadata (before condition)"); q("no", status=1)} 
  
## the input data matrices are already either voomed-transformed or logCPM, so do not need to redo this

## filter metadata based on condition
metadata = metadata_whole[which(metadata_whole$infection_status==condition),]

## factor variables
metadata$infection_status = factor(metadata$infection_status)
metadata$indiv_ID = factor(metadata$indiv_ID)
metadata$ethnicity = factor(metadata$ethnicity)

## filter voomed expression per condition
reads = reads_whole[,which(colnames(reads_whole) %in% metadata$data_name)]

## check again the coherence of samples order (should be 0)
length(which(colnames(reads)!=metadata$data_name))

if (length(which(colnames(reads)!=metadata$data_name)) != 0) {message("Different numbers of samples in data and metadata"); q("no", status=1)} else { 
  message(paste0(nrow(metadata), " samples for analysis"))}

## shift from sampleIDs to genotyping_IDs (there only will be one sample per genotype in every analysis)
## this removes the condition from the col names of reads and row names of metadata
colnames(reads)=metadata$indiv_ID
rownames(metadata)=metadata$indiv_ID

## recover alphabetical order with the new IDs
reads = reads[,order(colnames(reads))]
metadata = metadata[order(metadata$indiv_ID),]
## check to make sure length = 0
if (length(which(colnames(reads)!=metadata$indiv_ID)) != 0) {message("Different numbers of samples in data and metadata (after renaming)"); q("no", status=1)}

### build matrixEQTL input -- expression tables: regressing out first n PCs from reads
## empirically remove PCs from the phenotype data

# Identify number of PCs to regress out using the elbow method (https://github.com/heatherjzhou/PCAForQTL)

prcompResult<-prcomp(t(reads),center=TRUE,scale.=TRUE)
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
resultRunBE<-PCAForQTL::runBE(reads,B=20,alpha=0.05,mc.cores=1)
K_elbow<-resultRunElbow
K_BE<-resultRunBE$numOfPCsChosen
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE"),values=c(K_elbow,K_BE),
                         titleText=paste("screeplot for", condition, celltype, "cells, using pseudobulk", pseudobulk, "\n remove", K_elbow, "PCs using the elbow method", sep=" "))
ggplot2::ggsave(paste("./outputs/tmpdir/screeplot", celltype, condition, pseudobulk, "jpg", sep="."),width=16,height=11,unit="cm")

# use the number of PCs indicated using the elbow method
pc_set = c(1:K_elbow)

## regress those out
pca_rm <- function(input_data, pc_set) {
    pca = prcomp(t(input_data), na.action = na.omit)
    new = input_data
    new = apply(new, 1, FUN = function(x){return(lm(as.numeric(x) ~ -1 + pca$x[, as.numeric(pc_set)])$resid)})
    new = t(new)
    colnames(new) = colnames(input_data)
    rownames(new) = rownames(input_data)
    return(new)
}
expression = pca_rm(reads, pc_set)

## perform genotype PC analysis and clean genotype data
## set the column names of gtypes to be the samples
samples=colnames(gtypes)
gtypes_pca=data.frame(snp_id=rownames(gtypes),gtypes)

## NOTE: IF YOU RUN INTO AN ERROR: RUN snpgdsClose(genofile) TO RESET FILE OPEN/CLOSE
# snpgdsClose(genofile)

## this creates a SNP genotype dataset from the gtypes_pca matrix
snpgdsCreateGeno(paste("./outputs/tmpdir/GDS_genotypes", celltype, pseudobulk, condition, "gds", sep="."),
                 genmat = as.matrix(gtypes_pca[, samples]),
                 sample.id = unique(samples),
                 snp.id = gtypes_pca$snp_id,
                 snpfirstdim=TRUE)

## this command tells you the total number of samples and SNPs in the .gds file
snpgdsSummary(paste("./outputs/tmpdir/GDS_genotypes", celltype, pseudobulk, condition, "gds", sep="."))

## load .gds file in as genofile
genofile <- snpgdsOpen(paste("./outputs/tmpdir/GDS_genotypes", celltype, pseudobulk, condition, "gds", sep="."))

system(paste("rm ./outputs/tmpdir/GDS_genotypes", celltype, pseudobulk, condition, "gds", sep="."))

## perform a PCA on genotype information
pca <- snpgdsPCA(genofile)
## subset the first 3 PCs
tab <- data.frame(sample.id = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],
                  PC3 = pca$eigenvect[,3],
                  PC4 = pca$eigenvect[,4],
                  PC5 = pca$eigenvect[,5],
                  stringsAsFactors = FALSE)

## create covariates table
## make sure that the pcs_genotypes file is properly labeled and ordered
pcs_genotypes = tab[which(tab$sample.id %in% metadata$indiv_ID),]
pcs_genotypes = pcs_genotypes[order(pcs_genotypes$sample.id),]
if (length(which(metadata$indiv_ID!=pcs_genotypes$sample.id)) != 0) {message("Metadata doesn't match genotype PCs"); q("no", status=1)}


metadata$PC1 = pcs_genotypes$PC1
metadata$PC2 = pcs_genotypes$PC2
metadata$PC3 = pcs_genotypes$PC3
metadata$PC4 = pcs_genotypes$PC4
metadata$PC5 = pcs_genotypes$PC5

tmp <- summary(lm(metadata$Afr_admixture ~ metadata$PC1))$r.squared
message(paste0("PC1 explains ", tmp*100, "% of the variance in African admixture"))

covariates = t(model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + age_scale + ncells, data = metadata))
covariates = covariates[2:nrow(covariates),] # remove row for the intercept

## check to make sure that rownames are correct and match gene_pos
expression = expression[which(rownames(expression) %in% genepos$Gene_ID),]
genepos = genepos[order(genepos$Gene_ID),]

## subset the genotypes file with the individual is present in each condition
genotypes = gtypes[,which(colnames(gtypes) %in% colnames(covariates))]
genotypes = genotypes[which(rownames(genotypes) %in% snpspos$snp),]
genotypes = genotypes[,order(colnames(genotypes))]

if (length(which(rownames(genotypes)!=snpspos$snp)) != 0) {message("Genotype and SNP data don't match"); q("no", status=1)}

#########################################
#### 5. check that input files match ####

## check that all inputs are in the correct order (should all be 0)
length(which(rownames(metadata)!=colnames(covariates)))
length(which(rownames(metadata)!=colnames(expression)))
length(which(rownames(metadata)!=colnames(genotypes)))
length(which(rownames(genotypes)!=snpspos$snp))

## gene-wise check -- this step removes any genes for which there is no expression data (should all be 0)
genepos_trimmed <- genepos[(genepos$Gene_ID %in% rownames(expression)),]
genepos <- genepos_trimmed
length(which(rownames(expression)!=genepos_trimmed$Gene_ID))
length(which(rownames(expression)!=genepos$Gene_ID))

## this output will tell you how many genes for which there is expression data (they should all have the same length)

if (length(genepos_trimmed$Gene_ID) != length(genepos$Gene_ID) | length(genepos$Gene_ID) != length(rownames(expression))) {message("Different numbers of genes"); q("no", status=1)} else { 
  message(paste0(length(genepos_trimmed$Gene_ID), " number of genes"))}


##############################################
#### 6. save matrixEQTL temp input files. ####

expression_file_name = paste("./outputs/tmpdir/expression", celltype, pseudobulk, condition, "txt", sep=".")
covariates_file_name = paste("./outputs/tmpdir/covariates", celltype, pseudobulk, condition, "txt", sep=".")
SNP_file_name = paste("./outputs/tmpdir/genotypes", celltype, pseudobulk, condition, "txt", sep=".")

## write files
write.table(genotypes,SNP_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(expression,expression_file_name, quote=F, sep="\t", row.names=TRUE)
write.table(covariates,covariates_file_name, quote=F, sep="\t", row.names=TRUE)

## in this loop iter = 0 runs the actual analyses, iters 1 to iters 10 run permutations for FDR correction
permuted_pvalues_folder = "outputs/round1_perms/"
for(iteration in 0:iterations){
    
  ###############################################################
  #### 7. permute genotype data (only for iterations > 0) #######
  
  if(iteration > 0){
    
    cols <- colnames(genotypes)
    cols.perm <- sample(cols)
    if(iteration==1){
      random_individuals_df = data.frame(cols.perm)
    }else{
      random_individuals_df = cbind(random_individuals_df,cols.perm)
    }
    genotypes <- genotypes[,cols.perm]
    colnames(genotypes) <- cols
    write.table(genotypes, SNP_file_name, sep="\t", quote = FALSE)
  }

    ##############################################################
    #### 8. prepare & run matrix EQTL ############################
    
    ## load phenotype data
    gene = SlicedData$new();
    gene$fileDelimiter = "\t";     # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values
    gene$fileSkipRows = 1;          # one row of column labels
    gene$fileSkipColumns = 1;       # one column of row labels
    gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
    
    ## load covariates
    cvrt = SlicedData$new();
    cvrt$fileDelimiter = "\t";      # the TAB character
    cvrt$fileOmitCharacters = "NA"; # denote missing values
    cvrt$fileSkipRows = 1;          # one row of column labels
    cvrt$fileSkipColumns = 1;       # one column of row labels
    if(length(covariates_file_name)>0) {
        cvrt$LoadFile(covariates_file_name);
    }
    
    ## load genotype data
    snps = SlicedData$new();
    snps$fileDelimiter = "\t";      # the TAB character
    snps$fileOmitCharacters = "NA";
    snps$fileOmitCharacters = "-9" # denote missing values
    snps$fileSkipRows = 1;          # one row of column labels
    snps$fileSkipColumns = 1;       # one column of row labels
    snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name)
    
    ## set MatrixeQTL options
    useModel = modelLINEAR
    output_file_name_cis = tempfile()
    pvOutputThreshold_cis = 1
    pvOutputThreshold = 0;
    errorCovariance = numeric()
    cisDist = 1e5
    output_file_name = tempfile()
    output_file_name_cis = tempfile()
    
    ## begin MatrixeQTL
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
#### 9. write output files with all sites ##########################
    
    unlink(output_file_name_cis);
    
    if(iteration == 0){
        write.table(me$cis$eqtls, file = paste("./outputs/round1_real/results", celltype, pseudobulk, condition, "txt", sep=".")) }
    else{
        write.table(me$cis$eqtls, file = paste("./outputs/round1_perms/results", iteration, celltype, pseudobulk, condition, "txt", sep="."))
    }
    
}

