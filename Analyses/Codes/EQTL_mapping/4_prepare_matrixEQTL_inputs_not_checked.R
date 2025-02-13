
###
### Get set of usable genes, their positions and the read file that contains only them
###
library(openxlsx)

reads=read.table("inputs/processed/ready/reads_ready.txt")
cols=read.table("inputs/processed/ready/metadata_ready.txt")
dictionary=read.table("inputs/processed/ready/dictionary_ready.txt")
positions_g=read.table("inputs/processed/EQTL_mapping/further_stuff_to_keep/matrixEQTL/gene_positions.txt")
snp_pos=read.table("inputs/processed/EQTL_mapping/Genotypes/genotype_matrix_all_MAF_10_missing_50.012.pos")
snp_pos_g=read.table("inputs/processed/EQTL_mapping/further_stuff_to_keep/matrixEQTL/SNP_positions.txt")
individuals=read.table("inputs/processed/EQTL_mapping/Genotypes/genotype_matrix_all_MAF_10_missing_50.012.indv")
individuals_meta=read.xlsx("inputs/processed/genotype_matrix_building_files/genotype_IDs_tuned.xlsx")
## Check that the dictionary is coherent with the reads file:
length(which(rownames(reads)==dictionary$hgnc_symbol))
# 10833
dim(reads)
# 10833   751
## It is ok.

library("biomaRt")
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

pos = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = as.character(dictionary$ensembl_gene_id) , mart = human, attributesL = c("chromosome_name","start_position","end_position"), martL = human, uniqueRows=T)

dim(pos)
# 10760     4
length(unique(pos$Gene.stable.ID))
# 10760

## Select only autosomes
pos=pos[which(pos$Chromosome.scaffold.name %in% c(as.character(c(1:22)))),]
nrow(pos)
# 10357

dictionary=dictionary[which(dictionary$ensembl_gene_id %in% pos$Gene.stable.ID),]
dictionary=dictionary[order(dictionary$ensembl_gene_id),]
pos=pos[order(pos$Gene.stable.ID),]
length(which(pos$Gene.stable.ID==dictionary$ensembl_gene_id))
# 10357
dim(pos)
# 10357     4

pos$hgnc_symbol=dictionary$hgnc_symbol
pos=pos[,c(1,5,2,3,4)]
colnames(pos)=c("ensembl_gene_id","Gene_ID","chromosome","S1","S2")

reads_sub=reads[which(rownames(reads) %in% pos$Gene_ID),]
reads_sub=reads_sub[order(rownames(reads_sub)),]

pos$chromosome=paste0("chr",pos$chromosome)
pos=pos[order(pos$Gene_ID),]
length(which(rownames(reads_sub)!=pos$Gene_ID))

pos_ready=pos[,2:5]

dim(reads_sub)
# 10357   751
dim(pos_ready)
# 10357     4
length(which(rownames(reads_sub)!=pos_ready$Gene_ID))
# 0
length(which(rownames(reads_sub)==pos_ready$Gene_ID))
# 10357

system("mkdir -p inputs/processed/EQTL_mapping/further_stuff_to_keep")
write.table(pos,"inputs/processed/EQTL_mapping/further_stuff_to_keep/gene_positions_w_ensembl_IDs.txt")
write.table(pos_ready,"inputs/processed/EQTL_mapping/gene_positions.txt")
write.table(reads_sub,"inputs/processed/EQTL_mapping/global_expression.txt")
##This is a copy untouched of what we had before, but still
write.table(cols,"inputs/processed/EQTL_mapping/global_metadata.txt")

###
### Write SNP positions file: I will not match exactly the format we had in the HGs (no ref allele info on ID).
###

snp_pos$snp=paste0(snp_pos[,1],"_",snp_pos[,2])
snp_pos=snp_pos[,c(3,1,2)]
colnames(snp_pos)=c("snp","chr","pos")
snp_pos$chr=paste0("chr",snp_pos$chr)

write.table(snp_pos,"inputs/processed/EQTL_mapping/SNP_positions.txt")

###
###
###

individuals$order=c(1:nrow(individuals))
individuals=individuals[order(individuals$V1),]
individuals_meta=individuals_meta[order(individuals_meta$Gencove_external_IDs),]
length(which(individuals_meta$Gencove_external_IDs != individuals$V1))
individuals$Individual=individuals_meta$Individual
individuals=individuals[order(individuals$order),]

gtypes <- fread("inputs/processed/EQTL_mapping/Genotypes/genotype_matrix_all_MAF_10_missing_50 copy_012.txt",header=FALSE)
gtypes=data.frame(t(gtypes))

gtypes=gtypes[c(2:nrow(gtypes)),]
colnames(gtypes)=as.character(individuals$Individual)
# AHora intento esto, que falla
#rownames(gtypes)=snp_pos$snp
#Error in `row.names<-.data.frame`(`*tmp*`, value = value) :
#duplicate 'row.names' are not allowed
#In addition: Warning message:
#non-unique value when setting 'row.names': ‘8_37679934’
## Hay algún SNP reportado en la misma posición.

length(unique(snp_pos$snp))
# 6159737: sólo uno repe!!
which(duplicated(snp_pos$snp))
# 5672103
snp_pos$snp[5672103]
# "8_37679934"
snp_pos[which(snp_pos$snp=="8_37679934"),]
# snp  chr      pos
# 5672102 8_37679934 chr8 37679934
# 5672103 8_37679934 chr8 37679934
t(gtypes[c(5672102:5672103),])
## Viendo esto se ve que los dos SNPs son distintos: los quito y ya.

snp_pos_clean=snp_pos[-c(5672102:5672103),]

gtypes=gtypes[-c(5672102:5672103),]
rownames(gtypes)=snp_pos_clean$snp

for(i in 1:100){print(length(which(!gtypes[,i] %in% c(0,1,2))))}
## There are no missing entries.

write.table(gtypes,"inputs/processed/EQTL_mapping/genotypes.txt")
write.table(snp_pos_clean,"inputs/processed/EQTL_mapping/SNP_positions.txt")

set_to_choose=sample(1:nrow(gtypes),1000)

gtypes_chunk=gtypes[set_to_choose,]
snp_pos_chunk=snp_pos_clean[set_to_choose,]

write.table(gtypes_chunk,"inputs/processed/EQTL_mapping/further_stuff_to_keep/chunk_genotypes.txt")
write.table(snp_pos_chunk,"inputs/processed/EQTL_mapping/further_stuff_to_keep/chunk_SNP_positions.txt")

### Test bis samples

length(which(gtypes[,which(colnames(gtypes)=="I_8")]!=gtypes[,which(colnames(gtypes)=="I_8_bis")]))
# 2842507
length(which(gtypes[,which(colnames(gtypes)=="I_C29")]!=gtypes[,which(colnames(gtypes)=="I_C29_bis")]))
# 2841145
length(which(gtypes[,which(colnames(gtypes)=="I_36")]!=gtypes[,which(colnames(gtypes)=="I_36_bis")]))
# 2846519

###
### Generate sets of individuals to subset for each analysis
###

individuals_no_dup=individuals_meta$Individual[which(is.na(individuals_meta$Comment) | individuals_meta$Comment=="ID_corrected")]
individuals_clean=individuals_meta$Individual[which(is.na(individuals_meta$Comment))]



## 94 and 71 guys, respectively.

###
### Set de muestras e individuos.
###



get_indexes_exp=function(condition,times,inds,name){
    
    samples=rownames(cols)[which(cols$Treatment==condition & cols$TimePoint %in% times & cols$Individual %in% inds)]
    guys=cols$Individual[which(cols$Treatment==condition & cols$TimePoint %in% times & cols$Individual %in% inds)]
    
    print(paste0("Samples: ",length(samples),"(",length(unique(samples)),") guys:",length(guys),"(",length(unique(guys)),")"))
    
    out=data.frame(samples=samples,guys=guys)
    system("mkdir -p inputs/processed/EQTL_mapping/samples_sets_tables")
    write.table(out,paste0("inputs/processed/EQTL_mapping/samples_sets_tables/",name,".txt"))
    return(out)
}

tab_NI_2_94=get_indexes_exp("NI","02h",individuals_no_dup,"tab_NI_2_94")
# "Samples: 90(90) guys:90(90)"
tab_TB_2_94=get_indexes_exp("TB","02h",individuals_no_dup,"tab_TB_2_94")
# "Samples: 91(91) guys:91(91)"
tab_NI_20_94=get_indexes_exp("NI",c("18h","20h"),individuals_no_dup,"tab_NI_20_94")
# "Samples: 90(90) guys:90(90)"
tab_TB_20_94=get_indexes_exp("TB",c("18h","20h"),individuals_no_dup,"tab_TB_20_94")
# "Samples: 88(88) guys:88(88)"
tab_NI_48_94=get_indexes_exp("NI","48h",individuals_no_dup,"tab_NI_48_94")
# "Samples: 85(85) guys:85(85)"
tab_TB_48_94=get_indexes_exp("TB","48h",individuals_no_dup,"tab_TB_48_94")
# "Samples: 85(85) guys:85(85)"
tab_NI_72_94=get_indexes_exp("NI","72h",individuals_no_dup,"tab_NI_72_94")
# "Samples: 82(82) guys:82(82)"
tab_TB_72_94=get_indexes_exp("TB","72h",individuals_no_dup,"tab_TB_72_94")
# "Samples: 81(81) guys:81(81)"

tab_NI_2_71=get_indexes_exp("NI","02h",individuals_clean,"tab_NI_2_71")
tab_TB_2_71=get_indexes_exp("TB","02h",individuals_clean,"tab_TB_2_71")
tab_NI_20_71=get_indexes_exp("NI",c("18h","20h"),individuals_clean,"tab_NI_20_71")
tab_TB_20_71=get_indexes_exp("TB",c("18h","20h"),individuals_clean,"tab_TB_20_71")
tab_NI_48_71=get_indexes_exp("NI","48h",individuals_clean,"tab_NI_48_71")
tab_TB_48_71=get_indexes_exp("TB","48h",individuals_clean,"tab_TB_48_71")
tab_NI_72_71=get_indexes_exp("NI","72h",individuals_clean,"tab_NI_72_71")
tab_TB_72_71=get_indexes_exp("TB","72h",individuals_clean,"tab_TB_72_71")
# "Samples: 67(67) guys:67(67)"
# "Samples: 68(68) guys:68(68)"
# "Samples: 68(68) guys:68(68)"
# "Samples: 66(66) guys:66(66)"
# "Samples: 63(63) guys:63(63)"
# "Samples: 63(63) guys:63(63)"
# "Samples: 62(62) guys:62(62)"
# "Samples: 59(59) guys:59(59)"


get_indexes_resp=function(times,inds,name){
    
    guys_NI=cols$Individual[which(cols$Treatment=="NI" & cols$TimePoint %in% times & cols$Individual %in% inds)]
    guys_TB=cols$Individual[which(cols$Treatment=="TB" & cols$TimePoint %in% times & cols$Individual %in% inds)]
    
    guys=guys_NI[which(guys_NI %in% guys_TB)]
    
    print(paste0("Guys:",length(guys),"(",length(unique(guys)),")"))
    
    out=data.frame(guys=guys)
    system("mkdir -p inputs/processed/EQTL_mapping/samples_sets_tables")
    write.table(out,paste0("inputs/processed/EQTL_mapping/samples_sets_tables/",name,"_resp.txt"))
    return(out)
}

times="02h"
inds=individuals_no_dup
name="tab_2_94"

tab_2_94=get_indexes_resp("02h",individuals_no_dup,"tab_2_94")
tab_20_94=get_indexes_resp(c("18h","20h"),individuals_no_dup,"tab_20_94")
tab_48_94=get_indexes_resp("48h",individuals_no_dup,"tab_48_94")
tab_72_94=get_indexes_resp("72h",individuals_no_dup,"tab_72_94")

# "Guys:90(90)"
# "Guys:87(87)"
# "Guys:81(81)"
# "Guys:74(74)"

tab_2_71=get_indexes_resp("02h",individuals_clean,"tab_2_71")
tab_20_71=get_indexes_resp(c("18h","20h"),individuals_clean,"tab_20_71")
tab_48_71=get_indexes_resp("48h",individuals_clean,"tab_48_71")
tab_72_71=get_indexes_resp("72h",individuals_clean,"tab_72_71")

# "Guys:67(67)"
# "Guys:66(66)"
# "Guys:60(60)"
# "Guys:55(55)"

individuals_set=list(no_dup=individuals_no_dup,clean=individuals_clean)
save(individuals_set,file="inputs/processed/EQTL_mapping/samples_sets_tables/individuals_sets.Rdata")
write.table(individuals,"inputs/processed/EQTL_mapping/samples_sets_tables/individuals_order")

###
### Write scripts
###
scripter=function(condition,times,setup,iters,debug=FALSE,n_pcs=0){
    
    tiempos="02:00:00"
    memo="10G"
    name=paste0(condition,"_",times[length(times)],"_",setup,"_",iters,"_",n_pcs,"_deb",as.numeric(debug))
    system("mkdir -p codes/scripts_EQTLs_fdr")
    script_name=paste0("codes/scripts_EQTLs_fdr/",name,".sh")
    sink(script_name)
    cat("#!/bin/bash\n")
    cat("#SBATCH --account=def-barreiro\n")
    cat(paste0("#SBATCH --time=",tiempos,"\n"))
    cat(paste0("#SBATCH --job-name=",name,"\n"))
    cat("#SBATCH --output=logs_fdr/%x-%j.out\n")
    cat(paste0("#SBATCH --mem=",memo,"\n"))
    cat("module load nixpkgs/16.09  gcc/7.3.0 r/3.4.4\n")
    
    cat(paste0("R CMD BATCH --no-save --no-restore '--args debug=",debug," condition=\"",condition,"\" setup=",setup," time=",times[length(times)]," n_pcs=",n_pcs," iterations=",iters,"' codes/6_eQTL_FDR_corrections.R logs_fdr/log_",name,".out\n"))
    sink()
}




scripter(condition="NI",times=2,setup=94,iters=10)
scripter(condition="TB",times=2,setup=94,iters=10)
scripter(condition="NI",times=c(18,20),setup=94,iters=10)
scripter(condition="TB",times=c(18,20),setup=94,iters=10)
scripter(condition="NI",times=48,setup=94,iters=10)
scripter(condition="TB",times=48,setup=94,iters=10)
scripter(condition="NI",times=72,setup=94,iters=10)
scripter(condition="TB",times=72,setup=94,iters=10)



