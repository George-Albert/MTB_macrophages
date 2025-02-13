###
### 1. Download gencove files´urls, & Prepare list of commands to run in order to download & tabix them: one row per guy.
### Generates: /inputs/processed/genotype_matrix_building_files/genotypes_matrix_building_commands.txt  (Reproducible)

{
    library(openxlsx)
    
    tuned=read.xlsx("inputs/processed/genotype_matrix_building_files/genotype_IDs_tuned.xlsx")

    ### Write txt column table with gencove IDs to retrieve the presigned URLs to the cvf files:

    sink("inputs/processed/genotype_matrix_building_files/Gencove_external_IDs.txt")
    for(id in tuned$Gencove_external_IDs)
    cat(id,"\n")
    sink()

    system("gencove project raw_data --api-key f28cad3e-71e1-407c-b89a-f8644752b440 1067 inputs/processed/genotype_matrix_building_files/Gencove_external_IDs.txt inputs/processed/genotype_matrix_building_files/urls_ludo_Jun4.csv")

    urls=read.csv("inputs/processed/genotype_matrix_building_files/urls_ludo_Jun4.csv")

    urls_useful=urls[,c(1,2,5)]
    urls_useful$vcf_url_s3=as.character(urls_useful$vcf_url_s3)
    tuned=tuned[order(tuned$Gencove_Sample_ID),]
    urls_useful=urls_useful[order(urls_useful$sample_id),]
    length(which(tuned$Gencove_Sample_ID!=urls_useful$sample_id))
    urls_useful$new_name=tuned$Individual

    urls_useful$download_gz_command=paste0("wget -O inputs/gzs/",urls_useful$new_name,".vcf.gz ", "\"",urls_useful$vcf_url_s3,"\"")
    urls_useful$tabix_command=paste0("tabix -p vcf inputs/gzs/",urls_useful$new_name,".vcf.gz")

    write.table(urls_useful,"inputs/processed/genotype_matrix_building_files/genotypes_matrix_building_commands.txt")

}

###
### 2. Download & tabix vcf files.
### Generates files in inputs/gzs/ named individual.vcf.gz  (100 of these, Beware: I didn´t keep these)

{
    ## Make sure you load the following modules before running the loop.

    ## module load nixpkgs/16.09  gcc/7.3.0 r/3.4.4
    ## module load nixpkgs/16.09  intel/2018.3 vcftools
    ## module load nixpkgs/16.09  gcc/7.3.0 r/3.4.4
    ## module load bamtools/2.4.1
    ## module load htslib/1.9

    commands=read.table("inputs/processed/genotype_matrix_building_files/genotypes_matrix_building_commands.txt")

    for(i in 1:nrow(commands))
    {
        print(i)
        system(as.character(commands$download_gz_command[i]))
        system(as.character(commands$tabix_command[i]))
    }
}
###
### 3. Write scripts with code to remove indels, per chromosome, and then merge chromosomes.
### The scripts generated go to codes/process_genotype_scripts
### Each script process one chromosome: it reads the 100 inputs/gzs/individual.vcf.gz, removing indels, and then:
### First: save each chunk: (individual-chromosome), in inputs/processed_chr_chunks
### Then: merge all guys for each chromosome in merged_chr, to generate one file per chromosome for all samples.
### This may feel weird but it was the fastest way to compute this I found.


{
    ## First we write a script per chromosome:
    system("mkdir -p codes/process_genotype_scripts")
    system("mkdir -p inputs/files_to_merge_lists")
    system("mkdir -p inputs/merged_chrs")
    system("mkdir -p inputs/processed_chr_chunks")

    chr_groups=paste0("/ch_",c(1:22))

    for(i in 1:22)
    {
        print(paste("Chromosome ",i))
        
        
        sink(paste0("inputs/files_to_merge_lists/chromosome_",i,"_files.txt"))
        
        for(j in 1:100)
        {
            cat(paste0("inputs/processed_chr_chunks",chr_groups[i],"/",as.character(commands$new_name)[j],".SNPs_only_ch_",i,".recode.vcf.gz "))
        }
        sink()
        
        
        sink(paste0("codes/process_genotype_scripts/chromosome_",i,".sh"))
        cat("#!/bin/bash\n")
        cat("#SBATCH --account=def-barreiro\n")
        cat("#SBATCH --time=10:00:00\n")
        cat(paste0("#SBATCH --job-name=ch_",i,"\n"))
        cat("#SBATCH --output=logs/%x-%j.out\n")
        cat("#SBATCH --mem=50G\n")
        cat("module load nixpkgs/16.09  intel/2018.3 vcftools htslib/1.9\n")

        cat(paste0("mkdir -p inputs/processed_chr_chunks",chr_groups[i],"\n"))
        cat("for guy in ",as.character(commands$new_name),"\n")
        cat("do\n")
        #cat("echo \"Welcome $guy times\"\n")
        cat(paste0("vcftools  --gzvcf  inputs/gzs/$guy.vcf.gz  --remove-indels --chr ",i,"  --recode --recode-INFO-all --out  inputs/processed_chr_chunks",chr_groups[i],"/$guy.SNPs_only_ch_",i,"\n"))
        cat(paste0("bgzip < inputs/processed_chr_chunks",chr_groups[i],"/$guy.SNPs_only_ch_",i,".recode.vcf > inputs/processed_chr_chunks",chr_groups[i],"/$guy.SNPs_only_ch_",i,".recode.vcf.gz\n"))
        cat(paste0("rm inputs/processed_chr_chunks",chr_groups[i],"/$guy.SNPs_only_ch_",i,".recode.vcf\n"))
        cat(paste0("tabix -p vcf inputs/processed_chr_chunks",chr_groups[i],"/$guy.SNPs_only_ch_",i,".recode.vcf.gz\n"))
        cat("done\n")
        cat(paste0("vcf-merge $(cat inputs/files_to_merge_lists/chromosome_",i,"_files.txt) | bgzip -c > inputs/merged_chrs/merged_chr",i,".vcf.gz\n"))
        #cat(paste0("rm -rf inputs/processed_chr_chunks",chr_groups[i],"\n"))
        sink()
    }
}

system("mkdir -p inputs/merged_chrs")

## Now run them; to do so you may need to close R to ensure these modules are loaded,.

    ## module load intel/2018.3 r/3.4.4 vcftools
    ## module load nixpkgs/16.09  gcc/7.3.0 r/3.4.4
    ## module load nixpkgs/16.09  intel/2018.3 vcftools
    ## module load nixpkgs/16.09  gcc/7.3.0 r/3.4.4
    ## module load bamtools/2.4.1
    ## module load htslib/1.9


## Each of these scripts below submit a job to execute a group of chromosomes. For example, group 1 has chrs 1 and 2:

#!/bin/bash
#SBATCH --account=def-someuser
#SBATCH --time=00:01:00
#SBATCH --job-name=group_1
#SBATCH --output=%x-%j.out
#sh codes/process_genotype_scripts/chromosome_1.sh
#sh codes/process_genotype_scripts/chromosome_2.sh


system("sh codes/process_genotype_scripts/master_codes/group_1.sh")
system("sh codes/process_genotype_scripts/master_codes/group_2.sh")
system("sh codes/process_genotype_scripts/master_codes/group_3.sh")
system("sh codes/process_genotype_scripts/master_codes/group_4.sh")
system("sh codes/process_genotype_scripts/master_codes/group_5.sh")

## And now you concatenate all chromosomes into the final vcf file, called genotypes_all

debug=0
ejecutar=function(texto){
    if(debug==1){print(texto)}else{system(texto)}}

command="vcf-concat inputs/merged_chrs/merged_chr1.vcf.gz inputs/merged_chrs/merged_chr2.vcf.gz inputs/merged_chrs/merged_chr3.vcf.gz inputs/merged_chrs/merged_chr4.vcf.gz inputs/merged_chrs/merged_chr5.vcf.gz inputs/merged_chrs/merged_chr6.vcf.gz inputs/merged_chrs/merged_chr7.vcf.gz inputs/merged_chrs/merged_chr8.vcf.gz inputs/merged_chrs/merged_chr9.vcf.gz inputs/merged_chrs/merged_chr10.vcf.gz inputs/merged_chrs/merged_chr11.vcf.gz inputs/merged_chrs/merged_chr12.vcf.gz inputs/merged_chrs/merged_chr13.vcf.gz inputs/merged_chrs/merged_chr14.vcf.gz inputs/merged_chrs/merged_chr15.vcf.gz inputs/merged_chrs/merged_chr16.vcf.gz inputs/merged_chrs/merged_chr17.vcf.gz inputs/merged_chrs/merged_chr18.vcf.gz inputs/merged_chrs/merged_chr19.vcf.gz inputs/merged_chrs/merged_chr20.vcf.gz inputs/merged_chrs/merged_chr21.vcf.gz inputs/merged_chrs/merged_chr22.vcf.gz | bgzip > inputs/processed/genotypes_all.vcf.gz"
ejecutar(command)

### Here apply MAF and missing values filters to get the file genotypes_all_MAF_10_missing_50.vcf

system("sh codes/process_genotype_scripts/master_codes/filters.sh")


### THis is the step I skipped before creating the 012 matrix: annotation, to get rs SNP IDs,
### I would follow the suggestion by finswimmer here: https://www.biostars.org/p/300169/


### From the last vcf file (filtered & annotated, or just fileterd), create the 012 matrix.
system("mkdir -p inputs/processed/EQTL_mapping/Genotypes")
system("sh codes/process_genotype_scripts/master_codes/create_matrix.sh")
