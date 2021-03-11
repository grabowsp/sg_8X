# Script for calculating heterozygosity from tetrasomic genotype VCFs
#   Calculates both raw and adjusted (corrected) heterozygosity for all samples

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

args = commandArgs(trailingOnly = T)
rundir_args <- commandArgs(trailingOnly = F)

file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)

gen.basename <- sub('/diversity_measur', '', script.basename)

### LOAD PACKAGES ###

library(data.table)
library(R.utils)

# gen_funct_file <- '/global/homes/g/grabowsp/tools/sg_8X/general_r_tools/general_functions.r'
gen_funct_file <- file.path(gen.basename, 
  'general_r_tools/general_functions.r')
source(gen_funct_file)

### LOAD INPUTS ###

vcf_file_in <- args[1]
#vcf_file_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_tet_vcfs/Chr01K.tetrasomic.CDS.natv2.vcf_01'
#vcf_file_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/GW.natv2.ancestral.tet.100k.0001.vcf.gz'

vcf_header_file <- args[2]
#vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_tet_vcfs/CDS.tetrasomic.natv2.vcf.sampheader.txt'

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

first_line <- readLines(vcf_file_in, n = 1)
if(substr(first_line, 1,2) == '##'){
  vcf <- fread(vcf_file_in, skip = '#CHROM')
  colnames(vcf)[1] <- 'CHROM'
} else{
  vcf <- fread(vcf_file_in)
  colnames(vcf) <- vcf_header
}

### SET OUTPUTs ###

out_dir <- args[3]
out_dir <- add_slash(out_dir)

out_pre <- args[4]
# out_pre <- 'Chr01.geobig.CDS'

out_total <- paste(out_dir, out_pre, '.HETvals.rds', sep = '')

### SET VARIABLES ###
geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
oct_dosage_vec <- c('2', '1.5', '1', '0.5', '0')
#tet_dosage_vec <- c('2', '1', '1', '1', '0')

################

first_samp <- which(colnames(vcf) == 'FORMAT')+1

samp_df <- as.data.frame(vcf[, c(first_samp:ncol(vcf)), with = F], 
  stringsAsFactors = F)

for(i in seq(ncol(samp_df))){
  samp_df[, i] <- unlist(lapply(strsplit(samp_df[, i], split = ':'), 
    function(x) x[[1]]))
} 

samp_df[samp_df == './.'] <- NA

for(j in seq(length(geno_vec))){
  samp_df[samp_df == geno_vec[j]] <- oct_dosage_vec[j]
}

n_hom <- apply(samp_df, 2, function(x) sum(x == '2' | x == '0', na.rm = T))
n_het <- apply(samp_df, 2, function(x) sum(x == '0.5' | x == '1' | x == '1.5', 
  na.rm = T))

per_het <- n_het/(n_het+n_hom)

n_dose1 <- apply(samp_df, 2, function(x) sum(x == '0.5' | x == '1.5', 
  na.rm = T))
n_dose2 <- apply(samp_df, 2, function(x) sum(x == '1', na.rm = T))

corrected_het_per <- ((n_dose2*2/3)+(n_dose1*0.5)) / (n_het+n_hom)

compiled_list <- list()
compiled_list[['per_het']] <- per_het
compiled_list[['corrected_per_het']] <- corrected_het_per
compiled_list[['doseage_1v2_ratio']] <- n_dose1/n_dose2
compiled_list[['N_SNPs']] <- n_het+n_hom

saveRDS(compiled_list, out_total)

quit(save = 'no')


