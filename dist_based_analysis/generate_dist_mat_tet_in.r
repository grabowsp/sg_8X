# Script for generating a distance matrix from tetrasomic genotype VCF 
#   and using different types of genotypes; 
#   VCF must have header for this to work
#  
# Distance matrices are based on tetrasomic dosage genotypes:
#  RRRR (homozygous REF) = 2.0; AAAA (homozygous ALT) = 0.0
#  RRRA = 1.5, RRAA = 1.0, RAAA = 0.5
# If using disomic genotypes, then any HET is converted to 1.0
###############
# Arguments
# [1]: (character) filename of the VCF; must contain header
# [2]: (character) genotypes to be used for comparisons:
#			'diploid' = use disomic genotypes for all samples;
#			'tetraploid' = use tetrasomic genotypes for all samples;
#			'polyploid' = use ploidy-appropriate genotypes - 
#				4X use disomic, 6X/8X use tetrasomic
# [3]: (character) directory where distance matrix will be saved
###############

###############
# OUTPUT
# List with 4 elements:
# [['nSNPS']] = the number of SNPs used to calculate the distance
# [['euclidean_dist']] = euclidean distance matrix based on dosage genotypes
# [['manhattan_dist']] = manhattan distance matrix based on dosage genotypes
# [['n_NAs']] = the number of NA genotypes in each sample; NA's are removed
#			from distance calculations and normalized to nSNPs, so
#			these values can help determine if/how missing data
#			may be affecting results
###############

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

args = commandArgs(trailingOnly = TRUE)

### LOAD PACKAGES ###
library(data.table)
library(R.utils)

### LOAD INPUTS ###
vcf_in <- args[1]
#vcf_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/GW.natv2.ancestral.tet.100k.0001.vcf.gz'

vcf <- fread(vcf_in, skip = '#CHROM')
colnames(vcf)[1] <- 'CHROM'

samp_meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

comp_type <- args[2]

if(length(setdiff(c('diploid', 'tetraploid', 'polyploid'), comp_type)) >2){
  print('Error: comp_type args[2] not recognize; use diploid, tetraploid, or polyploid')
  quit(save = 'no')
}

### SET OUTPUTS ###
out_dir <- args[3]
#out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist'

# add slash to directory if not already there
out_dir_last_char <- rev(unlist(strsplit(out_dir, split = '')))[1]
if(out_dir_last_char != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

out_suffix <- paste('_', comp_type, '_DistMat.rds', sep = '')

#out_suffix <- '_tetrasomic_DistMat.rds'

file_pre <- gsub('.vcf', '', rev(unlist(strsplit(vcf_in, split = '/')))[1])

out_file <- paste(out_dir, file_pre, out_suffix, sep = '')
# out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/natv2_ancest_100k_tet_DistMat.rds'
### SET VARIABLES ###


##################
tet_libs <- intersect(colnames(vcf), 
  c(samp_meta[NQUIRE_PLOIDY == '4X', VCF_NAME], 'ANCESTRAL_TET_GENO', 
    'REF_TET_GENO'))
oct_libs <- intersect(colnames(vcf),
  samp_meta[NQUIRE_PLOIDY == '8X'|NQUIRE_PLOIDY == '6X', VCF_NAME])

#geno_df <- vcf[, c(10:ncol(vcf))]

oct_df <- as.data.frame(vcf[, c(oct_libs), with = F], stringsAsFactors = F)
tet_df <- as.data.frame(vcf[, c(tet_libs), with = F], stringsAsFactors = F)

for(i in oct_libs){
  oct_df[, i] <- unlist(lapply(strsplit(oct_df[, i], split = ':'), 
    function(x) x[[1]]))
}

for(j in tet_libs){
  tet_df[, j] <- unlist(lapply(strsplit(tet_df[, j], split = ':'),
    function(x) x[[1]])) 
}

oct_df[oct_df == './.'] <- NA
tet_df[tet_df == './.'] <- NA

geno_vec <- c('4/0', '3/1', '2/2', '1/3', '0/4')
oct_dosage_vec <- c('2', '1.5', '1', '0.5', '0')
tet_dosage_vec <- c('2', '1', '1', '1', '0')

if(comp_type == 'diploid'){
  for(i in seq(length(geno_vec))){
    oct_df[oct_df == geno_vec[i]] <- tet_dosage_vec[i]
    tet_df[tet_df == geno_vec[i]] <- tet_dosage_vec[i]   
  }
}

if(comp_type == 'tetraploid'){
  for(i in seq(length(geno_vec))){
    oct_df[oct_df == geno_vec[i]] <- oct_dosage_vec[i]
    tet_df[tet_df == geno_vec[i]] <- oct_dosage_vec[i]
  }
}

if(comp_type == 'polyploid'){
  for(i in seq(length(geno_vec))){
    oct_df[oct_df == geno_vec[i]] <- oct_dosage_vec[i]
    tet_df[tet_df == geno_vec[i]] <- tet_dosage_vec[i]
  }
}

for(i in seq(ncol(oct_df))){
  oct_df[, i] <- as.numeric(oct_df[, i])
}

for(i in seq(ncol(tet_df))){
  tet_df[, i] <- as.numeric(tet_df[, i])
}

geno_df <- cbind(oct_df, tet_df)

geno_mat <- matrix(unlist(geno_df), ncol = nrow(geno_df), 
  byrow = T)

rownames(geno_mat) <- colnames(geno_df)

dist_euc <- dist(geno_mat, diag = T, upper = T, method = 'euclidean')
dist_man <- dist(geno_mat, diag = T, upper = T, method = 'manhattan')

nSNPs <- nrow(vcf)

n_NAs <- apply(geno_mat, 1, function(x) sum(is.na(x)))

tot_list <- list()

tot_list[['nSNPs']] <- nSNPs
tot_list[['euclidean_dist']] <- dist_euc
tot_list[['manhattan_dist']] <- dist_man
tot_list[['n_NAs']] <- n_NAs

saveRDS(tot_list, file = out_file)

quit(save = 'no')

