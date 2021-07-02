# Looking at genotypes in the potential inversion of Chr08N

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

vcf_func_file <- '/home/grabowsky/tools/workflows/sg_8X/introgression/vcf_analysis_functions.r'
source(vcf_func_file)

### INPUT DATA ###
#vcf_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz'
#vcf <- fread(vcf_file)

# Atlantic results
atl_vcf_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_MWvATL.vcf'
#atl_vcf_in <- read.table(atl_vcf_file, header = F, stringsAsFactors = F)
#atl_vcf_head_tmp <- system(paste('grep CHR ', atl_vcf_file, sep = ''), 
#  intern = T)
#atl_vcf_head_2 <- sub('#CHROM', 'CHROM', atl_vcf_head_tmp)
#atl_vcf_head_3 <- unlist(strsplit(atl_vcf_head_2, split = '\t'))
#colnames(atl_vcf_in) <- atl_vcf_head_3
#atl_vcf_in <- data.table(atl_vcf_in)
atl_vcf_in <- read_vcf(atl_vcf_file)

atl_fval_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/',
  'MW8X_all_Fval_ATL_into_MW_allsnps.txt', sep = '')
atl_fvals <- fread(atl_fval_file)

atl_rda_res_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/',
  'Results_Atlantic_Midwest_8x_4x_rda.csv', sep = '')
atl_rda_res <- fread(atl_rda_res_file)

# Pop structure Results file
res_tab_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_tab_file)

# metadata file
meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

### SET VARIABLES ###
sd_cut <- 3
##########

# Atlantic
atl_outlier_cut <- mean(atl_fvals$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(atl_fvals$F_subgrp_v_pop1, na.rm = T))
atl_outlier_inds <- which(atl_fvals$F_subgrp_v_pop1 > atl_outlier_cut)

atl_fvals[, SNP_NAME := paste(atl_fvals$CHR, atl_fvals$POS, sep = '_')]
atl_rda_overlap_inds <- which(atl_fvals$SNP_NAME %in% atl_rda_res$snp)

atl_snp_names <- atl_fvals[intersect(atl_rda_overlap_inds, atl_outlier_inds), 
  SNP_NAME]

atl_geno_1 <- vcf_SNP_genotypes(atl_vcf_in, snp_name = atl_snp_names[1])
atl_geno_2 <- vcf_SNP_genotypes(atl_vcf_in, snp_name = atl_snp_names[2])

atl_overlap_genos <- vcf_many_SNP_genotypes(atl_vcf_in, 
  snp_name_vec = atl_snp_names)

atl_vcf_samp_vec <- colnames(atl_vcf_in)[c(10:ncol(atl_vcf_in))]
atl_vcf_subgrp_vec <- rep(as.character(NA), times = length(atl_vcf_samp_vec))
atl_vcf_ploidy_vec <- rep(as.character(NA), times = length(atl_vcf_samp_vec))
for(asn in seq(length(atl_vcf_samp_vec))){
  tmp_res_ind <- which(res_tab$samp_name == atl_vcf_samp_vec[asn])
  atl_vcf_subgrp_vec[asn] <- res_tab$subgrp_v2[tmp_res_ind]
  atl_vcf_ploidy_vec[asn] <- res_tab$ploidy[tmp_res_ind]
}
names(atl_vcf_subgrp_vec) <- atl_vcf_samp_vec
names(atl_vcf_ploidy_vec) <- atl_vcf_samp_vec

# look at breakdown of genotype classes by subgrp
atl_snp1_subgrps <- get_portion_genotype_1SNP(atl_geno_1, 
  group_info = atl_vcf_subgrp_vec, return_percent = F)

atl_snp1_ploidy <- get_portion_genotype_1SNP(atl_geno_1,
  group_info = atl_vcf_ploidy_vec, return_percent = F)

atl_overlap_subgrps <- get_portion_genotypes(geno_tab = atl_overlap_genos,
  group_info = atl_vcf_subgrp_vec)

atl_overlap_ploidy <- get_portion_genotypes(geno_tab = atl_overlap_genos,
  group_info = atl_vcf_ploidy_vec)

atl_Chr08N_snps <- atl_snp_names[grep('Chr08N', atl_snp_names)]

# 34134827-21077373=13,057,454 = 13 Mb region

atl_chr08N_fval_info <- atl_fvals[intersect(
  intersect(atl_rda_overlap_inds, atl_outlier_inds), 
    which(atl_fvals$CHR == 'Chr08N')), ]
# 31078673-22589884=8,488,789 = 8.5 Mb region where frequency is essentially
#   ATL frequency

atl_Chr08N_snps_short <- atl_Chr08N_snps[3:21]

atl_chr08N_short_genos <- vcf_many_SNP_genotypes(atl_vcf_in,
  snp_name_vec = atl_Chr08N_snps_short)
atl_chr08N_short_subgrps <- get_portion_genotypes(
  geno_tab = atl_chr08N_short_genos, group_info = atl_vcf_subgrp_vec)
atl_chr08N_short_ploidy <- get_portion_genotypes(
  geno_tab = atl_chr08N_short_genos, group_info = atl_vcf_ploidy_vec)

atl_chr08N_short_genos[, which(atl_vcf_subgrp_vec == 'MW_01_hi'), with = F]

# Look at Chr03K and Chr03N peaks

atl_chr03K_fval_info <- atl_fvals[intersect(
  intersect(atl_rda_overlap_inds, atl_outlier_inds),
    which(atl_fvals$CHR == 'Chr03K')), ]
# 52026630-51871367 = 155,263 = 155kb
atl_Chr03K_snps <- atl_snp_names[grep('Chr03K', atl_snp_names)]
atl_chr03K_genos <- vcf_many_SNP_genotypes(atl_vcf_in,
  snp_name_vec = atl_Chr03K_snps)
atl_chr03K_subgrps <- get_portion_genotypes(
  geno_tab = atl_chr03K_genos, group_info = atl_vcf_subgrp_vec)
atl_chr03K_ploidy <- get_portion_genotypes(
  geno_tab = atl_chr03K_genos, group_info = atl_vcf_ploidy_vec)

atl_chr03K_genos[, which(atl_vcf_subgrp_vec == 'MW_01'), with = F]

atl_chr03N_fval_info <- atl_fvals[intersect(
  intersect(atl_rda_overlap_inds, atl_outlier_inds),
    which(atl_fvals$CHR == 'Chr03N')), ]
# 56452006-56281436 = 170,570 = 170kb
atl_Chr03N_snps <- atl_snp_names[grep('Chr03N', atl_snp_names)]
atl_chr03N_genos <- vcf_many_SNP_genotypes(atl_vcf_in,
  snp_name_vec = atl_Chr03N_snps)
atl_chr03N_subgrps <- get_portion_genotypes(
  geno_tab = atl_chr03N_genos, group_info = atl_vcf_subgrp_vec)
atl_chr03N_ploidy <- get_portion_genotypes(
  geno_tab = atl_chr03N_genos, group_info = atl_vcf_ploidy_vec)

atl_chr03N_genos[, which(atl_vcf_subgrp_vec == 'MW_01'), with = F]


