# Steps for looking at genotypes in the potential inversion of Chr08N

## Overview
* Load Hi-Fst SNP genotypes for MW
* Find SNPs with high introgression and RDA scores
* Look at distribution of genotypes at these SNPs in MW
* Generate Hi-Fst SNPs for all samples (including those not in the paper)
* Look at distribution of genotypes in all samples

## Look at Candidate SNPs in MW 
```
### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
#vcf_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz'
#vcf <- fread(vcf_file)

# Atlantic results
atl_vcf_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_MWvATL.vcf'
atl_vcf_in <- read.table(atl_vcf_file, header = F, stringsAsFactors = F)
atl_vcf_head_tmp <- system(paste('grep CHR ', atl_vcf_file, sep = ''), 
  intern = T)
atl_vcf_head_2 <- sub('#CHROM', 'CHROM', atl_vcf_head_tmp)
atl_vcf_head_3 <- unlist(strsplit(atl_vcf_head_2, split = '\t'))
colnames(atl_vcf_in) <- atl_vcf_head_3
atl_vcf_in <- data.table(atl_vcf_in)

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

##########




```





