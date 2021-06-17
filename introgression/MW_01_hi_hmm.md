# Steps for HMM to examine MW_01_hi, the highly-introgressed MW8X-West samples

## Steps
* 1. Generate sample set list
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_01_hmm_samps.txt`
* 2. Generate VCF
* 3. Generate sample "map"
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_01_hmm_samp_map.txt`
## Files we need
* List of all samples to include in VCF
  * MW4X training samples
  * Gulf training samples
  * MW_01
  * MW_01_hi
* VCF of hiFst MWvsGULF SNPs
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMWhiFst.v3.disomic.CDS.vcf.gz`

## Output
* VCF
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW01_HMM_samps.vcf`
* Map of samples
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_01_hmm_samp_map_v2.txt`

## Generate Sample List and Map
* List = sample names for VCFtools
* Map = information for John
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_info_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt', sep = '')
res_info <- fread(res_info_file)

gulf_train_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/Gulf_libnames.txt', sep = '')
gulf_train_names <- fread(gulf_train_file, header = F)

mw_train_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/Midwest_libnames.txt', sep = '')
mw_train_names <- fread(mw_train_file, header = F)

meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

admix_res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/GW_50k_geobig.3.results.txt'
admix_res <- fread(admix_res_file)

### SET OUTPUT ###
out_dir <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/',
  'MW_analysis/', sep = '')

mw01_hmm_samp_out <- paste(out_dir, 'MW_01_hmm_samps.txt', sep = '')
#mw01_hmm_map_out <- paste(out_dir, 'MW_01_hmm_samp_map.txt', sep = '')
mw01_hmm_map_out <- paste(out_dir, 'MW_01_hmm_samp_map_v2.txt', sep = '')

#####################
mw01_reg_samps <- res_info[which(res_info$subgrp_v2 == 'MW_01'), samp_name]
mw01_hi_samps <- res_info[which(res_info$subgrp_v2 == 'MW_01_hi'), samp_name]

tot_samp_list <- c(gulf_train_names$V1, mw_train_names$V1, mw01_reg_samps, 
  mw01_hi_samps)

mw01_map <- data.table(GROUP = as.character(NA), SAMP = tot_samp_list)
mw01_map[which(mw01_map$SAMP %in% mw_train_names$V1), GROUP := 'MW_train']
mw01_map[which(mw01_map$SAMP %in% gulf_train_names$V1), GROUP := 'GULF_train']
mw01_map[which(mw01_map$SAMP %in% mw01_reg_samps), GROUP := 'MW01_reg']
mw01_map[which(mw01_map$SAMP %in% mw01_hi_samps), GROUP := 'MW01_hi']

mw_ancestry_vec <- c()
for(i in tot_samp_list){
  ri_ind <- which(admix_res$V1 == i)
  tmp_anc <- admix_res$V4[ri_ind]
  mw_ancestry_vec <- c(mw_ancestry_vec, tmp_anc)
}

gulf_ancestry_vec <- c()
for(i in tot_samp_list){
  ri_ind <- which(admix_res$V1 == i)
  tmp_anc <- admix_res$V3[ri_ind]
  gulf_ancestry_vec <- c(gulf_ancestry_vec, tmp_anc)
}

mw01_map[, MW_ancestry := mw_ancestry_vec]

fwrite(mw01_map[ , list(SAMP)], file = mw01_hmm_samp_out, col.names = F)
fwrite(mw01_map, file = mw01_hmm_map_out)

res_info[which(res_info$samp_name %in% mw01_hi_samps), full_admix_k3_MW]


```

## Generate VCF
```
bash
source activate bioinformatics_env

VCF_IN=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMWhiFst.v3.disomic.CDS.vcf.gz

POS_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_keep_pos.txt

SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_01_hmm_samps.txt

OUT_FILE=MW01_HMM_samps.vcf

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/

cd $OUT_DIR

vcftools --gzvcf $VCF_IN -c --keep $SAMP_FILE \
--positions $POS_FILE --recode --recode-INFO-all > $OUT_FILE

```
