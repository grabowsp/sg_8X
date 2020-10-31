# Generate Consolidated File Containing nQuire results and ploidy designation
#  for all samples

# on NERSC
# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

# on HA
# source activate R_analysis

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)
library(patchwork)

### LOAD DATA ###
# on HA
#data_dir <- '/home/t4c1/WORK/grabowsk/data/switchgrass/nquire_8X/'
# on NERSC
data_dir <- '/global/homes/g/grabowsp/data/switchgrass/nquire_8X/'

nquire_res_file_short <- 'Ploidy_338Samples_nquire.txt'
nquire_res_file <- paste(data_dir, nquire_res_file_short, sep = '')
nquire_res <- fread(nquire_res_file)
#nquire_res <- read.table(nquire_res_file, header = T, sep = '\t', 
#  stringsAsFactors = F)

# old nQuire res file on Cori
old_nquire_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total.txt',
  sep = '')
old_nquire_res <- fread(old_nquire_file)

coverage_res_file_short <- 'Pvirg_1070G_reads_mapped.txt'
coverage_res_file <- paste(data_dir, coverage_res_file_short, sep = '')
coverage_res <- fread(coverage_res_file)
#coverage_res <- read.table(nquire_res_file, header = T, sep = '\t',
#  stringsAsFactors = F)

# on NERSC
samp_metadata_file <- paste('/global/homes/g/grabowsp/data/switchgrass/', 
  'reseq_metadata/Reseq_Metadata_Sept_2019_Edited_for_R.tsv', sep = '')
# on HA
# samp_metadata_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_9_3_2019_tmp_for_R.txt'

samp_metadata <- fread(samp_metadata_file)

### SET OUTPUT ###

out_file <- '/global/homes/g/grabowsp/data/switchgrass/nquire_8X/consolidated_nQuire_results.txt'
#################

new_meta <- data.table(LIB = coverage_res$LIB, PLANT_ID = as.character(NA), 
  ALIAS = as.character(NA),
  GOOD_PAIRED_READS = coverage_res$Properly_paired_reads, 
  APPROX_BP_MAPPED = coverage_res$Approx_bp_mapped, 
  APPROX_COVERAGE = coverage_res$Approx_genome_coverage_mapped,
  COVERAGE_CATEGORY = coverage_res$Category)

for(i in which(new_meta$LIB %in% samp_metadata$LIBRARY)){
  tmp_ind <- which(samp_metadata$LIBRARY == new_meta$LIB[i])
  new_meta[i, PLANT_ID := samp_metadata$PLANT_ID[tmp_ind]]
#  new_meta$PLANT_ID[i] <- samp_metadata$PLANT_ID[tmp_ind]
}

for(j in seq(nrow(nquire_res))){
  tmp_ind <- which(new_meta$LIB == nquire_res$LIB[j])
  tmp_ID <- nquire_res$Sample[j]
  if(is.na(new_meta[tmp_ind, PLANT_ID])){
    if(gsub('_', '.', tmp_ID) %in% samp_metadata$PLANT_ID){
      new_meta[tmp_ind, PLANT_ID := gsub('_', '.', tmp_ID)]
      new_meta[tmp_ind, ALIAS := tmp_ID]
    } else{
      new_meta[tmp_ind, PLANT_ID := tmp_ID]
    }
  } else if((tmp_ID != new_meta[tmp_ind, PLANT_ID]) & 
             is.na(new_meta[tmp_ind, ALIAS])){
    new_meta[tmp_ind, ALIAS := tmp_ID]
  }  
}

nquire_ploidy_cols <- c('nquire_cds_20', 'nquire_full_20')

nquire_ploidy_cols <- toupper(nquire_ploidy_cols)

new_meta[, (nquire_ploidy_cols) := as.character(NA)]

nquire_res_extra_cols <- c('NQUIRE_CATEGORY', 'OLD_META_PLOIDY')

new_meta[, (nquire_res_extra_cols) := as.character(NA)]

#

nquire_res_libs <- unique(nquire_res$LIB)

nq_cds <- which(nquire_res$region == 'cds')
nq_full <- which(nquire_res$region == 'genome')
nq_20 <- which(nquire_res$depth == 20)

for(nrl in nquire_res_libs){
  nm_ind <- which(new_meta$LIB == nrl)
  new_meta[nm_ind, NQUIRE_CDS_20 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
  new_meta[nm_ind, NQUIRE_FULL_20 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_full, nq_20))]]
  new_meta[nm_ind, NQUIRE_CATEGORY := nquire_res$Category[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
  new_meta[nm_ind, OLD_META_PLOIDY := nquire_res$ploidy_metadat[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
}

old_20_cols <- c('d_dip_portion_20', 'd_tri_portion_20', 'd_tet_portion_20')
old_20_calls <- paste((apply(old_nquire_res[, ..old_20_cols],  
  1, which.min) + 1) * 2, 'X', sep = '')

old_nquire_res[, ploidy_20 := old_20_calls]

new_meta[, NQUIRE_NSNPS_20 := as.numeric(NA)]

for(onl in seq(nrow(old_nquire_res))){
  new_ind <- which(new_meta$LIB == old_nquire_res$samp[onl])
  new_meta[new_ind, NQUIRE_CDS_20 := old_nquire_res$ploidy_20[onl]]
  new_meta[new_ind, NQUIRE_NSNPS_20 := old_nquire_res$nSNPS_20[onl]]
}

new_meta[is.na(OLD_META_PLOIDY), OLD_META_PLOIDY := '4x']

# change all ploidy to have uppercase X
nquire_ploid_cols <- c('NQUIRE_CDS_20', 'NQUIRE_FULL_20', 'OLD_META_PLOIDY')

for(npc in nquire_ploid_cols){
  tmp_p <- toupper(unlist(new_meta[, npc, with = F]))
  new_meta[, (npc) := tmp_p]
}

# Based on plot, want at least 20 coverage to ensure that have 500K SNPs
## however, can have 500k SNPs at coverage less than 20

dup_libs <- new_meta$PLANT_ID[which(duplicated(new_meta$PLANT_ID))]

#####
# Generate metadata file with PLANT_ID as unique column

dup_minor_lib <- c()
for(dup in dup_libs){
  dup_minor_lib <- c(dup_minor_lib, new_meta[PLANT_ID == dup & 
    APPROX_COVERAGE == min(new_meta[PLANT_ID == dup, APPROX_COVERAGE]), LIB])
}

dup_minor_inds <- c()
for(dmi in dup_minor_lib){
  dup_minor_inds <- c(dup_minor_inds, which(new_meta$LIB == dmi))
}

meta_3 <- new_meta[-dup_minor_inds, ]
meta_3[, MINOR_LIB := as.character(NA)]
meta_3[, MINOR_LIB_COV := as.numeric(NA)]
meta_3[, MINOR_LIB_NQUIRE_CDS_20 := as.character(NA)]

for(dmi in dup_minor_inds){
  meta_3[PLANT_ID == new_meta$PLANT_ID[dmi], MINOR_LIB := new_meta$LIB[dmi]]
  meta_3[PLANT_ID == new_meta$PLANT_ID[dmi], 
    MINOR_LIB_COV := new_meta$APPROX_COVERAGE[dmi]]
  meta_3[PLANT_ID == new_meta$PLANT_ID[dmi], 
    MINOR_LIB_NQUIRE_CDS_20 := new_meta$NQUIRE_CDS_20[dmi]]
}

meta_3[, NQUIRE_PLOIDY := NQUIRE_CDS_20]

# add confidence measure for nQuire ploidy
## 1 = high, 2 = coverage below 20, 3 = CDS and FULL don't match, 4 = 2&3
meta_3[, NQUIRE_CONFIDENCE := as.numeric(NA)]
meta_3[APPROX_COVERAGE >= 20, NQUIRE_CONFIDENCE := 1]
meta_3[APPROX_COVERAGE < 20, NQUIRE_CONFIDENCE := 2]
meta_3[NQUIRE_CDS_20 != NQUIRE_FULL_20, NQUIRE_CONFIDENCE := 3]
meta_3[APPROX_COVERAGE < 20 & NQUIRE_CDS_20 != NQUIRE_FULL_20, 
  NQUIRE_CONFIDENCE := 4]

fwrite(meta_3, out_file, sep = '\t', na = 'NA', quote = F)

quit(save = 'no')


