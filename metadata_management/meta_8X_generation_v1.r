# Generate first version of 8X-specific metadata

# Consolidate info from old metadata, coverage info, newest nQuire results, 
#   old nQuire results, old MNP results and the names in the VCFs

# on NERSC
# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

# on HA
# source activate R_analysis

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)
library(bit64)

### LOAD DATA ###
# John's metadata
samp_meta_file <- paste('/global/homes/g/grabowsp/data/switchgrass/',
  'reseq_metadata/', 'genotype.metadata.May2020.rds', sep = '')
samp_meta <- data.table(readRDS(samp_meta_file))

# Jason's metadata file
j_samp_metadata_file <- paste('/global/homes/g/grabowsp/data/switchgrass/',
  'reseq_metadata/PVDIV_Master_Metadata_File_7-1-20.txt', sep = '')
jsamp_meta <- fread(j_samp_metadata_file)

# Sujan's nQuire results for new libraries and previous 8X
nquire_res_file <- paste('/global/homes/g/grabowsp/data/switchgrass/', 
  'nquire_8X/', 'Ploidy_338Samples_nquire.txt', sep = '')
nquire_res <- fread(nquire_res_file)

# previous nQuire res file for 4X and previous 8X libs
old_nquire_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/',
  'sg_nquire/sg_CDS_nquire/sg_reseq_nQuire_results_summary_total.txt',
  sep = '')
old_nquire_res <- fread(old_nquire_file)

# Info about coverage for old and new libraries
coverage_res_file <- paste('/global/homes/g/grabowsp/data/switchgrass/', 
  'nquire_8X/', 'Pvirg_1070G_reads_mapped.txt', sep = '')
coverage_res <- fread(coverage_res_file)

# Old ploidy calls using nQuire and MNPs
old_ploidy_file <- paste('/global/cscratch1/sd/grabowsp/sg_ploidy/', 
  'ploidy_calling/', 'sg_ploidy_results_v3.0.txt', sep = '')
old_ploidy <- fread(old_ploidy_file)

vcf_samp_name_file <- paste('/global/cscratch1/sd/grabowsp/sg_8X_scratch/', 
  'orig_sujan_files/', 'tot_VCF_samp_ID.txt', sep = '')
vcf_samps_raw <- unlist(fread(vcf_samp_name_file))
vcf_samps <- vcf_samps_raw[c(10:length(vcf_samps_raw))]

### SET OUTPUT ###

out_file <- paste('/global/homes/g/grabowsp/data/switchgrass/metadata_8X/', 
  'sg_8X_metadata_v1.0.csv', sep = '') 
#################

### Start metadata using coverage info from Sujan - these are the libraries
#     that he is using to generate the VCFs
new_meta <- data.table(LIB = coverage_res$LIB, PLANT_ID = as.character(NA), 
  ALIAS = as.character(NA),
  GOOD_PAIRED_READS = coverage_res$Properly_paired_reads, 
  APPROX_BP_MAPPED = coverage_res$Approx_bp_mapped, 
  APPROX_COVERAGE = coverage_res$Approx_genome_coverage_mapped,
  COVERAGE_CATEGORY = coverage_res$Category)

# Add AP13 as reference - John has this in his metadata
new_meta <- add_row(new_meta, LIB = 'REF', PLANT_ID = 'AP13', ALIAS = NA, 
  GOOD_PAIRED_READS = NA, APPROX_BP_MAPPED = NA, APPROX_COVERAGE = NA,
  COVERAGE_CATEGORY = NA)

### Add PLANT_ID from old libraries using old metadata
for(i in which(new_meta$LIB %in% samp_meta$LIBRARY)){
  tmp_ind <- which(samp_meta$LIBRARY == new_meta$LIB[i])
  new_meta[i, PLANT_ID := samp_meta$PLANT_ID[tmp_ind]]
}

### add PLANT_ID for new libraries using nQuire results from Sujan
for(j in seq(nrow(nquire_res))){
  tmp_ind <- which(new_meta$LIB == nquire_res$LIB[j])
  tmp_ID <- nquire_res$Sample[j]
  if(is.na(new_meta[tmp_ind, PLANT_ID])){
    if(gsub('_', '.', tmp_ID) %in% samp_meta$PLANT_ID){
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

### Add Sujans nQuire results
nquire_ploidy_cols <- c('nquire_cds_20', 'nquire_full_20')
nquire_ploidy_cols <- toupper(nquire_ploidy_cols)
new_meta[, (nquire_ploidy_cols) := as.character(NA)]
nquire_res_extra_cols <- 'OLD_META_PLOIDY'
# Don't need to use NQUIRE_CATEGORY because is the same as COVERAGE_CATEGORY
new_meta[, (nquire_res_extra_cols) := as.character(NA)]

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
  new_meta[nm_ind, OLD_META_PLOIDY := nquire_res$ploidy_metadat[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
}

### Add previous nQuire results
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

### CONSOLIDATE DUPLICATED PLANT_ID
dup_libs <- new_meta$PLANT_ID[which(duplicated(new_meta$PLANT_ID))]

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

### Add info from Jason's metadata

old_meta_char_cols <- c('ACC', 'COLLECTION_TYPE', 'COLLECTION_METHOD', 
  'SOURCE_ID_1', 'SOURCE_ID_2', 'COUNTRY', 'STATE', 'COUNTY',
  'NOTE_LATLONG', 'LOCALITY', 'TAXON', 'COLL_DATE',
  'PLD_SNP')

old_meta_num_cols <- c('LATITUDE', 'LONGITUDE', 'ELEVATION')

meta_3[, (old_meta_char_cols) := as.character(NA)]
meta_3[, (old_meta_num_cols) := as.numeric(NA)]
meta_3[, NAT_PAPER_LIB := as.character(NA)]

for(msn in intersect(meta_3$PLANT_ID, jsamp_meta$PLANT_ID)){
  nm_ind <- which(meta_3$PLANT_ID == msn)
  old_ind <- which(jsamp_meta$PLANT_ID == msn)
  tmp_char_vec <- jsamp_meta[old_ind, ..old_meta_char_cols]
  tmp_num_vec <- jsamp_meta[old_ind, ..old_meta_num_cols]
  meta_3[nm_ind, (old_meta_char_cols) := tmp_char_vec]
  meta_3[nm_ind, (old_meta_num_cols) := tmp_num_vec]
  meta_3[nm_ind, NAT_PAPER_LIB := jsamp_meta$LIB_CLIMATE[old_ind]]
}

meta_3[PLANT_ID == 'AP13', NAT_PAPER_LIB := 'Y']

### Add PLANTED_2018 and LIB_NOTES (I prefer John's formatting for this)
jl_meta_cols <- c('PLANTED_2018', 'LIB_NOTES')
meta_3[, (jl_meta_cols) := as.character(NA)]

for(msn in intersect(meta_3$PLANT_ID, samp_meta$PLANT_ID)){
  nm_ind <- which(meta_3$PLANT_ID == msn)
  old_ind <- which(samp_meta$PLANT_ID == msn)
  tmp_jl_vals <- samp_meta[old_ind, ..jl_meta_cols]
  meta_3[nm_ind, (jl_meta_cols) := tmp_jl_vals]
}

### Add info from previous MNP work
meta_3[, MNP_CDS_PLOIDY := as.character(NA)]
meta_3[, MNP_GENIC_PLOIDY := as.character(NA)]
meta_3[, OLD_TOT_PLOIDY := as.character(NA)]
meta_3[, OLD_TOT_PLOIDY_CONFIDENCE := as.numeric(NA)]

for(msn in intersect(meta_3$PLANT_ID, old_ploidy$PLANT_ID)){
  nm_ind <- which(meta_3$PLANT_ID == msn)
  old_ind <- which(old_ploidy$PLANT_ID == msn)
  meta_3[nm_ind, MNP_CDS_PLOIDY := old_ploidy$cds_mnp_ploidy[old_ind]]
  meta_3[nm_ind, MNP_GENIC_PLOIDY := old_ploidy$genic_mnp_ploidy[old_ind]]
  meta_3[nm_ind, OLD_TOT_PLOIDY := old_ploidy$total_ploidy_2[old_ind]]
  meta_3[nm_ind, OLD_TOT_PLOIDY_CONFIDENCE := 
    old_ploidy$tot_ploid_confidence[old_ind]]
}

### Edit the COLLECTION_TYPE column for new LIBs and based on notes
ct_tmp_ind_1 <- which(meta_3$COLLECTION_TYPE == 'Unknown')
meta_3[ct_tmp_ind_1, COLLECTION_TYPE := 'Contaminant']
meta_3[which(meta_3$PLANT_ID == 'AP13'), COLLECTION_TYPE := 'Reference']
mx_inds <- which(is.na(meta_3$COLLECTION_TYPE))
meta_3[mx_inds, COLLECTION_TYPE := 'Mexican Samples']

# Add names used in VCF
meta_3[, VCF_NAME = as.character(NA)]
mete_3[, MINOR_LIB_VCF_NAME = as.character(NA)]
for(vsamp in vcf_samps){
  if(vsamp %in% nquire_res$Sample){
    test_lib <- nquire_res$LIB[which(nquire_res$Sample == vsamp)][1]
    if(test_lib %in% meta_3$MINOR_LIB){
      tmp_ind <- which(meta_3$MINOR_LIB == test_lib)
      meta_3[tmp_ind, MINOR_LIB_VCF_NAME := vsamp]
    } else{
      tmp_ind <- which(meta_3$LIB == test_lib)
      meta_3[tmp_ind, VCF_NAME := vsamp]
    }
  } else if(vsamp %in% meta_3$PLANT_ID){
    tmp_ind <- which(meta_3$PLANT_ID == vsamp)
    meta_3[tmp_ind, VCF_NAME := vsamp]
  } else if(vsamp %in% meta_3$ALIAS){
    tmp_ind <- which(meta_3$ALIAS == vsamp)
    meta_3[tmp_ind, VCF_NAME := vsamp]
  }
}

### Fix one VCF sample name that's not in the PLANT_ID or ALIAS line
missing_name <- setdiff(vcf_samps, c(meta_3$VCF_NAME, 
  meta_3$MINOR_LIB_VCF_NAME))
# [1] "9001-3_BN389-69S"
tmp_ind <- grep('9001-3', meta_3$PLANT_ID)
meta_3[tmp_ind, VCF_NAME := missing_name]
meta_3[tmp_ind, ALIAS := missing_name]

### Re-order columns
new_name_ord <- c('PLANT_ID', 'LIB', 'ALIAS', 'VCF_NAME', 'ACC')
new_info_ord <- c('NAT_PAPER_LIB', 'PLANTED_2018', 'COVERAGE_CATEGORY', 
  'COLLECTION_TYPE', 'COLLECTION_METHOD', 'TAXON', 'SOURCE_ID_1', 
  'SOURCE_ID_2', 'LIB_NOTES', 'COLL_DATE')
new_geo_ord <- c('LATITUDE', 'LONGITUDE', 'ELEVATION', 'NOTE_LATLONG', 
  'COUNTRY', 'STATE', 'COUNTY', 'LOCALITY')
new_stats_ord <- c('APPROX_COVERAGE', 'APPROX_BP_MAPPED', 'GOOD_PAIRED_READS')
new_minlib_ord <- c('MINOR_LIB', 'MINOR_LIB_COV', 'MINOR_LIB_VCF_NAME')
new_ploidy_ord <- c('OLD_META_PLOIDY', 'PLD_SNP', 'NQUIRE_CDS_20', 
  'NQUIRE_FULL_20', 'NQUIRE_NSNPS_20', 'NQUIRE_PLOIDY', 'NQUIRE_CONFIDENCE', 
  'MINOR_LIB_NQUIRE_CDS_20', 'MNP_CDS_PLOIDY', 'MNP_GENIC_PLOIDY', 
  'OLD_TOT_PLOIDY', 'OLD_TOT_PLOIDY_CONFIDENCE')

new_ord_tot <- c(new_name_ord, new_info_ord, new_geo_ord, new_stats_ord, 
  new_minlib_ord, new_ploidy_ord)

meta_4 <- meta_3[, ..new_ord_tot]

fwrite(meta_4, out_file, na = 'NA')

#meta_4[, .N, by = COLLECTION_TYPE]
#      COLLECTION_TYPE   N
#1:    Mexican Samples  42
#2: Natural Collection 821
#3: Breeding Selection 121
#4:           Cultivar  73
#5:        Contaminant   1
#6:          Reference   1

# meta_4[, .N, by = .(COLLECTION_TYPE, NQUIRE_PLOIDY)]
#       COLLECTION_TYPE NQUIRE_PLOIDY   N
# 1:    Mexican Samples            4X  40
# 2:    Mexican Samples            8X   2
# 3: Natural Collection            6X   5
# 4: Breeding Selection            8X  14
# 5:           Cultivar            8X  29
# 6: Natural Collection            8X 204
# 7:           Cultivar            4X  43
# 8: Breeding Selection            6X   2
# 9: Natural Collection            4X 612
#10: Breeding Selection            4X 105
#11:           Cultivar            6X   1
#12:        Contaminant            8X   1
#13:          Reference          <NA>   1




quit(save = 'no')


