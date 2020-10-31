# Initial exploration of nquire results

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

nquire_ploidy_cols <- c('nquire_cds_20', 'nquire_full_20', 'nquire_cds_35', 
  'nquire_full_35', 'nquire_cds_50', 'nquire_full_50')

nquire_ploidy_cols <- toupper(nquire_ploidy_cols)

new_meta[, (nquire_ploidy_cols) := as.character(NA)]

nquire_res_extra_cols <- c('NQUIRE_CATEGORY', 'OLD_META_PLOIDY')

new_meta[, (nquire_res_extra_cols) := as.character(NA)]

#

nquire_res_libs <- unique(nquire_res$LIB)

nq_cds <- which(nquire_res$region == 'cds')
nq_full <- which(nquire_res$region == 'genome')

nq_20 <- which(nquire_res$depth == 20)
nq_35 <- which(nquire_res$depth == 35)
nq_50 <- which(nquire_res$depth == 50)

for(nrl in nquire_res_libs){
  nm_ind <- which(new_meta$LIB == nrl)
  new_meta[nm_ind, NQUIRE_CDS_20 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
  new_meta[nm_ind, NQUIRE_FULL_20 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_full, nq_20))]]
  new_meta[nm_ind, NQUIRE_CDS_35 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_35))]]
  new_meta[nm_ind, NQUIRE_FULL_35 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_full, nq_35))]]
  new_meta[nm_ind, NQUIRE_CDS_50 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_50))]]
  new_meta[nm_ind, NQUIRE_FULL_50 := nquire_res$Ploidy_nQuire[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_full, nq_50))]]
  new_meta[nm_ind, NQUIRE_CATEGORY := nquire_res$Category[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
  new_meta[nm_ind, OLD_META_PLOIDY := nquire_res$ploidy_metadat[
    intersect(which(nquire_res$LIB == nrl), intersect(nq_cds, nq_20))]]
}



old_20_cols <- c('d_dip_portion_20', 'd_tri_portion_20', 'd_tet_portion_20')
old_20_calls <- paste((apply(old_nquire_res[, ..old_20_cols],  
  1, which.min) + 1) * 2, 'X', sep = '')

old_50_cols <- sub('20', '50', old_20_cols)
old_50_calls <- paste((apply(old_nquire_res[, ..old_50_cols],
  1, which.min) + 1) * 2, 'X', sep = '')

old_nquire_res[, ploidy_20 := old_20_calls]
old_nquire_res[, ploidy_50 := old_50_calls]

new_meta[, NQUIRE_NSNPS_20 := as.numeric(NA)]
new_meta[, NQUIRE_NSNPS_50 := as.numeric(NA)]

for(onl in seq(nrow(old_nquire_res))){
  new_ind <- which(new_meta$LIB == old_nquire_res$samp[onl])
  new_meta[new_ind, NQUIRE_CDS_20 := old_nquire_res$ploidy_20[onl]]
  new_meta[new_ind, NQUIRE_CDS_50 := old_nquire_res$ploidy_50[onl]]
  new_meta[new_ind, NQUIRE_NSNPS_20 := old_nquire_res$nSNPS_20[onl]]
  new_meta[new_ind, NQUIRE_NSNPS_50 := old_nquire_res$nSNPS_50[onl]]
}

# change all ploidy to have uppercase X
nquire_ploid_cols <- colnames(new_meta)[grep('_CDS_|_FULL_', 
  colnames(new_meta))]

for(npc in nquire_ploid_cols){
  tmp_p <- toupper(unlist(new_meta[, npc, with = F]))
  new_meta[, (npc) := tmp_p]
}

new_meta[, OLD_META_PLOIDY := toupper(OLD_META_PLOIDY)]
new_meta[is.na(OLD_META_PLOIDY), OLD_META_PLOIDY := '4X']

# I'm curious about relationship between coverage and nSNPS;
## In previous analysis, nQuire was very relible with >500k SNPs
gg_nsnps <- ggplot(new_meta[NQUIRE_NSNPS_20 >= 0,], 
  aes(x = APPROX_COVERAGE, y = NQUIRE_NSNPS_20)) +
  geom_point() +
  ggtitle('Seq coverage vs nSNPs with 20+ coverage')

cov_v_nsnp_pdf <- paste('/global/cscratch1/sd/grabowsp/sg_8X_scratch/', 
  'nquire_scratch/', 'coverage_vs_nSNPS20.pdf', sep = '')

pdf(file = cov_v_nsnp_pdf, width = 5, height = 5)
gg_nsnps
dev.off()

# Based on plot, want at least 20 coverage to ensure that have 500K SNPs
## however, can have 500k SNPs at coverage less than 20

dup_libs <- new_meta$PLANT_ID[which(duplicated(new_meta$PLANT_ID))]

new_meta[PLANT_ID %in% dup_libs,]
# all duplicated PLANT_IDs have new LIBs with coverage > 20

low_cov_libs <- new_meta[APPROX_COVERAGE < 20, LIB]
# 50 total libraries with coverage below 20

low_cov_libs_2 <- setdiff(low_cov_libs, new_meta[PLANT_ID %in% dup_libs, LIB])
# 38 Plants with coverage below 20

low_cov_libs_3 <- setdiff(low_cov_libs_2, 
  new_meta[NQUIRE_NSNPS_20 >= 5e5, LIB])
# 19 Plants with coverage below 19 and nSNPs at 20depth < 500k or unknown

new_meta[LIB %in% low_cov_libs_3, ]
# 8 MX samples
# 1 new sample - 8X
# 8 LIBs used in Nature paper: 7 assigned 4X, 1 assigned 6X using nQuire
##  IJBQ = J593.A, which is assigned 6X by nQuire. Not sure what to do
##         about that one
# 2 previous 8X libs:
##  IIDE = J581.A, there are two other J581 samples, BUT should retain IIDE
##          for pop structure analysis: IIDE looks 8X, and there is one other
##          8X and one 4X, so should use IIDE to make sure the other 8X is
##          correct pop (though can't corroborate the 4X if from different
##          genepool.
##  IIDF = J045.A, there is one other J045, also with lower coverage, but with
##         enough nSNPs_20 to safely assign as 4X. Not sure what to do about
##         this one. May want to try resequencing again?

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


# CONTINUE FROM HERE

meta_3[, CDS20_NUM := as.numeric(sub('X', '', NQUIRE_CDS_20))]
meta_3[, OLD_AGREE := NQUIRE_CDS_20 == OLD_META_PLOIDY]
meta_3[, CDS_FULL_AGREE := NQUIRE_CDS_20 == NQUIRE_FULL_20]

gg_new_old_agree <- ggplot(meta_3[COVERAGE_CATEGORY == 'new_164g', ], 
  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = OLD_AGREE)) +
#  geom_point() + 
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) + 
  ylab('Ploidy from nQuire') + 
  ggtitle(paste('nQuire Ploidy Calls for 164 new samples\nand if ', 
    'agrees with old metadata', sep = ''))

gg_new_full_agree <- ggplot(meta_3[COVERAGE_CATEGORY == 'new_164g', ],
  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = CDS_FULL_AGREE)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) +
  ylab('Ploidy from nQuire') +
  ggtitle(paste('nQuire Ploidy Calls for 164 new samples\nand if ',
    'CDS and FULL results agree', sep = ''))

gg_mx_old_agree <- ggplot(meta_3[COVERAGE_CATEGORY == 'mexican', ],
  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = OLD_AGREE)) +
#  geom_point() + 
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) +
  ylab('Ploidy from nQuire') +
  ggtitle(paste('nQuire Ploidy Calls for 42 Mexican samples\nand if ',
    'agrees with old metadata', sep = ''))

gg_mx_full_agree <- ggplot(meta_3[COVERAGE_CATEGORY == 'mexican', ],
  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = CDS_FULL_AGREE)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) +
  ylab('Ploidy from nQuire') +
  ggtitle(paste('nQuire Ploidy Calls for 42 Mexican samples\nand if ',
    'CDS and FULL results agree', sep = ''))

meta_3[, OLD_AGREE_2 := as.character(OLD_AGREE)]
meta_3[is.na(OLD_AGREE), OLD_AGREE_2 := 'MISSING']
meta_3[, FULL_AGREE_2 := as.character(CDS_FULL_AGREE)]
meta_3[is.na(CDS_FULL_AGREE), FULL_AGREE_2 := 'MISSING']

gg_all_old_agree <- ggplot(meta_3,
#  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = OLD_AGREE_2, 
#    shape = COVERAGE_CATEGORY)) +
  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = COVERAGE_CATEGORY)) + 
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) +
  ylab('Ploidy from nQuire') +
  ggtitle(paste('nQuire Ploidy Calls for 1058 samples\nby sample ',
    'category', sep = ''))

gg_all_full_agree <- ggplot(meta_3,
#  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = FULL_AGREE_2, 
#    shape = COVERAGE_CATEGORY)) +
  aes(x = APPROX_COVERAGE, y = CDS20_NUM, color = COVERAGE_CATEGORY)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) +
  ylab('Ploidy from nQuire') +
  ggtitle(paste('nQuire Ploidy Calls for 1058 samples\nand if ',
    'CDS and FULL results agree', sep = ''))

########

new164_ploidy_v1.0_pdf <- paste('/global/cscratch1/sd/grabowsp/sg_8X_scratch/',
  'nquire_scratch/', 'new164_ploidy_vs_coverage_v1.pdf', sep = '')

pdf(file = new164_ploidy_v1.0_pdf, width = 5.5, height = 10)
gg_new_old_agree / gg_new_full_agree
dev.off()

full_ploidy_v1_pdf <- paste('/global/cscratch1/sd/grabowsp/sg_8X_scratch/',
  'nquire_scratch/', 'all_ploidy_vs_coverage_v1.pdf', sep = '')

pdf(file = full_ploidy_v1_pdf, width = 12, height = 15)
(gg_new_old_agree + gg_new_full_agree) / (gg_mx_old_agree + 
  gg_mx_full_agree) / (gg_all_old_agree)
dev.off()

#######

# Barplots

meta_tab_1 <- meta_3[, .N, by = .(NQUIRE_CDS_20, COVERAGE_CATEGORY)]

gg_bar_1 <- ggplot(meta_tab_1, aes(x = NQUIRE_CDS_20, y = N, 
  fill = COVERAGE_CATEGORY)) +
  geom_bar(stat = 'identity') + 
  ggtitle('nQuire-based ploidy for all 1058 samples\nby category')

meta_tab_2 <- meta_3[, .N, 
  by = .(NQUIRE_CDS_20, COVERAGE_CATEGORY, OLD_AGREE_2)]

gg_bar_4X <- ggplot(meta_tab_2[NQUIRE_CDS_20 == '4X', ], 
  aes(x = COVERAGE_CATEGORY, y = N, fill = OLD_AGREE_2)) +
  geom_bar(stat = 'identity') + 
  ggtitle(paste('4X samples by category\nand agreement', 
    ' with old metadata', sep = '')) +
  guides(fill=guide_legend(title='Expected?')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gg_bar_6X <- ggplot(meta_tab_2[NQUIRE_CDS_20 == '6X', ], 
  aes(x = COVERAGE_CATEGORY, y = N, fill = OLD_AGREE_2)) +
  geom_bar(stat = 'identity') +
  ggtitle(paste('6X samples by category\nand agreement', 
    ' with old metadata', sep = '')) +
  guides(fill=guide_legend(title='Expected?')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gg_bar_8X <- ggplot(meta_tab_2[NQUIRE_CDS_20 == '8X', ],
  aes(x = COVERAGE_CATEGORY, y = N, fill = OLD_AGREE_2)) +
  geom_bar(stat = 'identity') +
  ggtitle(paste('8X samples by category\nand agreement',
    ' with old metadata', sep = '')) +
  guides(fill=guide_legend(title='Expected?')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

gg_bar_new <- ggplot(meta_tab_2[COVERAGE_CATEGORY == 'new_164g', ], 
  aes(x = NQUIRE_CDS_20, y = N, fill = OLD_AGREE_2)) +
  geom_bar(stat = 'identity') +
  ggtitle(paste('nQuire-based ploidy in 164 new samples\nand agreement', 
    ' with old metadata', sep = '')) +
  guides(fill=guide_legend(title='Expected?'))

ploidy_bar_1 <- paste('/global/cscratch1/sd/grabowsp/sg_8X_scratch/',
  'nquire_scratch/', 'ploidy_by_category_bar_1.pdf', sep = '')

pdf(ploidy_bar_1, width = 8, height = 9)
(gg_bar_1 / gg_bar_new) | (gg_bar_4X /gg_bar_6X /gg_bar_8X)
dev.off()


quit(save = 'no')


