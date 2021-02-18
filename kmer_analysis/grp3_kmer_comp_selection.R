# Script to use geographic distance to select samples for comparing to Grp3
#   samples

### LOAD PACKAGES ###
library(tidyverse)
library('data.table')
#install.packages('geodist')
library(geodist)

### INPUT DATA ####
natv2_name_file <- '/Users/grabowsk/Analysis/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_name_file, header = F)$V1

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

sample_dir <- '/Users/grabowsk/data/sg_8X_analysis/'

atlantic_file <- paste(sample_dir, 'Atlantic_ref_samps_40.txt', sep = '')
atl_samps <- fread(atlantic_file, header = T)[, ID]

gulf_file <- paste(sample_dir, 'Gulf_ref_samps_40.txt', sep = '')
gulf_samps <- fread(gulf_file, header = T)[, ID]

mw_file <- paste(sample_dir, 'Midwest_ref_samps_40.txt', sep = '')
mw_samps <- fread(mw_file, header = T)[, ID]

candidate_samp_file <- paste(sample_dir, 
                             'Introgression_target_individuals.txt', sep = '')
candidate_samps <- fread(candidate_samp_file)

grp3_file <- paste(sample_dir, 'Grp4_subgroups.csv', sep = '')
grp3_samps <- fread(grp3_file)

admix_res_file <- '/Users/grabowsk/data/sg_8X_analysis/GW_50k_geobig.3.Q'
admix_res_0 <- fread(admix_res_file)
# 1 = Atlantic, 2 = Gulf, 3 = Midwest

admix_name_file <- '/Users/grabowsk/data/sg_8X_analysis/GW_50k_geobig.fam'
admix_names <- fread(admix_name_file)$V1

admix_res <- cbind(admix_names, admix_res_0)
colnames(admix_res) <- c('SAMP_NAME', 'ATL', 'GULF', 'MW')

### SET OUTPUTS ###

final_grp3_tab_out <- '/Users/grabowsk/data/sg_8X_analysis/kmer_grp3_comp_libs_v1.csv'

### SET VARIABLES ###
# the number of samples for each test group
n_test_samps <- 20

##############

# Select Groups

grp_1_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp1', ID],
                        natv2_samps)
grp_2_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp2', ID],
                        natv2_samps)
grp_4_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp4', ID],
                        natv2_samps)
grp_5_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp5', ID],
                        natv2_samps)
grp_3_1_filt <- intersect(grp3_samps[Sub_Group == 'subgrp1', ID],
                          natv2_samps)
grp_3_2_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp2', ID],
                           natv2_samps)
grp_3_3_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp3', ID],
                           natv2_samps)

grp_3_4_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp4', ID],
                           natv2_samps)

# Get gene pool names
atl_4X_samps <- intersect(natv2_samps, intersect(
  samp_meta[NQUIRE_PLOIDY == '4X', VCF_NAME], admix_res[ATL > 0.9, SAMP_NAME]
))
atl_all_samps <- intersect(natv2_samps, admix_res[ATL > 0.9, SAMP_NAME])

mw_4X_samps <- intersect(natv2_samps, intersect(
  samp_meta[NQUIRE_PLOIDY == '4X', VCF_NAME], admix_res[MW > 0.9, SAMP_NAME]
))
mw_all_samps <- intersect(natv2_samps, admix_res[MW > 0.9, SAMP_NAME])

gulf_4X_samps <- intersect(natv2_samps, intersect(
  samp_meta[NQUIRE_PLOIDY == '4X', VCF_NAME], admix_res[GULF > 0.9, SAMP_NAME]
))
gulf_all_samps <- intersect(natv2_samps, admix_res[GULF > 0.9, SAMP_NAME])

###

meta_nat_inds <- which(samp_meta$VCF_NAME %in% natv2_samps)

geo_tab <- samp_meta[meta_nat_inds, c('LONGITUDE', 'LATITUDE')]

geo_dist_mat <- geodist(geo_tab, measure = 'geodesic')
rownames(geo_dist_mat) <- colnames(geo_dist_mat) <- samp_meta$VCF_NAME[
  meta_nat_inds]

grp3_2_inds <- which(rownames(geo_dist_mat) %in% grp_3_2_filt)
grp3_3_inds <- which(rownames(geo_dist_mat) %in% grp_3_3_filt)

atl_4_inds <- which(rownames(geo_dist_mat) %in% atl_4X_samps)
mw_4_inds <- which(rownames(geo_dist_mat) %in% mw_4X_samps)
gulf_4_inds <- which(rownames(geo_dist_mat) %in% gulf_4X_samps)

atl_all_inds <- which(rownames(geo_dist_mat) %in% atl_all_samps)
mw_all_inds <- which(rownames(geo_dist_mat) %in% mw_all_samps)
gulf_all_inds <- which(rownames(geo_dist_mat) %in% gulf_all_samps)

# start sample table
grp3_comp_tab <- data.table(SAMP_GROUP = 'grp3_2', 
                            SAMP_NAMES = grp_3_2_filt)

### Select closest 4X samples to grp3_2 for each gene pool

grp3_2_atl4_dist <- geo_dist_mat[grp3_2_inds, atl_4_inds]
grp3_2_atl4_dist_2 <- apply(grp3_2_atl4_dist, 2, mean)
# sort(grp3_2_atl4_dist_2, decreasing = F)
grp3_2_atl4_names <- names(sort(grp3_2_atl4_dist_2, 
                                decreasing = F))[seq(n_test_samps)]
grp3_2_atl4_control <- names(sort(grp3_2_atl4_dist_2, 
                                  decreasing = T))[seq(n_test_samps)]

grp3_2_mw4_dist <- geo_dist_mat[grp3_2_inds, mw_4_inds]
grp3_2_mw4_dist_2 <- apply(grp3_2_mw4_dist, 2, mean)
# sort(grp3_2_mw4_dist_2, decreasing = F)
grp3_2_mw4_names <- names(sort(grp3_2_mw4_dist_2, 
                               decreasing = F))[seq(n_test_samps)]
grp3_2_mw4_control <- names(sort(grp3_2_mw4_dist_2, 
                               decreasing = T))[seq(n_test_samps)]

grp3_2_gulf4_dist <- geo_dist_mat[grp3_2_inds, gulf_4_inds]
grp3_2_gulf4_dist_2 <- apply(grp3_2_gulf4_dist, 2, mean)
grp3_2_gulf4_names <- names(sort(grp3_2_gulf4_dist_2, 
                                 decreasing = F))[seq(n_test_samps)]
grp3_2_gulf4_control <- names(sort(grp3_2_gulf4_dist_2, 
                                 decreasing = T))[seq(n_test_samps)]

grp3_comp_tab <- rbindlist(list(grp3_comp_tab, 
                           data.table(SAMP_GROUP = 'grp3_2_ATL_4X_NEAR',
                                      SAMP_NAMES = grp3_2_atl4_names),
                           data.table(SAMP_GROUP = 'grp3_2_MW_4X_NEAR',
                                      SAMP_NAMES = grp3_2_mw4_names),
                           data.table(SAMP_GROUP = 'grp3_2_GULF_4X_NEAR',
                                      SAMP_NAMES = grp3_2_gulf4_names),
                           data.table(SAMP_GROUP = 'grp3_2_ATL_4X_FAR',
                                      SAMP_NAMES = grp3_2_atl4_control),
                           data.table(SAMP_GROUP = 'grp3_2_MW_4X_FAR',
                                      SAMP_NAMES = grp3_2_mw4_control),
                           data.table(SAMP_GROUP = 'grp3_2_GULF_4X_FAR',
                                      SAMP_NAMES = grp3_2_gulf4_control)
                           ))

### Select closest samples of any ploidy to grp3_2 for each gene pool
grp3_2_atlall_dist <- geo_dist_mat[grp3_2_inds, atl_all_inds]
grp3_2_atlall_dist_2 <- apply(grp3_2_atlall_dist, 2, mean)
# sort(grp3_2_atl4_dist_2, decreasing = F)
grp3_2_atlall_names <- names(sort(grp3_2_atlall_dist_2, 
                                decreasing = F))[seq(n_test_samps)]
# same as 4X-only samps
grp3_2_atlall_control <- names(sort(grp3_2_atlall_dist_2, 
                                  decreasing = T))[seq(n_test_samps)]
# same as 4X-only samps

grp3_2_mwall_dist <- geo_dist_mat[grp3_2_inds, mw_all_inds]
grp3_2_mwall_dist_2 <- apply(grp3_2_mwall_dist, 2, mean)
# sort(grp3_2_atl4_dist_2, decreasing = F)
grp3_2_mwall_names <- names(sort(grp3_2_mwall_dist_2, 
                                  decreasing = F))[seq(n_test_samps)]
# 6 differences from 4X-only samps
grp3_2_mwall_control <- names(sort(grp3_2_mwall_dist_2, 
                                 decreasing = T))[seq(n_test_samps)]
# completely different from 4X-only samps

grp3_2_gulfall_dist <- geo_dist_mat[grp3_2_inds, gulf_all_inds]
grp3_2_gulfall_dist_2 <- apply(grp3_2_gulfall_dist, 2, mean)
# sort(grp3_2_atl4_dist_2, decreasing = F)
grp3_2_gulfall_names <- names(sort(grp3_2_gulfall_dist_2, 
                                 decreasing = F))[seq(n_test_samps)]
# same as 4X-only samps
grp3_2_gulfall_control <- names(sort(grp3_2_gulfall_dist_2, 
                                   decreasing = T))[seq(n_test_samps)]
# 1 difference from 4X-only samps

grp3_comp_tab <- rbindlist(list(grp3_comp_tab, 
                                data.table(SAMP_GROUP = 'grp3_2_ATL_NEAR',
                                           SAMP_NAMES = grp3_2_atlall_names),
                                data.table(SAMP_GROUP = 'grp3_2_ATL_FAR',
                                           SAMP_NAMES = grp3_2_atlall_control),
                                data.table(SAMP_GROUP = 'grp3_2_MW_NEAR',
                                           SAMP_NAMES = grp3_2_mwall_names),
                                data.table(SAMP_GROUP = 'grp3_2_MW_FAR',
                                           SAMP_NAMES = grp3_2_mwall_control),
                                data.table(SAMP_GROUP = 'grp3_2_GULF_NEAR',
                                           SAMP_NAMES = grp3_2_gulfall_names),
                                data.table(SAMP_GROUP = 'grp3_2_GULF_FAR',
                                           SAMP_NAMES = grp3_2_gulfall_control)
                                ))

### Group3_3

grp3_comp_tab <- rbind(grp3_comp_tab,
                       data.table(SAMP_GROUP = 'grp3_3',
                                  SAMP_NAMES = grp_3_3_filt))

grp3_3_atl4_dist <- geo_dist_mat[grp3_3_inds, atl_4_inds]
grp3_3_atl4_dist_2 <- apply(grp3_3_atl4_dist, 2, mean)
grp3_3_atl4_names <- names(sort(grp3_3_atl4_dist_2, 
                                decreasing = F))[seq(n_test_samps)]
grp3_3_atl4_control <- names(sort(grp3_3_atl4_dist_2, 
                                  decreasing = T))[seq(n_test_samps)]

grp3_3_mw4_dist <- geo_dist_mat[grp3_3_inds, mw_4_inds]
grp3_3_mw4_dist_2 <- apply(grp3_3_mw4_dist, 2, mean)
grp3_3_mw4_names <- names(sort(grp3_3_mw4_dist_2, 
                               decreasing = F))[seq(n_test_samps)]
grp3_3_mw4_control <- names(sort(grp3_3_mw4_dist_2, 
                                 decreasing = T))[seq(n_test_samps)]

grp3_3_gulf4_dist <- geo_dist_mat[grp3_3_inds, gulf_4_inds]
grp3_3_gulf4_dist_2 <- apply(grp3_3_gulf4_dist, 2, mean)
grp3_3_gulf4_names <- names(sort(grp3_3_gulf4_dist_2, 
                                 decreasing = F))[seq(n_test_samps)]
grp3_3_gulf4_control <- names(sort(grp3_3_gulf4_dist_2, 
                                   decreasing = T))[seq(n_test_samps)]

grp3_comp_tab <- rbindlist(list(grp3_comp_tab, 
                                data.table(SAMP_GROUP = 'grp3_3_ATL_4X_NEAR',
                                           SAMP_NAMES = grp3_3_atl4_names),
                                data.table(SAMP_GROUP = 'grp3_3_MW_4X_NEAR',
                                           SAMP_NAMES = grp3_3_mw4_names),
                                data.table(SAMP_GROUP = 'grp3_3_GULF_4X_NEAR',
                                           SAMP_NAMES = grp3_3_gulf4_names),
                                data.table(SAMP_GROUP = 'grp3_3_ATL_4X_FAR',
                                           SAMP_NAMES = grp3_3_atl4_control),
                                data.table(SAMP_GROUP = 'grp3_3_MW_4X_FAR',
                                           SAMP_NAMES = grp3_3_mw4_control),
                                data.table(SAMP_GROUP = 'grp3_3_GULF_4X_FAR',
                                           SAMP_NAMES = grp3_3_gulf4_control)
))

### Select samples of any ploidy to grp3_3 for each gene pool
grp3_3_atlall_dist <- geo_dist_mat[grp3_3_inds, atl_all_inds]
grp3_3_atlall_dist_2 <- apply(grp3_3_atlall_dist, 2, mean)
grp3_3_atlall_names <- names(sort(grp3_3_atlall_dist_2, 
                                  decreasing = F))[seq(n_test_samps)]
# same as 4X-only samps
grp3_3_atlall_control <- names(sort(grp3_3_atlall_dist_2, 
                                    decreasing = T))[seq(n_test_samps)]
# same as 4X-only samps

grp3_3_mwall_dist <- geo_dist_mat[grp3_3_inds, mw_all_inds]
grp3_3_mwall_dist_2 <- apply(grp3_3_mwall_dist, 2, mean)
grp3_3_mwall_names <- names(sort(grp3_3_mwall_dist_2, 
                                 decreasing = F))[seq(n_test_samps)]
# 8 differences from 4X-only samps
grp3_3_mwall_control <- names(sort(grp3_3_mwall_dist_2, 
                                   decreasing = T))[seq(n_test_samps)]
# 9 differences from 4X-only samps

grp3_3_gulfall_dist <- geo_dist_mat[grp3_3_inds, gulf_all_inds]
grp3_3_gulfall_dist_2 <- apply(grp3_3_gulfall_dist, 2, mean)
grp3_3_gulfall_names <- names(sort(grp3_3_gulfall_dist_2, 
                                   decreasing = F))[seq(n_test_samps)]
# same as 4X-only samps
grp3_3_gulfall_control <- names(sort(grp3_3_gulfall_dist_2, 
                                     decreasing = T))[seq(n_test_samps)]
# same as 4X-only samps
grp3_comp_tab <- rbindlist(list(grp3_comp_tab, 
                                data.table(SAMP_GROUP = 'grp3_3_ATL_NEAR',
                                           SAMP_NAMES = grp3_3_atlall_names),
                                data.table(SAMP_GROUP = 'grp3_3_ATL_FAR',
                                           SAMP_NAMES = grp3_3_atlall_control),
                                data.table(SAMP_GROUP = 'grp3_3_MW_NEAR',
                                           SAMP_NAMES = grp3_3_mwall_names),
                                data.table(SAMP_GROUP = 'grp3_3_MW_FAR',
                                           SAMP_NAMES = grp3_3_mwall_control),
                                data.table(SAMP_GROUP = 'grp3_3_GULF_NEAR',
                                           SAMP_NAMES = grp3_3_gulfall_names),
                                data.table(SAMP_GROUP = 'grp3_3_GULF_FAR',
                                           SAMP_NAMES = grp3_3_gulfall_control)
))

# add LIBRARY ID
grp3_comp_tab[, LIB := as.character(NA)]
for(i in seq(nrow(grp3_comp_tab))){
  tmp_meta_ind <- which(samp_meta$VCF_NAME == grp3_comp_tab$SAMP_NAMES[i])
  grp3_comp_tab[i, LIB := samp_meta$LIB[tmp_meta_ind]]
}

fwrite(grp3_comp_tab, file = final_grp3_tab_out)


#quit(save = 'no')