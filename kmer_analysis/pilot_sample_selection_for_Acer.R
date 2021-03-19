# Choose pilot set of 8X and 4X samples for Acer to use for kmer analysis
# Want 30 samples total: 15 4X and 15 8X
## for 8X: 3 Grp1, 3 Grp2, 3 Grp3, 3 Grp4, 3 Grp5
## For 4X: 3 pure ATL, 3 pure Gulf, 3 pure MW, 2 ATL/MW, 2 ATL/Gulf, 2 Gulf/MW

### LOAD PACKAGES ###
library(tidyverse)
library('data.table')

### INPUT DATA ###
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
acer_tab_out <- '/Users/grabowsk/data/sg_8X_analysis/kmer_acer_pilot_samps.csv'

### SET VARIABLES ###


#########

# samples from Introgression analysis
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


# grp1 samples

grp1_sub <- c('J072_A', 'J593.A', 'J192.L1')
atl_4X <- c('J640_A', 'J576.B', 'J671.A')
atl_gulf <- c('J186.A', 'J613.B')
# atl_mw <- c(), not great options
grp2_sub <- c('J346.A', 'J417_A', 'J560.C')
grp4_sub <- c('J474.A', 'J485.B', 'J214_A')
mw_4X <- c('J349.B', 'J470.A', 'J412.A')
# mw_atl <- c(), not great options
# mw_gulf <- c(), not great 4X options
grp3_sub <- c('J468_A', 'J579.B', 'J200_A')
grp5_sub <- c('J221_A', 'J462.A', 'J188.L1')
gulf_4X <- c('J314.A', 'J336.A', 'J273.A', 'J303.A')
gulf_mw <- c('J477.A')
gulf_atl <- c('J501.C', 'J463.A')

acer_pilot_tab <- rbindlist(list(data.table(SAMP_GROUP = 'grp1_ATL_8X',
                                            SAMP_NAMES = grp1_sub),
                                 data.table(SAMP_GROUP = 'ATL_4X',
                                            SAMP_NAMES = atl_4X),
                                 data.table(SAMP_GROUP = 'ATL_GULF_4X',
                                            SAMP_NAMES = atl_gulf),
                                 data.table(SAMP_GROUP = 'grp2_MW_8X',
                                            SAMP_NAMES = grp2_sub),
                                 data.table(SAMP_GROUP = 'grp4_MW_8X',
                                            SAMP_NAMES = grp4_sub),
                                 data.table(SAMP_GROUP = 'MW_4X',
                                            SAMP_NAMES = mw_4X),
                                 data.table(SAMP_GROUP = 'grp3_MIX_8X',
                                            SAMP_NAMES = grp3_sub),
                                 data.table(SAMP_GROUP = 'grp5_GULF_8X',
                                            SAMP_NAMES = grp5_sub),
                                 data.table(SAMP_GROUP = 'GULF_4X',
                                            SAMP_NAMES = gulf_4X),
                                 data.table(SAMP_GROUP = 'GULF_MW_4X',
                                            SAMP_NAMES = gulf_mw),
                                 data.table(SAMP_GROUP = 'GULF_ATL_4X',
                                            SAMP_NAMES = gulf_atl)
                                 ))

acer_pilot_tab[, PLOIDY := as.character(NA)]
acer_pilot_tab[grep('8X', acer_pilot_tab$SAMP_GROUP), PLOIDY := '8X']
acer_pilot_tab[grep('4X', acer_pilot_tab$SAMP_GROUP), PLOIDY := '4X']

acer_pilot_tab[, LIB := as.character(NA)]
for(i in seq(nrow(acer_pilot_tab))){
  tmp_m_ind <- which(samp_meta$VCF_NAME == acer_pilot_tab$SAMP_NAMES[i])
  acer_pilot_tab[i, LIB := samp_meta$LIB[tmp_m_ind]]
}

fwrite(acer_pilot_tab, file = acer_tab_out)

# quit(save = 'no')