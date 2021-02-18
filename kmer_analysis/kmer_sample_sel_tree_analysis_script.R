# Script for generating sample file for Adam to analyze 4X/8X with Kmers
# Run locally because had trouble getting all packages loaded on Cori

### LOAD PACKAGES ###
# install.packages('tidyverse')
library(tidyverse)
# install.packages('BiocManager')
# BiocManager::install('ggtree')
library(ggtree)
library('tidytree')
#install.packages('data.table')
library('data.table')

### INPUT DATA ###
# Trees used for selecting samples
# full NJ tree
nj_tree_file <- '/Users/grabowsk/data/sg_8X_analysis/allsamp_NJ_tree.newick'
nj_2 <- read.tree(nj_tree_file)
nj_3 <- as_tibble(nj_2)
as.treedata(nj_3)

# Atlantic NJ tree
atl_tree_file <- '/Users/grabowsk/data/sg_8X_analysis/atlantic_nj_tree.newick'
atl_2 <- read.tree(atl_tree_file)
atl_3 <- as_tibble(atl_2)
as.treedata(atl_3)

# MW NJ tree
mw_tree_file <- '/Users/grabowsk/data/sg_8X_analysis/midwest_nj_tree.newick'
mw_2 <- read.tree(mw_tree_file)
mw_3 <- as_tibble(mw_2)
as.treedata(mw_3)

# Gulf NJ tree
gulf_tree_file <- '/Users/grabowsk/data/sg_8X_analysis/gulf_nj_tree.newick'
gulf_2 <- read.tree(gulf_tree_file)
gulf_3 <- as_tibble(gulf_2)
as.treedata(gulf_3)

# sample files
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

### SET OUTPUT ###

final_tab_out <- '/Users/grabowsk/data/sg_8X_analysis/kmer_4x8x_libs_v1.csv'

### MAKE FUNCTIONS ###

get_filt_names <- function(samp_names, tot_filt_set){
  # Function for pulling out only the "samp_names" that are part of
  #  tot_filt_set; also remove NA's
  tmp_names <- intersect(tot_filt_set, samp_names[-which(is.na(samp_names))])
  return(tmp_names)
}

##########
 
### Group 1 Sample Sets ###
grp_1_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp1', ID],
                        natv2_samps)
# Atlantic ancestry; most have > 10% Gulf ancestry, but not all...
# 13/13 are in filtered sample set
# J192.L1 is separate from other grp_1 samps and clusters with other samples
#  from same ACC
# J593.A is probably 6X and is separate from other grp_1 samps and clusters
#   with other sample from same ACC
# J180.L2 clusters within its own ACC

## Start total table ##
tot_samp_tab <- data.table(SAMP_GROUP = 'grp_1_8X_all', 
                           SAMP_NAMES = grp_1_filt)

## Group 1 8X subsets ##
grp_1_1_8X <- setdiff(grp_1_filt, c('J192.L1', 'J593.A', 'J180.L2'))
# FL 8X with Atlantic background
# 55-85% Atlantic, rest mainly Gulf
tot_samp_tab <- rbind(tot_samp_tab, 
                      data.table(SAMP_GROUP = 'grp_1_1_8X', 
                                 SAMP_NAMES = grp_1_1_8X))
#
grp_1_2_8X <- 'J192.L1'
# NC Pure Atlantic, clusters with same ACC
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_2_8X',
                                 SAMP_NAMES = grp_1_2_8X))
#
grp_1_3_8X <- 'J593.A'
# NC almost pure Atlantic (0.4% Gulf), clusters with same ACC
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_3_8X',
                                 SAMP_NAMES = grp_1_3_8X))
#
grp_1_4_8X <- 'J180.L2'
# FL Altantic with 2% Gulf
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_4_8X',
                                 SAMP_NAMES = grp_1_4_8X))
########
# Choose 4X comparison samples for 8X test samples
## Comparison samples for full grp_1 ##
grp_1_full_4X_tmp <- unique(
  offspring(atl_3, MRCA(atl_3, 'J593.B', 'J680.A')$node)$label)
grp_1_all_4X <- sample(setdiff(get_filt_names(grp_1_full_4X_tmp, natv2_samps),
                         grp_1_filt), size = 15)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_all_4X',
                                 SAMP_NAMES = grp_1_all_4X))
## Group 1_1 comparison samples ##
grp_1_1_4X_1 <- setdiff(intersect(natv2_samps,
                          c(samp_meta[ACC == 'J185', VCF_NAME],
                            samp_meta[ACC == 'J186', VCF_NAME])),
                        grp_1_filt)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_1_4X_1',
                                 SAMP_NAMES = grp_1_1_4X_1))
# Group 1_1 Atlantic samples
grp_1_1_4X_t2 <- unique(
  offspring(atl_3, MRCA(atl_3, 'J577.B', 'J019.L1')$node)$label) 
grp_1_1_4X_t2 <- setdiff(get_filt_names(grp_1_1_4X_t2, natv2_samps),
                         grp_1_filt)
grp_1_1_4X_2 <- c(grp_1_1_4X_1, sample(grp_1_1_4X_t2, size = 10))
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_1_4X_2',
                                 SAMP_NAMES = grp_1_1_4X_2))
# Group 1_1 Gulf samples
grp_1_1_4X_tmp <- unique(
  offspring(nj_3, MRCA(nj_3, 'J499.C', 'J456.A')$node)$label)
grp_1_1_4X_tmp <- setdiff(get_filt_names(grp_1_1_4X_tmp, natv2_samps), 
                          grp_1_filt)
grp_1_1_4X_3 <- c(grp_1_1_4X_1, sample(grp_1_1_4X_tmp, size = 10))
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_1_4X_3',
                                 SAMP_NAMES = grp_1_1_4X_3))
# Group 1_1 Atlantic and Gulf samples
grp_1_1_4X_4 <- c(grp_1_1_4X_1, sample(grp_1_1_4X_t2, size = 6), 
                  sample(grp_1_1_4X_tmp, size = 6))
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_1_4X_4',
                                 SAMP_NAMES = grp_1_1_4X_4))
## Group 1_2 comparison samples ##
grp_1_2_4X_1 <- setdiff(intersect(natv2_samps,
                                  c(samp_meta[ACC == 'J192', VCF_NAME],
                                    samp_meta[ACC == 'J667', VCF_NAME],
                                    samp_meta[ACC == 'J193', VCF_NAME])),
                        grp_1_filt)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_2_4X_1',
                                 SAMP_NAMES = grp_1_2_4X_1))
grp_1_2_4X_tmp <- unique(
  offspring(nj_3, MRCA(nj_3, 'J193.A', 'J576.C')$node)$label)
grp_1_2_4X_tmp <- setdiff(get_filt_names(grp_1_2_4X_tmp, natv2_samps), 
                          grp_1_filt)
grp_1_2_4X_2 <- sample(grp_1_2_4X_tmp, size = 10)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_2_4X_2',
                                 SAMP_NAMES = grp_1_2_4X_2))
## Group 1_3 comparions samples ##
grp_1_3_4X_1 <- setdiff(intersect(natv2_samps,
                                  c(samp_meta[ACC == 'J593', VCF_NAME],
                                    samp_meta[ACC == 'J677', VCF_NAME],
                                    samp_meta[ACC == 'J595', VCF_NAME])),
                        grp_1_filt)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_3_4X_1',
                                 SAMP_NAMES = grp_1_3_4X_1))
grp_1_3_4X_tmp <- unique(
  offspring(nj_3, MRCA(nj_3, 'J595.B', 'J204.A')$node)$label)
grp_1_3_4X_tmp <- setdiff(get_filt_names(grp_1_3_4X_tmp, natv2_samps), 
                          grp_1_filt)
grp_1_3_4X_2 <- sample(grp_1_3_4X_tmp, size = 10)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_3_4X_2',
                                 SAMP_NAMES = grp_1_3_4X_2))
## Group 1_4 comparison samples ##
grp_1_4_4X_1 <- setdiff(intersect(natv2_samps,
                                  c(samp_meta[ACC == 'J180', VCF_NAME],
                                    samp_meta[ACC == 'J181', VCF_NAME],
                                    samp_meta[ACC == 'J182', VCF_NAME],
                                    samp_meta[ACC == 'J019', VCF_NAME],
                                    samp_meta[ACC == 'J164', VCF_NAME])),
                        grp_1_filt)
# will only use one group for 1_3 because next closest cluster includes an 8X
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_4_4X_1',
                                 SAMP_NAMES = grp_1_4_4X_1))
## Group 1 control - Atlantic Cluster ##
atlantic_test <- offspring(nj_3, MRCA(nj_3, 'J527.C', 'J644.C')$node)$label
# intersect(atlantic_test, grp_1_filt)
# character(0) - can use this as a test group
atlantic_test <- get_filt_names(atlantic_test, natv2_samps)
grp_1_4X_ATL_control <- sample(atlantic_test, size = 15)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_4X_ATL_control',
                                 SAMP_NAMES = grp_1_4X_ATL_control))
## Group 1 control - MW cluster
mw_test <- offspring(nj_3, MRCA(nj_3, 'J581.C', 'J368.A')$node)$label
mw_test <- get_filt_names(mw_test, natv2_samps)
grp_1_4X_MW_control <- sample(mw_test, size = 15)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_1_4X_MW_control',
                                 SAMP_NAMES = grp_1_4X_MW_control))

#########

### Group 2 Sample Sets ###
grp_2_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp2', ID],
                        natv2_samps)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_2_8X_all',
                                 SAMP_NAMES = grp_2_filt))

## Grp 2 comparison samples
grp_2_4X_test <- offspring(mw_3, MRCA(mw_3, 'J471_B', 'J658.B')$node)$label
grp_2_4X_1 <- get_filt_names(grp_2_4X_test, natv2_samps)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_2_4X',
                                 SAMP_NAMES = grp_2_4X_1))

## Grp 2 MW control
grp_2_4X_MW_tmp <- offspring(mw_3, MRCA(mw_3, 'J582.C', 'J420.B')$node)$label
grp_2_4X_MW_tmp <- get_filt_names(grp_2_4X_MW_tmp, natv2_samps)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_2_4X_MW_control',
                                 SAMP_NAMES = grp_2_4X_MW_tmp))
## Grp 2 Atlantic control
grp_2_4X_ATL_tmp <- atlantic_test
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_2_4X_ATL_control',
                                 SAMP_NAMES = grp_2_4X_ATL_tmp))

### Group 4 Sample Sets ###
grp_4_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp4', ID],
                        natv2_samps)
grp_4_1_8X <- setdiff(grp_4_filt, 'J400.A')
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_4_1_8X',
                                 SAMP_NAMES = grp_4_1_8X))
grp_4_2_8X <- 'J400.A'
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_4_2_8X',
                                 SAMP_NAMES = grp_4_2_8X))
## Grp 4 comparison sets ##
grp_4_1_4X <- grp_2_4X_1
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_4_1_4X',
                                 SAMP_NAMES = grp_4_1_4X))
grp_4_2_4X_tmp <- offspring(mw_3, MRCA(mw_3, 'J400.B', 'J359.A')$node)$label
grp_4_2_4X_tmp <- setdiff(get_filt_names(grp_4_2_4X_tmp, natv2_samps), 
                          grp_4_2_8X)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_4_2_4X',
                                 SAMP_NAMES = grp_4_2_4X_tmp))

grp_4_4X_MW_tmp <- grp_2_4X_MW_tmp
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_4_4X_MW_control',
                                 SAMP_NAMES = grp_4_4X_MW_tmp))
grp_4_4X_ATL_tmp <- grp_2_4X_ATL_tmp
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_4_4X_ATL_control',
                                 SAMP_NAMES = grp_4_4X_ATL_tmp))

### Group 5 sample sets ###
grp_5_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp5', ID],
                        natv2_samps)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_8X_all',
                                 SAMP_NAMES = grp_5_filt))
grp_5_1_8X <- grp_5_filt[
  grep('J021|J236|J304|J302|J221|J213|J328|J309', grp_5_filt)]
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_1_8X',
                                 SAMP_NAMES = grp_5_1_8X))
grp_5_2_8X <- grp_5_filt[grep('J464|J462|J495', grp_5_filt)]
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_2_8X',
                                 SAMP_NAMES = grp_5_2_8X))
grp_5_3_8X <- 'J188.L1'
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_3_8X',
                                 SAMP_NAMES = grp_5_3_8X))
grp_5_4_8X <- 'J249.A'
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_4_8X',
                                 SAMP_NAMES = grp_5_4_8X))

## Comparison groups for Grp 5 ##
grp_5_1_4X_tmp <- offspring(gulf_3, 
                            MRCA(gulf_3, 'J309_A', 'J021.B')$node)$label
grp_5_1_4X_tmp <- setdiff(get_filt_names(grp_5_1_4X_tmp, natv2_samps), 
                          grp_5_filt)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_1_4X',
                                 SAMP_NAMES = grp_5_1_4X_tmp))

grp_5_2_4X_tmp <- offspring(gulf_3, 
                            MRCA(gulf_3, 'J482.A', 'J465.B')$node)$label
grp_5_2_4X_tmp <- setdiff(get_filt_names(grp_5_2_4X_tmp, natv2_samps),
                          grp_5_filt)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_2_4X',
                                 SAMP_NAMES = grp_5_2_4X_tmp))

grp_5_3_4X_tmp <- offspring(gulf_3, 
                            MRCA(gulf_3, 'J191.A', 'J190.A')$node)$label
grp_5_3_4X_tmp <- setdiff(get_filt_names(grp_5_3_4X_tmp, natv2_samps),
                          grp_5_filt)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_3_4X',
                                 SAMP_NAMES = grp_5_3_4X_tmp))

grp_5_4_4X_tmp <- offspring(gulf_3, 
                            MRCA(gulf_3, 'J237.A', 'J272.A')$node)$label
grp_5_4_4X_tmp <- setdiff(get_filt_names(grp_5_4_4X_tmp, natv2_samps),
                          grp_5_filt)
grp_5_4_4X <- sample(grp_5_4_4X_tmp, size = 10)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_4_4X',
                                 SAMP_NAMES = grp_5_4_4X))

grp_5_gulf_tmp <- offspring(gulf_3, 
                            MRCA(gulf_3, 'J323.A', 'J335.A')$node)$label
grp_5_gulf_tmp <- get_filt_names(grp_5_gulf_tmp, natv2_samps)
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_4X_GULF_control',
                                 SAMP_NAMES = grp_5_gulf_tmp))

grp_5_MW_control <- grp_1_4X_MW_control
tot_samp_tab <- rbind(tot_samp_tab,
                      data.table(SAMP_GROUP = 'grp_5_4X_MW_control',
                                 SAMP_NAMES = grp_5_MW_control))

### Assign library names to sample names
tot_samp_tab[, LIB := as.character(NA)]
for(i in seq(nrow(tot_samp_tab))){
  tmp_meta_ind <- which(samp_meta$VCF_NAME == tot_samp_tab$SAMP_NAMES[i])
  tot_samp_tab[i, LIB := samp_meta$LIB[tmp_meta_ind]]
}

fwrite(tot_samp_tab, final_tab_out)
