# Selecting updated training sets for the first draft of the manuscript
#
## run on HA

### SET ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
# population structure results file
res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_file)

# metadata file
meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

### SET OUTPUTS ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/'

mw_train_full_file <- paste(out_dir, 
  'MW_firstdraft_train_all60.txt', sep = '')
mw_train_filt_file <- paste(out_dir, 
  'MW_firstdraft_train_filt40.txt', sep = '')

atl_train_full_file <- paste(out_dir, 
  'ATL_firstdraft_train_all60.txt', sep = '')
atl_train_filt_file <- paste(out_dir, 
  'ATL_firstdraft_train_filt40.txt', sep = '')

gulf_train_full_file <- paste(out_dir,
  'GULF_firstdraft_train_all60.txt', sep = '')
gulf_train_filt_file <- paste(out_dir,
  'GULF_firstdraft_train_filt40.txt', sep = '')

### SET VARIABLES ###
# how many total samples to choose from each genepool
n_gp <- 60

# how many filtered samples to choose from each genepool
n_filt <- 40

# maximum ATL latitude to minimize presence of MW introgressions into 
#  northern ATL in the training set
max_atl_lat <- 40.2

#######################

# MW training samples
top_pc1 <- order(res_tab$full_pc01_stand, decreasing = T)

top_pc1_inds <- top_pc1[1:n_gp]

# make sure is all 4X
sum(res_tab[top_pc1_inds, ploidy] != '4X') == 0
#[1] TRUE

# make sure is all 99+% MW
sum(res_tab[top_pc1_inds, full_admix_k3_MW] < 0.99) == 0
#[1] TRUE

# how many subgrps
length(unique(res_tab[top_pc1_inds, subgrp_v2]))
# 4
res_tab[top_pc1_inds, .N, by = subgrp_v2]
# is mainly MW_10 (40 of 60)

tmp_mw_names <- res_tab[top_pc1_inds, samp_name]

# how many accessions
length(unique(samp_meta[VCF_NAME %in% tmp_mw_names, UNI_ACC]))
# 34 accessions

# how many states
length(unique(samp_meta[VCF_NAME %in% tmp_mw_names, STATE]))
# 8 states

# these 60 seem good
tot_mw_names <- tmp_mw_names

# Select 40 across accessions
sub_mw_names <- c()
for(acc_name in unique(samp_meta[VCF_NAME %in% tot_mw_names, UNI_ACC])){
  tmp_names <- intersect(tot_mw_names, samp_meta[UNI_ACC == acc_name, VCF_NAME])
  tmp_choice <- sample(tmp_names, size = 1)
  sub_mw_names <- c(sub_mw_names, tmp_choice) 
  }

mw_leftover <- setdiff(tot_mw_names, sub_mw_names)

length(unique(samp_meta[VCF_NAME %in% mw_leftover, UNI_ACC]))
# 15
length(unique(samp_meta[VCF_NAME %in% mw_leftover, STATE]))
# 7

for(stn in sample(unique(samp_meta[VCF_NAME %in% mw_leftover, STATE]), 
  size = 6)){
  tmp_names <- intersect(mw_leftover, samp_meta[STATE == stn, VCF_NAME])
  tmp_choice <- sample(tmp_names, size = 1)
  sub_mw_names <- c(sub_mw_names, tmp_choice)
}

length(sub_mw_names)
# 40

fwrite(data.table(samp = tot_mw_names), file = mw_train_full_file,
  sep = '\t', row.names = F, col.names = F)
fwrite(data.table(samp = sub_mw_names), file = mw_train_filt_file,
  sep = '\t', row.names = F, col.names = F)

#########

# Atlantic training samples

north_names <- union(samp_meta[LATITUDE >= 40.2, VCF_NAME], 
  samp_meta[is.na(STATE), VCF_NAME])
north_inds <- which(res_tab$samp_name %in% north_names)

bot_pc1 <- order(res_tab$full_pc01_stand, decreasing = F)
atl_1_samp_pos <- rep(NA, times = length(bot_pc1))
atl_1_samp_pos[bot_pc1] <- seq(length(bot_pc1))

top_pc2 <- order(res_tab$full_pc02_stand, decreasing = T)
atl_2_samp_pos <- rep(NA, times = length(top_pc2))
atl_2_samp_pos[top_pc2] <- seq(length(top_pc2))

atl_joint_score <- atl_1_samp_pos + atl_2_samp_pos
best_atl_joint <- order(atl_joint_score, decreasing = F)

tot_admix_atl_inds <- which(res_tab$full_admix_k3_ATL > 0.99)

# choose samples that 1) have top "joint" score, 2) are NOT northern 
#   3) 99% ATL according to admix
tmp_atl_inds <- intersect(setdiff(best_atl_joint, north_inds), 
  tot_admix_atl_inds)[1:n_gp]

# make sure is all 4X
sum(res_tab[tmp_atl_inds, ploidy] != '4X') == 0
# [1] TRUE

tmp_atl_names <- res_tab[tmp_atl_inds, samp_name]

# how many accessions
length(unique(samp_meta[VCF_NAME %in% tmp_atl_names, UNI_ACC]))
# 25 accessions

# how many states
length(unique(samp_meta[VCF_NAME %in% tmp_atl_names, STATE]))
# 6 states

# these 60 seem good
tot_atl_names <- tmp_atl_names

# Select 40 across accessions
sub_atl_names <- c()
for(acc_name in unique(samp_meta[VCF_NAME %in% tot_atl_names, UNI_ACC])){
  tmp_names <- intersect(tot_atl_names, 
    samp_meta[UNI_ACC == acc_name, VCF_NAME])
  tmp_choice <- sample(tmp_names, size = 1)
  sub_atl_names <- c(sub_atl_names, tmp_choice)
  }

atl_leftover <- setdiff(tot_atl_names, sub_atl_names)

40-length(sub_atl_names)
#15

length(unique(samp_meta[VCF_NAME %in% atl_leftover, UNI_ACC]))
# 22
length(unique(samp_meta[VCF_NAME %in% atl_leftover, STATE]))
# 5

atl_left_acc <- unique(samp_meta[VCF_NAME %in% atl_leftover, UNI_ACC])
atl_left_state <- unique(samp_meta[VCF_NAME %in% atl_leftover, STATE])

for(stn in atl_left_state){
  left_in_state <- intersect(atl_leftover, samp_meta[STATE == stn, VCF_NAME])
  for(acc_name in sample(unique(
   samp_meta[VCF_NAME %in% left_in_state, UNI_ACC]), size = 2)){
    tmp_names <- intersect(left_in_state,
      samp_meta[UNI_ACC == acc_name, VCF_NAME])
    tmp_choice <- sample(tmp_names, size = 1)
    sub_atl_names <- c(sub_atl_names, tmp_choice)
    atl_left_acc <- setdiff(atl_left_acc, acc_name)
    }
  }

atl_leftover <- setdiff(tot_atl_names, sub_atl_names)

40-length(sub_atl_names)
# 5

for(acc_name in sample(atl_left_acc, size = 5)){
  tmp_names <- intersect(atl_leftover, 
    samp_meta[UNI_ACC == acc_name, VCF_NAME])
  tmp_choice <- sample(tmp_names, size = 1)
  sub_atl_names <- c(sub_atl_names, tmp_choice)
}

length(sub_atl_names)
# 40

fwrite(data.table(samp = tot_atl_names), file = atl_train_full_file,
  sep = '\t', row.names = F, col.names = F)
fwrite(data.table(samp = sub_atl_names), file = atl_train_filt_file,
  sep = '\t', row.names = F, col.names = F)

##########

# Gulf Training Samples

bot_pc2 <- order(res_tab$full_pc02_stand, decreasing = F)

bot_pc2_inds <- bot_pc2[1:n_gp]

# make sure is all 4X
sum(res_tab[bot_pc2_inds, ploidy] != '4X') == 0
#[1] FALSE

# make sure is all 99+% GULF
sum(res_tab[bot_pc2_inds, full_admix_k3_GULF] < 0.99) == 0
#[1] FALSE

# need to find overlap with 4X and ADMIX-GULF
tmp_gulf_filt_inds <- intersect(bot_pc2, 
  intersect(which(res_tab$ploidy == '4X'), 
  which(res_tab$full_admix_k3_GULF > 0.99)))[1:n_gp]

sum(res_tab[tmp_gulf_filt_inds, ploidy] != '4X') == 0
# [1] TRUE

sum(res_tab[tmp_gulf_filt_inds, full_admix_k3_GULF] < 0.99) == 0
# TRUE

summary(res_tab[tmp_gulf_filt_inds, full_pc02_stand])
# a little higher than ideal, but it is what it is
summary(res_tab[bot_pc2_inds, full_pc02_stand])
summary(res_tab[full_admix_k3_GULF > 0.99, full_pc02_stand])

# NOT fully "LOWLAND" in the K=2 Admix results is  the most worrisome 
#  aspect of using GULF samples, but nothing we can do about it
summary(res_tab[tmp_gulf_filt_inds, full_admix_k2_LOW])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.8335  0.8695  0.8935  0.8899  0.9078  0.9236

# how many subgrps
length(unique(res_tab[tmp_gulf_filt_inds, subgrp_v2]))
# 8
res_tab[tmp_gulf_filt_inds, .N, by = subgrp_v2]
# largely GULF_14 and GULF_01, but others, as well

tmp_gulf_names <- res_tab[tmp_gulf_filt_inds, samp_name]

# how many accessions
length(unique(samp_meta[VCF_NAME %in% tmp_gulf_names, UNI_ACC]))
# 37 accessions

# how many states
length(unique(samp_meta[VCF_NAME %in% tmp_gulf_names, STATE]))
# 3 states 
#   1 is SC, which is from samples that were supposed to be removed because
#     is probably ALAMO, however, I'm fine with them being included here
#     because I do think they accurately reflect a "pure" gulf genotype

# these 60 seem good
tot_gulf_names <- tmp_gulf_names

# Select 40 across accessions
sub_gulf_names <- c()
for(acc_name in unique(samp_meta[VCF_NAME %in% tot_gulf_names, UNI_ACC])){
  tmp_names <- intersect(tot_gulf_names,
    samp_meta[UNI_ACC == acc_name, VCF_NAME])
  tmp_choice <- sample(tmp_names, size = 1)
  sub_gulf_names <- c(sub_gulf_names, tmp_choice)
  }

gulf_leftover <- setdiff(tot_gulf_names, sub_gulf_names)

40-length(sub_gulf_names)
#3

length(unique(samp_meta[VCF_NAME %in% gulf_leftover, UNI_ACC]))
# 13
length(unique(samp_meta[VCF_NAME %in% gulf_leftover, STATE]))
# 3

gulf_left_state <- unique(samp_meta[VCF_NAME %in% gulf_leftover, STATE])

for(stn in gulf_left_state){
  left_in_state <- intersect(gulf_leftover, samp_meta[STATE == stn, VCF_NAME])
  for(acc_name in sample(unique(
   samp_meta[VCF_NAME %in% left_in_state, UNI_ACC]), size = 1)){
    tmp_names <- intersect(left_in_state,
      samp_meta[UNI_ACC == acc_name, VCF_NAME])
    tmp_choice <- sample(tmp_names, size = 1)
    sub_gulf_names <- c(sub_gulf_names, tmp_choice)
    }
  }

40-length(sub_gulf_names)
# 0

fwrite(data.table(samp = tot_gulf_names), file = gulf_train_full_file,
  sep = '\t', row.names = F, col.names = F)
fwrite(data.table(samp = sub_gulf_names), file = gulf_train_filt_file,
  sep = '\t', row.names = F, col.names = F)

quit(save = 'no')





