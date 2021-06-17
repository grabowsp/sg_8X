# adjust subgroups to include MW_01_hi and subgroup groupings for plotting

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###

samp_info_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'sg_8X_result_tabs/', 'natv2filt_res_tab_v4.0.txt', sep = '')
samp_info <- fread(samp_info_file)

train_samp_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/', 'introgression_samps_and_groups.txt', sep = '')
train_samps <- fread(train_samp_file)

mw_g_result_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_GULF_into_MW_10SNP_window_results.rds', sep = '')
mw_g_res_list <- readRDS(mw_g_result_file)

# SNP postion info
mw_v_g_pos_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MW_v_GULF_keep_pos.txt', sep = '')
mw_g_pos <- fread(mw_v_g_pos_file)

### SET OUTPUT ###
new_samp_info_out <- samp_info_file <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'sg_8X_result_tabs/', 'natv2filt_res_tab_v4.1.txt', sep = '')

### SET VARIABLES ###
### Introgression weights for calculating "score"
#pop1_only_weight <- -1
#pop1_mainly_weight <- -0.5
pop2_only_weight <- 1
pop2_mainly_weight <- 0.5

intro_wt_vec <- c()
intro_wt_vec$cat_e <- 2*pop2_only_weight
intro_wt_vec$cat_f <- 2*pop2_mainly_weight
intro_wt_vec$cat_g <- pop2_only_weight
intro_wt_vec$cat_h <- pop2_mainly_weight
intro_wt_vec$cat_i <- pop2_only_weight
intro_wt_vec$cat_j <- pop2_mainly_weight
intro_wt_vec$cat_k <- pop2_only_weight
intro_wt_vec$cat_l <- pop2_mainly_weight
intro_wt_vec[c('cat_a', 'cat_b', 'cat_c', 'cat_d', 'cat_m', 'cat_n')] <- 0

#######

mw_g_gw_list <- list()

for(mws in names(mw_g_res_list)){
  tmp_tab <- data.table(mw_g_pos,
    geno_state = mw_g_res_list[[mws]]$geno_state_vec,
    geno_cat = mw_g_res_list[[mws]]$geno_cat_vec,
    rda_geno = mw_g_res_list[[mws]]$rda_geno_vec)
  tmp_score_vec <- rep(NA, times = nrow(tmp_tab))
  for(cn in unique(tmp_tab$geno_cat)){
    tmp_score_vec[which(tmp_tab$geno_cat == cn)] <- intro_wt_vec[[cn]]
  }
  tmp_tab[ , geno_score := tmp_score_vec]
  mw_g_gw_list[[mws]] <- tmp_tab
}

# Split off MW_01_hi
mean_mw_v_g_rda <- unlist(lapply(mw_g_gw_list, function(x) mean(x$rda_geno)))
mw8X_hi_names <- names(mean_mw_v_g_rda)[which(mean_mw_v_g_rda > 0.21)]

mw8X_names <- samp_info[
  grep('MW_01|MW_02|MW_03|MW_04|MW_05', samp_info$sub_grp), samp_name]

# add MW_01_hi to subgrp_v2
samp_info[, subgrp_v2 := samp_info$sub_grp]
samp_info[which(samp_info$samp_name %in% mw8X_hi_names), 
  subgrp_v2 := 'MW_01_hi']

# make MW_grp_1, which makes 4 MW groups: MW_4X, Mw_4X_intro, MW_8X, MW_8X_hi
samp_info[, MW_grp_1 := as.character(NA)]
samp_info[grep('MW_', samp_info$subgrp_v2), MW_grp_1 := 'MW_4X']
samp_info[grep('MW_06|MW_07', samp_info$subgrp_v2), MW_grp_1 := 'MW_4X_intro']
samp_info[which(samp_info$samp_name %in% mw8X_names), MW_grp_1 := 'MW_8X']
samp_info[which(samp_info$samp_name %in% mw8X_hi_names), MW_grp_1 := 'MW_8X_hi']

fwrite(samp_info, file = new_samp_info_out, sep = '\t')


