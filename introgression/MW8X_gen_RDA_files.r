# Generate files for RDA analysis of introgressed regions in MW-8X and
#  subgroups

# Groups to generate:
#  All MW-8X
#  MW_01
#  MW_03

# Comparisons
# MW vs GULF
# MW vs ATL

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATE ###

# info about samples
samp_info_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'sg_8X_result_tabs/', 'natv2filt_res_tab_v4.0.txt', sep = '')
samp_info <- fread(samp_info_file)

# introgression analysis results
mw_v_g_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/', 
  'MW_8X_GULF_into_MW_10SNP_window_results.rds', sep = '')
mw_v_g_list <- readRDS(mw_v_g_file)

mw_v_a_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/', 
  'MW_8X_ATL_into_MW_10SNP_window_results.rds', sep = '')
mw_v_a_list <- readRDS(mw_v_a_file)

# allele state file
mw_v_g_state_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/MW_v_GULF_allele_states.txt', sep = '')
mw_v_g_state <- fread(mw_v_g_state_file)

mw_v_a_state_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/MW_v_ATL_allele_states.txt', sep = '')
mw_v_a_state <- fread(mw_v_a_state_file)

### SET OUTPUT ###
out_dir <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', sep = '')

mw8X_mw_v_g_rda_out <- paste(out_dir, 'MW8X_all_GULF_into_MW.RDA_genos.rds',
  sep = '')
mw01_mw_v_g_rda_out <- paste(out_dir, 'MW01_GULF_into_MW.RDA_genos.rds',
  sep = '')
mw03_mw_v_g_rda_out <- paste(out_dir, 'MW03_GULF_into_MW.RDA_genos.rds',
  sep = '')

mw8X_mw_v_a_rda_out <- paste(out_dir, 'MW8X_all_ATL_into_MW.RDA_genos.rds',
  sep = '')
mw01_mw_v_a_rda_out <- paste(out_dir, 'MW01_ATL_into_MW.RDA_genos.rds',
  sep = '')
mw03_mw_v_a_rda_out <- paste(out_dir, 'MW03_ATL_into_MW.RDA_genos.rds',
  sep = '')


##########

# Subgroup names
MW_01_names <- samp_info[sub_grp == 'MW_01', samp_name]
MW_03_names <- samp_info[sub_grp == 'MW_03', samp_name]

# MW vs GULF
full_mw_v_g_rda_mat <- matrix(
  unlist(lapply(mw_v_g_list, function(x) x[['rda_geno_vec']])), byrow = T,
  nrow = length(mw_v_g_list))

rownames(full_mw_v_g_rda_mat) <- names(mw_v_g_list)
mw_v_g_snp_names <- paste(mw_v_g_state$CHR, mw_v_g_state$POS, sep = '_')
colnames(full_mw_v_g_rda_mat) <- mw_v_g_snp_names
saveRDS(full_mw_v_g_rda_mat, file = mw8X_mw_v_g_rda_out)

mw_v_g_mw01_inds <- which(rownames(full_mw_v_g_rda_mat) %in% MW_01_names)
mw01_mw_v_g_rda_mat <- full_mw_v_g_rda_mat[mw_v_g_mw01_inds, ]
saveRDS(mw01_mw_v_g_rda_mat, file = mw01_mw_v_g_rda_out)

mw_v_g_mw03_inds <- which(rownames(full_mw_v_g_rda_mat) %in% MW_03_names)
mw03_mw_v_g_rda_mat <- full_mw_v_g_rda_mat[mw_v_g_mw03_inds, ]
saveRDS(mw03_mw_v_g_rda_mat, file = mw03_mw_v_g_rda_out)

# MW vs ATL
full_mw_v_a_rda_mat <- matrix(
  unlist(lapply(mw_v_a_list, function(x) x[['rda_geno_vec']])), byrow = T,
  nrow = length(mw_v_a_list))

rownames(full_mw_a_g_rda_mat) <- names(mw_v_a_list)
mw_v_a_snp_names <- paste(mw_v_a_state$CHR, mw_v_a_state$POS, sep = '_')
colnames(full_mw_v_a_rda_mat) <- mw_v_a_snp_names
saveRDS(full_mw_v_a_rda_mat, file = mw8X_mw_v_a_rda_out)

mw_v_a_mw01_inds <- which(rownames(full_mw_v_a_rda_mat) %in% MW_01_names)
mw01_mw_v_a_rda_mat <- full_mw_v_a_rda_mat[mw_v_a_mw01_inds, ]
saveRDS(mw01_mw_v_a_rda_mat, file = mw01_mw_v_a_rda_out)

mw_v_a_mw03_inds <- which(rownames(full_mw_v_a_rda_mat) %in% MW_03_names)
mw03_mw_v_a_rda_mat <- full_mw_v_a_rda_mat[mw_v_a_mw03_inds, ]
saveRDS(mw03_mw_v_a_rda_mat, file = mw03_mw_v_a_rda_out)

quit(save = 'no')






