# Look at patterns and overlaps of pop-level/F-val results for
#  MW groups

# GOALS
# WITHIN-GROUP
## - overlap in window-based and SD-based outliers
# BETWEEN GROUPS
## - overlaps in window-based, SD-based, and both between groups
## - distinct SNPs in each set
# BETWEEN SNP sets
## - Window overlaps between GULF and ATL results
## - Top windows/Top SNPs in windows private to each group

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(gridExtra)
### INPUT DATA ###

data_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/'

samp_set_names <- c('MW8X_all', 'MW8X_reg', 'MW4X_intro', 'MW8X_hi', 
  'MW_01', 'MW_03')

pop2_names <- c('GULF', 'ATL')

all_window_suf <- '_into_MW_10SNP_allwindows.txt'
top_window_suf <- '_into_MW_10SNP_topwindows.txt'
all_snp_suf <- '_into_MW_allsnps.txt'
outlier_snp_suf <- '_into_MW_outliersnps.txt'
window_snp_suf <- '_into_MW_topwindowsnps.txt'
###########################
# Tally of results

gulf_outlier_snp_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_strict_Fval_', 'GULF', outlier_snp_suf, 
    sep = '')
  gulf_outlier_snp_list[[gn]] <- fread(tmp_file)
}
unlist(lapply(gulf_outlier_snp_list, function(x) nrow(x)))
#   MW8X_all   MW8X_reg MW4X_intro    MW8X_hi      MW_01      MW_03 
#         45         50          3          5          6         26

gulf_window_snp_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_strict_Fval_', 'GULF', window_snp_suf,   
    sep = '')
  gulf_window_snp_list[[gn]] <- fread(tmp_file)
}
unlist(lapply(gulf_window_snp_list, function(x) nrow(x)))
#   MW8X_all   MW8X_reg MW4X_intro    MW8X_hi      MW_01      MW_03 
#        137        145         45         32        125        137

atl_outlier_snp_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_strict_Fval_', 'ATL', outlier_snp_suf,   
    sep = '')
  atl_outlier_snp_list[[gn]] <- fread(tmp_file)
}
unlist(lapply(atl_outlier_snp_list, function(x) nrow(x)))
#   MW8X_all   MW8X_reg MW4X_intro    MW8X_hi      MW_01      MW_03 
#         10         11          2          0          2          8 

atl_window_snp_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_strict_Fval_', 'ATL', window_snp_suf,    
    sep = '')
  atl_window_snp_list[[gn]] <- fread(tmp_file)
}
unlist(lapply(atl_window_snp_list, function(x) nrow(x)))
#  MW8X_all   MW8X_reg MW4X_intro    MW8X_hi      MW_01      MW_03 
#       102        101         25         25         93         87 

###########################################
# Overlap within sample sets
gulf_outlier_snpname_list <- list()
gulf_window_snpname_list <- list()
gulf_overlap_list <- list()
for(gn in samp_set_names){
  tmp_out_tab <- gulf_outlier_snp_list[[gn]]
  tmp_outlier_snp_names <- paste(tmp_out_tab$CHR, tmp_out_tab$POS, sep = '_')
  gulf_outlier_snpname_list[[gn]] <- tmp_outlier_snp_names
  tmp_wind_tab <- gulf_window_snp_list[[gn]]
  tmp_wind_snp_names <- paste(tmp_wind_tab$CHR, tmp_wind_tab$POS, sep = '_')
  gulf_window_snpname_list[[gn]] <- tmp_wind_snp_names
  gulf_overlap_list[[gn]] <- intersect(tmp_outlier_snp_names, 
    tmp_wind_snp_names)
 }
unlist(lapply(gulf_overlap_list, length))
#   MW8X_all   MW8X_reg MW4X_intro    MW8X_hi      MW_01      MW_03 
#         14         13          0          0          2          8 

atl_outlier_snpname_list <- list()
atl_window_snpname_list <- list()
atl_overlap_list <- list()
for(gn in samp_set_names){
  tmp_out_tab <- atl_outlier_snp_list[[gn]]
  tmp_outlier_snp_names <- paste(tmp_out_tab$CHR, tmp_out_tab$POS, sep = '_')
  atl_outlier_snpname_list[[gn]] <- tmp_outlier_snp_names
  tmp_wind_tab <- atl_window_snp_list[[gn]]
  tmp_wind_snp_names <- paste(tmp_wind_tab$CHR, tmp_wind_tab$POS, sep = '_')
  atl_window_snpname_list[[gn]] <- tmp_wind_snp_names
  atl_overlap_list[[gn]] <- intersect(tmp_outlier_snp_names,
    tmp_wind_snp_names)
}
unlist(lapply(atl_overlap_list, length))
#   MW8X_all   MW8X_reg MW4X_intro    MW8X_hi      MW_01      MW_03 
#          2          2          0          0          1          1 

############################################
### Overlap between sample sets
# Gulf
table(table(unlist(gulf_outlier_snpname_list)))
#  1  2  3  4  5 
# 23 32 13  1  1
sum(table(table(unlist(gulf_outlier_snpname_list))))
# 70

table(table(unlist(gulf_window_snpname_list)))
#   1   2   3   4   5 
# 139  41  48  39  20
sum(table(table(unlist(gulf_window_snpname_list))))
# 287

table(table(c(unlist(gulf_outlier_snpname_list), 
  unlist(gulf_window_snpname_list))))
#   1   2   3   4   5   6   7   8 
# 153  57  54  41  19   4   4   2
sum(table(table(c(unlist(gulf_outlier_snpname_list), 
  unlist(gulf_window_snpname_list)))))
# [1] 334
# 334-153 = 181 SNPs in more than 1 result set

#which(table(unlist(gulf_outlier_snpname_list)) == 6)
#which(table(unlist(gulf_window_snpname_list)) == 6)

# Atlantic
table(table(unlist(atl_outlier_snpname_list)))
# 1 2 3 4 
# 6 4 5 1
sum(table(table(unlist(atl_outlier_snpname_list))))
# 16

table(table(unlist(atl_window_snpname_list)))
#  1  2  3  4  5 
# 66 16 45 20 24 
sum(table(table(unlist(atl_window_snpname_list))))
# 171

table(table(c(unlist(atl_outlier_snpname_list), 
  unlist(atl_window_snpname_list))))
#  1  2  3  4  5  6 
# 71 19 46 22 25  1
sum(table(table(c(unlist(atl_outlier_snpname_list), 
  unlist(atl_window_snpname_list)))))
# [1] 184
# 184-71 = 113 SNPs in more than 1 result set

# CONTINUE FROM HERE


##########################################################
### Comparison of MW8X_all and MW8X_reg
length(intersect(gulf_outlier_snpname_list[['MW8X_all']], 
  gulf_outlier_snpname_list[['MW8X_reg']]))
# 703
length(intersect(gulf_outlier_snpname_list[['MW8X_all']], 
  gulf_outlier_snpname_list[['MW8X_reg']]))/length(
  gulf_outlier_snpname_list[['MW8X_all']])
# [1] 0.9348404
# 93.4% of GULF "outlier" SNPs from MW8X_all are also included in MW8X_reg 
#   results

length(intersect(gulf_window_snpname_list[['MW8X_all']], 
  gulf_window_snpname_list[['MW8X_reg']]))
# 401
length(intersect(gulf_window_snpname_list[['MW8X_all']], 
  gulf_window_snpname_list[['MW8X_reg']])) / length(
  gulf_window_snpname_list[['MW8X_all']])
# [1] 0.926097
# 92.6% of GULF window SNPs from MW8X_all are also included in MW8X_reg results

# Atlantic
length(intersect(atl_outlier_snpname_list[['MW8X_all']], 
  atl_outlier_snpname_list[['MW8X_reg']]))
# 568
length(intersect(atl_outlier_snpname_list[['MW8X_all']],
  atl_outlier_snpname_list[['MW8X_reg']]))/length(
  atl_outlier_snpname_list[['MW8X_all']])
# [1] 0.9530201
# 95.3% of ATL outlier SNPs from MW8Xall also included in MW8X_reg results

length(intersect(atl_window_snpname_list[['MW8X_all']], 
  atl_window_snpname_list[['MW8X_reg']]))
# 355
length(intersect(atl_window_snpname_list[['MW8X_all']],
  atl_window_snpname_list[['MW8X_reg']])) / length(
  atl_window_snpname_list[['MW8X_all']])
# [1] 0.9244792
# 92.4% of ATL window SNPs from MW8X_all also included in MW8X_reg results
#  - is a bit lower than for "outlier" SNPs

############################################
### Look at distribution of F_val results
## SNPs
gulf_all_snp_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_Fval_', 'GULF', all_snp_suf, 
    sep = '')
  gulf_all_snp_list[[gn]] <- fread(tmp_file)
}

atl_all_snp_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_Fval_', 'ATL', all_snp_suf,    
    sep = '')
  atl_all_snp_list[[gn]] <- fread(tmp_file)
}

gulf_fval_dist_gg_list <- list()
for(gn in samp_set_names){
  tmp_tab <- gulf_all_snp_list[[gn]]
  tmp_sd <- sd(tmp_tab[, F_subgrp_v_pop1], na.rm = T)
  tmp_cut <- mean(tmp_tab[, F_subgrp_v_pop1], na.rm = T) + (3*tmp_sd)
  tmp_gg <- ggplot(tmp_tab) +
    geom_density(aes_string(x = 'F_subgrp_v_pop1'), fill = 'grey50') +
    xlab('F_val') +
    ggtitle(paste(gn, 'Gulf into MW SNP F_val\ndistribution', sep = ' ')) +
    geom_vline(xintercept = tmp_cut, color = 'red', linetype = 'dotted')
  gulf_fval_dist_gg_list[[gn]] <- tmp_gg
}

gulf_snp_fval_dist_pdf <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'GULFintoMW_SNP_Fval_distribution.pdf',
  sep = '')

pdf(gulf_snp_fval_dist_pdf, width = 4*3, height = 4*2)
grid.arrange(grobs = gulf_fval_dist_gg_list, nrow = 2)
dev.off()

atl_fval_dist_gg_list <- list()
for(gn in samp_set_names){
  tmp_tab <- atl_all_snp_list[[gn]]
  tmp_sd <- sd(tmp_tab[, F_subgrp_v_pop1], na.rm = T)
  tmp_cut <- mean(tmp_tab[, F_subgrp_v_pop1], na.rm = T) + (3*tmp_sd)
  tmp_gg <- ggplot(tmp_tab) +
    geom_density(aes_string(x = 'F_subgrp_v_pop1'), fill = 'grey50') +
    xlab('F_val') +
    ggtitle(paste(gn, 'Atlantic into MW SNP F_val\ndistribution', sep = ' ')) +
    geom_vline(xintercept = tmp_cut, color = 'red', linetype = 'dotted')
  atl_fval_dist_gg_list[[gn]] <- tmp_gg
}

atl_snp_fval_dist_pdf <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'ATLintoMW_SNP_Fval_distribution.pdf',
  sep = '')

pdf(atl_snp_fval_dist_pdf, width = 4*3, height = 4*2)
grid.arrange(grobs = atl_fval_dist_gg_list, nrow = 2)
dev.off()

## Windows
gulf_all_window_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_Fval_', 'GULF', all_window_suf,
    sep = '')
  gulf_all_window_list[[gn]] <- fread(tmp_file)
}

gulf_fval_window_gg_list <- list()
for(gn in samp_set_names){
  tmp_tab <- gulf_all_window_list[[gn]]
  tmp_sd <- sd(tmp_tab[, F_SUB_V_POP1], na.rm = T)
  tmp_cut <- mean(tmp_tab[, F_SUB_V_POP1], na.rm = T) + (3*tmp_sd)
  tmp_gg <- ggplot(tmp_tab) +
    geom_density(aes_string(x = 'F_SUB_V_POP1'), fill = 'grey50') +
    xlab('F_val') +
    ggtitle(paste(gn, 'Gulf into MW 10bp window\nF_val distribution', 
      sep = ' ')) +
    geom_vline(xintercept = tmp_cut, color = 'red', linetype = 'dotted')
  gulf_fval_window_gg_list[[gn]] <- tmp_gg
}

gulf_window_fval_dist_pdf <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'GULFintoMW_Window_Fval_distribution.pdf',
  sep = '')

pdf(gulf_window_fval_dist_pdf, width = 4*3, height = 4*2)
grid.arrange(grobs = gulf_fval_window_gg_list, nrow = 2)
dev.off()

atl_all_window_list <- list()
for(gn in samp_set_names){
  tmp_file <- paste(data_dir, gn, '_Fval_', 'ATL', all_window_suf,
    sep = '')
  atl_all_window_list[[gn]] <- fread(tmp_file)
}

atl_fval_window_gg_list <- list()
for(gn in samp_set_names){
  tmp_tab <- atl_all_window_list[[gn]]
  tmp_sd <- sd(tmp_tab[, F_SUB_V_POP1], na.rm = T)
  tmp_cut <- mean(tmp_tab[, F_SUB_V_POP1], na.rm = T) + (3*tmp_sd)
  tmp_gg <- ggplot(tmp_tab) +
    geom_density(aes_string(x = 'F_SUB_V_POP1'), fill = 'grey50') +
    xlab('F_val') +
    ggtitle(paste(gn, 'Atlantic into MW 10bp window\nF_val distribution', 
      sep = ' ')) +
    geom_vline(xintercept = tmp_cut, color = 'red', linetype = 'dotted')
  atl_fval_window_gg_list[[gn]] <- tmp_gg
}

atl_window_fval_dist_pdf <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'ATLintoMW_Window_Fval_distribution.pdf',
  sep = '')

pdf(atl_window_fval_dist_pdf, width = 4*3, height = 4*2)
grid.arrange(grobs = atl_fval_window_gg_list, nrow = 2)
dev.off()






