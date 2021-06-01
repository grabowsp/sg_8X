# Compare levels of inferred GULF and ATL introgression in samples to 
#   differentiate the sources of lowland alleles in different MW groups

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)

### INPUT DATE ###
samp_info_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'sg_8X_result_tabs/', 'natv2filt_res_tab_v4.0.txt', sep = '')
samp_info <- fread(samp_info_file)

train_samp_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/', 'introgression_samps_and_groups.txt', sep = '')
train_samps <- fread(train_samp_file)

# introgression allele results

mw_g_result_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_GULF_into_MW_10SNP_window_results.rds', sep = '')
mw_g_res_list <- readRDS(mw_g_result_file)

mw_a_result_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_ATL_into_MW_10SNP_window_results.rds', sep = '')
mw_a_res_list <- readRDS(mw_a_result_file)

# SNP postion info
mw_v_g_pos_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/', 'MW_v_GULF_keep_pos.txt', sep = '')
mw_g_pos <- fread(mw_v_g_pos_file)

mw_v_a_pos_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/', 'MW_v_ATL_keep_pos.txt', sep = '')
mw_a_pos <- fread(mw_v_a_pos_file)

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

########
# MW vs Gulf results
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

mean_mw_v_g_rda <- unlist(lapply(mw_g_gw_list, function(x) mean(x$rda_geno)))
mean_mw_v_g_score <- unlist(lapply(mw_g_gw_list, 
  function(x) mean(x$geno_score)))

# MW vs ATL results
mw_a_gw_list <- list()
for(mws in names(mw_a_res_list)){
  tmp_tab <- data.table(mw_a_pos,
    geno_state = mw_a_res_list[[mws]]$geno_state_vec,
    geno_cat = mw_a_res_list[[mws]]$geno_cat_vec,
    rda_geno = mw_a_res_list[[mws]]$rda_geno_vec)
  tmp_score_vec <- rep(NA, times = nrow(tmp_tab))
  for(cn in unique(tmp_tab$geno_cat)){
    tmp_score_vec[which(tmp_tab$geno_cat == cn)] <- intro_wt_vec[[cn]]
  }
  tmp_tab[ , geno_score := tmp_score_vec]
  mw_a_gw_list[[mws]] <- tmp_tab
}

mean_mw_v_a_rda <- unlist(lapply(mw_a_gw_list, function(x) mean(x$rda_geno)))
mean_mw_v_a_score <- unlist(lapply(mw_a_gw_list,
  function(x) mean(x$geno_score)))

sum(names(mw_a_gw_list) != names(mw_g_gw_list))
# 0 - same sample order list

# Generate sample table
mean_val_tab <- data.table(samp_name = names(mw_g_res_list),
  MWvG_rda = mean_mw_v_g_rda, 
  MWvG_score = mean_mw_v_g_score,
  MWvA_rda = mean_mw_v_a_rda,
  MWvA_score = mean_mw_v_a_score)

# assign genepools
mean_val_tab[, subgrp := as.character(NA)]
for(i in seq(nrow(mean_val_tab))){
  tmp_info_ind <- which(samp_info$samp_name == mean_val_tab$samp_name[i])
  mean_val_tab[i, subgrp := samp_info$sub_grp[tmp_info_ind]]
}

# NOTE: I think we need to split up MW_01 because of the high values for some
#  samples, but I want to wait and look at the MW_v_A results first
# On second thought, maybe not... There isn't the same clean break for the 
#   ATL results.
# Still probably worth doing

mean_val_tab[intersect(which(mean_val_tab$subgrp == 'MW_01'), 
  which(mean_val_tab$MWvG_rda > 0.21)), subgrp := 'MW_01_hi']

# condolidate groups - hi introgression 8X, reg 8X, 4X
mean_val_tab[, sub_2 := '4X']
mean_val_tab[grep('MW_01|MW_02|MW_03|MW_04|MW_05', mean_val_tab$subgrp), 
  sub_2 := '8X']
mean_val_tab[which(mean_val_tab$subgrp == 'MW_01_hi'), sub_2 := '8X_hi']

# Split out the 4X with signs of ingrogression
mean_val_tab[, sub_3 := mean_val_tab$sub_2]
mean_val_tab[grep('MW_06|MW_07', mean_val_tab$subgrp), 
  sub_3 := '4X_intro']


#tapply(mean_val_tab$MWvG_rda, mean_val_tab$subgrp, mean)

## Look at correlation of values

cor(mean_mw_v_g_rda, mean_mw_v_g_score)
#[1] 0.9962535
# should plot this

gg_cor <- ggplot(mean_val_tab, aes(x = MWvG_rda, y = MWvG_score, 
    color = subgrp)) +
  geom_point() +
  labs(title = 'MW_v_Gulf mean RDA vs MWvGulf mean score')

cor_plot_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_MWvGulf_rda_score_cor.pdf', sep = '')

pdf(cor_plot_out, height = 5, width = 5)
gg_cor
dev.off()

cor(mean_mw_v_a_rda, mean_mw_v_a_score)
# [1] 0.9972682

cor(mean_mw_v_a_rda, mean_mw_v_g_rda)
# [1] 0.9807359

lm_gva <- lm(MWvA_rda ~ MWvG_rda, 
  data = mean_val_tab[-which(mean_val_tab$subgrp == 'MW_01_hi'), ])

gg_g_a_rda <- ggplot(mean_val_tab, aes(x = MWvG_rda, y = MWvA_rda, 
  color = subgrp)) +
  geom_point() +
  geom_abline(intercept = lm_gva$coefficients[1], 
    slope = lm_gva$coefficients[2], color = 'gray50', linetype = 'dashed') +
  labs(title = 'MWvGULF rda vs MWvATL rda')

comp_score_plot_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_GulfvATL_rda.pdf', sep = '')

pdf(comp_score_plot_out, height = 5, width = 6)
gg_g_a_rda
dev.off()

### boxplots of rda (Pop2 allele counts) for each genepool
gg_g_rda <- ggplot(mean_val_tab, aes(x = subgrp, y = MWvG_rda)) +
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = 'Mean Gulf alleles in MW subgroups')

gulf_rda_box_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_Gulf_rda_box.pdf', sep = '')

pdf(gulf_rda_box_out, height = 5, width = 5)
gg_g_rda
dev.off()

gg_g_rda_2 <- ggplot(mean_val_tab, aes(x = sub_2, y = MWvG_rda)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('MW groups') +
  labs(title = 'Mean Gulf alleles in MW split by ploidy')

gulf_rda_box_out_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_Gulf_rda_box_v2.pdf', sep = '')

pdf(gulf_rda_box_out_2, height = 5, width = 5)
gg_g_rda_2
dev.off()

gg_g_rda_3 <- ggplot(mean_val_tab, aes(x = sub_3, y = MWvG_rda)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('MW groups') +
  labs(title = 'Mean Gulf alleles in MW split by ploidy')

gulf_rda_box_out_3 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_Gulf_rda_box_v3.pdf', sep = '')

pdf(gulf_rda_box_out_3, height = 5, width = 5)
gg_g_rda_3
dev.off()



### Atlantic
gg_a_rda <- ggplot(mean_val_tab, aes(x = subgrp, y = MWvA_rda)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = 'Mean Atlantic alleles in MW subgroups')

atl_rda_box_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_ATL_rda_box.pdf', sep = '')

pdf(atl_rda_box_out, height = 5, width = 5)
gg_a_rda
dev.off()

gg_a_rda_2 <- ggplot(mean_val_tab, aes(x = sub_2, y = MWvA_rda)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('MW groups') +
  labs(title = 'Mean Atlantic alleles in MW split by ploidy')

atl_rda_box_out_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_ATL_rda_box_v2.pdf', sep = '')

pdf(atl_rda_box_out_2, height = 5, width = 5)
gg_a_rda_2
dev.off()

gg_a_rda_3 <- ggplot(mean_val_tab, aes(x = sub_3, y = MWvA_rda)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('MW groups') +
  labs(title = 'Mean Atlantic alleles in MW split by ploidy')

atl_rda_box_out_3 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_ATL_rda_box_v3.pdf', sep = '')

pdf(atl_rda_box_out_3, height = 5, width = 5)
gg_a_rda_3
dev.off()


###
# Gulf rda vs ATL rda by genepool

mean_rda_dif_1 <- ((mean_val_tab$MWvG_rda - mean_val_tab$MWvA_rda)/
  mean_val_tab$MWvA_rda)

mean_val_tab[, rda_dif_1 := mean_rda_dif_1]


gg_rda_dif_box <- ggplot(mean_val_tab, aes(x = subgrp, y = rda_dif_1)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

rda_dif_box_plot <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_rda_diff_boxplot.pdf', sep = '')

pdf(rda_dif_box_plot, height = 5, width = 5)
gg_rda_dif_box
dev.off()

### Statisitcal tests for MW_01 vs MW_03

rda_dif_test <- t.test(
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_01'), rda_dif_1],
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_03'), rda_dif_1], 
  var.equal = F, alternative = 'greater')
# simple t-test 
# p-value = 4.8e-7 for test that MW_01 is greater (more GULF) than MW_03
# So, amount of "GULF" is different between the two, but not sure can say that
#  ATL has overall higher amount in MW_03

rda_non_param_test <- wilcox.test(
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_01'), rda_dif_1],
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_03'), rda_dif_1], 
  alternative = 'greater')
# Mann-Whitney U test; non-parameteric test
# p-value = 1.233e-9

### Statistical tests for MW_01 vs non-introgressed 4X
cont_4X_inds <- grep('MW_08|MW_09|MW_10|MW_12|MW_13|MW_14', mean_val_tab$subgrp)

rda_MW01v4X_test <-t.test(
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_01'), rda_dif_1],
  mean_val_tab[cont_4X_inds, rda_dif_1],
  var.equal = F)
# t-test for MW_01 vs non-introgressed 4X
# p-value = 0.002334
# 95 percent confidence interval for difference between means:
# 0.008598834 0.035304169
# sample estimates:
#   mean of x    mean of y; x=MW_01, y=4X_control
#  0.002480684 -0.019470817

rda_np_MW01v4X_test <- wilcox.test(
 mean_val_tab[which(mean_val_tab$subgrp == 'MW_01'), rda_dif_1],
  mean_val_tab[cont_4X_inds, rda_dif_1])
# Mann-Whitney U test; non-parameteric test
# p-value = 0.001665

### Statistical tests for MW_03 vs non-introgress 4X
rda_MW03v4X_test <-t.test(
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_03'), rda_dif_1],
  mean_val_tab[cont_4X_inds, rda_dif_1],
  var.equal = F)
# t-test for MW_03 vs non-introgressed 4X
# p-value = 3.63e-7
# 95 percent confidence interval:
#  -0.02693540 -0.01251704
# sample estimates:
#   mean of x   mean of y ; x=MW_03; y=4X_control 
# -0.03919704 -0.01947082

rda_np_MW03v4X_test <- wilcox.test(
  mean_val_tab[which(mean_val_tab$subgrp == 'MW_03'), rda_dif_1],
  mean_val_tab[cont_4X_inds, rda_dif_1])
# Mann-Whitney U test; non-parametric test
# p-value = 6.77e-7
### Look at MW_01_hi vs MW_01 in PCA results to see if distinct
mw_01_hi_names <- mean_val_tab[which(mean_val_tab$subgrp == 'MW_01_hi'), 
  samp_name]
mw_01_names <- mean_val_tab[which(mean_val_tab$subgrp == 'MW_01'), 
  samp_name]

summary(samp_info[which(samp_info$samp_name %in% mw_01_hi_names), 
  full_pc01_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7006  0.7352  0.7505  0.7507  0.7612  0.8032 

summary(samp_info[which(samp_info$samp_name %in% mw_01_names), full_pc01_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7786  0.8199  0.8789  0.8585  0.8938  0.9023 

summary(samp_info[which(samp_info$samp_name %in% mw_01_hi_names), 
  full_pc02_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5806  0.6312  0.6431  0.6416  0.6478  0.6883

summary(samp_info[which(samp_info$samp_name %in% mw_01_names), full_pc02_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.7036  0.7368  0.7683  0.7552  0.7767  0.7831 

# MW_01 and MW_01_hi do NOT overlap on PC2

### MW pca values
summary(samp_info[which(samp_info$samp_name %in% mw_01_hi_names), 
  MW_pc01_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.1214  0.1818  0.2228  0.2173  0.2483  0.3000 

summary(samp_info[which(samp_info$samp_name %in% mw_01_names),
  MW_pc01_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1574  0.2308  0.2711  0.2471  0.2773  0.3301 

summary(samp_info[which(samp_info$samp_name %in% mw_01_hi_names),
  MW_pc02_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1395  0.2260  0.2632  0.2528  0.2998  0.3252 

summary(samp_info[which(samp_info$samp_name %in% mw_01_names),
  MW_pc02_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1077  0.1840  0.2087  0.2038  0.2351  0.2562 

# no clear distinction between MW_01_hi and MW_01 in MW PCA, at least from 
#   what I can tell

### Look at number of homozygous genotypes in MW_03 vs MW_06 and MW_07

mw03_names <- samp_info[which(samp_info$sub_grp == 'MW_03'), samp_name]

tmp_ind <- which(names(mw_g_res_list) == mw03_names[1])

table(mw_g_res_list[[tmp_ind]]$rda_geno_vec)

mw06_names <- samp_info[which(samp_info$sub_grp == 'MW_06'), samp_name]

tmp_ind_2 <- which(names(mw_g_res_list) == mw06_names[1])
table(mw_g_res_list[[tmp_ind_2]]$rda_geno_vec)

mw07_names <- samp_info[which(samp_info$sub_grp == 'MW_07'), samp_name]
tmp_ind_3 <- which(names(mw_g_res_list) == mw07_names[1])
table(mw_g_res_list[[tmp_ind_3]]$rda_geno_vec)

### Look at distribution of training samples
samp_info[which(samp_info$samp_name %in% 
  train_samps[which(train_samps$GROUP == 'Midwest'), LIB]), .N, by = sub_grp]

#   sub_grp  N
#1:   MW_14  8
#2:   MW_09  5
#3:   MW_10 16
#4:   MW_08  4
#5:   MW_13  3
#6:   MW_12  1






