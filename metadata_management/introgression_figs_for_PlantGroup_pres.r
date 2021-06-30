# Figures of introgression results for JGI Plant Group presentation, Jul 1 2021

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)

### INPUT DATE ###
samp_info_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'sg_8X_result_tabs/', 'natv2filt_res_tab_v4.1.txt', sep = '')
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

### SET OUTPUTS ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/'


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
#mean_val_tab[, subgrp_v2 := as.character(NA)]
for(i in seq(nrow(mean_val_tab))){
  tmp_info_ind <- which(samp_info$samp_name == mean_val_tab$samp_name[i])
  mean_val_tab[i, subgrp_v2 := samp_info$subgrp_v2[tmp_info_ind]]
  mean_val_tab[i, ploidy := samp_info$ploidy[tmp_info_ind]]
}

# assign introgression classes
#length(which(mean_val_tab[,MWvG_rda] < 0.1))

# need to adjust
mean_val_tab[, Lowland_class := as.character(NA)]
mean_val_tab[MWvA_rda < 0.1, Lowland_class := 'Low']
mean_val_tab[
  which(mean_val_tab$MWvA_rda >= 0.1 & mean_val_tab$MWvG_rda < 0.22), 
  Lowland_class := 'Moderate']
mean_val_tab[MWvG_rda > 0.22, Lowland_class := 'High']

mean_val_tab[, Lowland_class := factor(mean_val_tab$Lowland_class, 
  levels = c('Low', 'Moderate', 'High'))]

### Plot of Gulf vs Atlantic scores to show introgression classes

lm_gva <- lm(MWvA_rda ~ MWvG_rda, 
  data = mean_val_tab[-which(mean_val_tab$Lowland_class == 'High'), ])

gg_atl_gulf_dotplot <- ggplot(mean_val_tab, aes(x = MWvG_rda, y = MWvA_rda, 
    color = Lowland_class)) +
  geom_point() +
  geom_abline(intercept = lm_gva$coefficients[1],
    slope = lm_gva$coefficients[2], color = 'gray50', linetype = 'dashed') +
  labs('Score for Gulf vs Atlantic alleles in Midwest Samples') +
  xlab('Gulf Alleles Score') +
  ylab('Atlantic Alleles Score')

GvA_dot_out <- paste(out_dir, 'GvA_lineplot_classes.pdf', sep = '')

pdf(GvA_dot_out, height = 4, width = 5)
gg_atl_gulf_dotplot
dev.off()

### boxplots of rda (labeled as "scores" here) for each class
# 3 classes (no ploidy info)
gg_a_rda <- ggplot(mean_val_tab, aes(x = Lowland_class, y = MWvA_rda, 
    fill = Lowland_class)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab('Lowland Allele Class') +
  ylab('Atlantic Alleles Score') +
  ylim(0, 0.3)

atl_box_out_1 <- paste(out_dir, 'Atl_allele_score_box_1.pdf', sep = '')

pdf(atl_box_out_1, height = 5, width = 5)
gg_a_rda
dev.off()

mean_val_tab[, .N, by = c('Lowland_class', 'ploidy')]

mean_val_tab[, Low_cyt_class := paste(mean_val_tab$Lowland_class, 
  mean_val_tab$ploidy, sep = '_')]

mean_val_tab[, Low_cyt_class := factor(mean_val_tab$Low_cyt_class, 
  levels = c('Low_4X', 'Low_8X', 'Moderate_4X', 'Moderate_8X', 'High_8X'))]

# Class x ploidy
gg_a_LowPloidy_box <- ggplot(mean_val_tab, aes(x = Low_cyt_class, y = MWvA_rda,
    fill = Lowland_class)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Lowland Allele Class') +
  ylab('Atlantic Alleles Score') +
  ylim(0, 0.3)


atl_boxplot_LowPloidy_out <- paste(out_dir, 
  'Atl_allele_score_LowPloidy_box.pdf', sep = '')

pdf(atl_boxplot_LowPloidy_out, height = 5, width = 6)
gg_a_LowPloidy_box
dev.off()

# Only ploidy
gg_atl_ploidy_box <- ggplot(mean_val_tab, aes(x = ploidy, y = MWvA_rda,
    fill = ploidy)) +
  geom_boxplot() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Cytotype') +
  ylab('Atlantic Alleles Score') +
  ylim(0, 0.3)

atl_boxplot_ploidy_out <- paste(out_dir, 
  'Atl_allele_score_ploidy_box.pdf', sep = '')

pdf(atl_boxplot_ploidy_out, height = 5, width = 4)
gg_atl_ploidy_box
dev.off()

### PCA of Midwest samples colored by introgression class

for(i in seq(nrow(mean_val_tab))){
  tmp_info_ind <- which(samp_info$samp_name == mean_val_tab$samp_name[i])
  samp_info[tmp_info_ind, Lowland_class := mean_val_tab$Lowland_class[i]]
}

samp_info[, Lowland_class := factor(samp_info$Lowland_class, 
  levels = c('Low', 'Moderate', 'High'))]

gg_pca_mw <- ggplot(samp_info[which(is.na(samp_info$Lowland_class) == F)], 
    aes(x = MW_pc01_raw, y = MW_pc02_raw, color = Lowland_class, 
    shape = ploidy)) +
  geom_point() +
  xlab('PC1') +
  ylab('PC2') +
  ggtitle('Midwest Genepool PC1 vs PC2')

mw_pca_LowlandClass_out <- paste(out_dir,
  'Midwest_PC1vPC2_LowlandClass_dotplot.pdf', sep = '')

pdf(mw_pca_LowlandClass_out, height = 5, width = 6)
gg_pca_mw
dev.off()


### Boxplots of Gulf vs Atlantic score

mean_rda_dif_1 <- ((mean_val_tab$MWvG_rda - mean_val_tab$MWvA_rda)/
  mean_val_tab$MWvA_rda)

mean_val_tab[, rda_dif_1 := mean_rda_dif_1]

# add groups for this

mean_val_tab[Low_cyt_class == 'Low_4X', EvW_class := 'Low_4X']
mean_val_tab[Low_cyt_class == 'Moderate_4X', EvW_class := 'Moderate_4X']
mean_val_tab[which(mean_val_tab$Low_cyt_class == 'Moderate_8X' & mean_val_tab$subgrp_v2 == 'MW_01'), EvW_class := 'MW8X-West']
mean_val_tab[Low_cyt_class == 'High_8X', EvW_class := 'MW8X-West-High']
mean_val_tab[which(mean_val_tab$Low_cyt_class == 'Moderate_8X' & mean_val_tab$subgrp_v2 == 'MW_03'), EvW_class := 'MW8X-East']

mean_val_tab[, EvW_class := factor(mean_val_tab$EvW_class, 
  levels = c('Low_4X', 'Moderate_4X', 'MW8X-East', 'MW8X-West', 
  'MW8X-West-High'))]

mean_val_tab_1 <- mean_val_tab[which(is.na(mean_val_tab$EvW_class) == F)]

mean_val_tab_1[, EvW_class := factor(EvW_class, 
  levels = c('Low_4X', 'Moderate_4X', 'MW8X-East', 'MW8X-West', 
  'MW8X-West-High'))]

gg_EvW_box <- ggplot(mean_val_tab[which(is.na(mean_val_tab$EvW_class) == F)], 
    aes(x = EvW_class, y = rda_dif_1, 
    fill = EvW_class)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Sample Group') +
  ylab('Atlantic Allele vs\nGulf Allele Signal')

EvW_box_1_out <- paste(out_dir, 'EastVsWest_Allele_boxplot.pdf', sep = '')

pdf(EvW_box_1_out, height = 5, width = 5)
gg_EvW_box
dev.off()

gg_EvW_box_2 <- ggplot(mean_val_tab[which(is.na(mean_val_tab$EvW_class) == F & 
    mean_val_tab$EvW_class != 'MW8X-West-High')],
    aes(x = EvW_class, y = rda_dif_1,
    fill = EvW_class)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab('Sample Group') +
  ylab('Atlantic Allele vs\nGulf Allele Signal')

EvW_box_2_out <- paste(out_dir, 'EastVsWest_Allele_boxplot_noHigh.pdf', 
  sep = '')

pdf(EvW_box_2_out, height = 5, width = 4)
gg_EvW_box_2
dev.off()










