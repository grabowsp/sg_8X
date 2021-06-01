# Compare levels of inferred GULF and ATL introgression in samples and
#   partition by chromosome

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)
library(gridExtra)

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

mw_g_chr_list <- lapply(mw_g_gw_list, function(x) 
  tapply(x$rda_geno, x$CHR, mean))
mw_g_chr_mat <- matrix(unlist(mw_g_chr_list), nrow = length(mw_g_chr_list),
  byrow = T)
colnames(mw_g_chr_mat) <- names(mw_g_chr_list[[1]])
rownames(mw_g_chr_mat) <- names(mw_g_chr_list)

# Find mean rda value when removing highest and lowest 3 chromosomes
cent_mean_list <- list()
for(i in seq(nrow(mw_g_chr_mat))){
  tmp_samp <- rownames(mw_g_chr_mat)[i]
  tmp_ord <- colnames(mw_g_chr_mat)[order(mw_g_chr_mat[i, ])]
  keep_chrs <- tmp_ord[c(4:15)]
  tmp_grep_string <- paste(keep_chrs, collapse = '|') 
  tmp_mean_inds <- grep(tmp_grep_string, mw_g_gw_list[[tmp_samp]]$CHR)
  tmp_cent_mean <- mean(mw_g_gw_list[[tmp_samp]][tmp_mean_inds, rda_geno])
  cent_mean_list[[tmp_samp]] <- tmp_cent_mean
}

cent_mean_vec <- unlist(cent_mean_list)
summary(cent_mean_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05051 0.06269 0.10983 0.11013 0.13339 0.32665 

tot_mean_vec <- unlist(lapply(mw_g_gw_list, function(x) mean(x[, rda_geno])))
summary(tot_mean_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05111 0.06232 0.11103 0.11040 0.13223 0.33200

# divide chromosomal mean values by "central mean"
mw_g_chr_mat_2 <- mw_g_chr_mat
for(j in seq(nrow(mw_g_chr_mat_2))){
  mw_g_chr_mat_2[j, ] <- mw_g_chr_mat_2[j, ] / cent_mean_vec[j]
}

# Generate data table with Chr values

mw_g_chr_tab <- data.table(samp_name = rownames(mw_g_chr_mat), mw_g_chr_mat_2,
  cent_mean = cent_mean_vec)
mw_g_chr_tab[, subgrp := as.character(NA)]
for(i in seq(nrow(mw_g_chr_tab))){
  tmp_info_ind <- which(samp_info$samp_name == mw_g_chr_tab$samp_name[i])
  mw_g_chr_tab[i, subgrp := samp_info$sub_grp[tmp_info_ind]]
}

# split off MW_01_hi
mean_mw_v_g_rda <- unlist(lapply(mw_g_gw_list, function(x) mean(x$rda_geno)))

mw_g_chr_tab[intersect(which(mw_g_chr_tab$subgrp == 'MW_01'),
  which(mean_mw_v_g_rda > 0.21)), subgrp := 'MW_01_hi']

mw_g_chr_tab[, sub_2 := '4X']
mw_g_chr_tab[grep('MW_01|MW_02|MW_03|MW_04|MW_05', mw_g_chr_tab$subgrp), 
  sub_2 := '8X']
mw_g_chr_tab[which(mw_g_chr_tab$subgrp == 'MW_01_hi'), sub_2 := '8X_hi']

mw_g_chr_tab[, sub_3 := mw_g_chr_tab$sub_2]
mw_g_chr_tab[grep('MW_06|MW_07', mw_g_chr_tab$subgrp), sub_3 := '4X_intro']

######
min_y_val <- min(mw_g_chr_mat_2)
max_y_val <- 2.5

#####
chr_names <- colnames(mw_g_chr_mat_2)
#####

# boxplot of Chr values
gg_chr01k_box <- ggplot(mw_g_chr_tab, aes(x = subgrp, y = Chr01K)) +
  geom_boxplot() +
  ylim(min_y_val, max_y_val) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = 'Mean Gulf alleles in MW subgroups vs\nGenome-wide average')

chr01k_box_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_Chr01K_Gulf_rda.pdf', sep = '')

pdf(chr01k_box_out, height = 5, width = 5)
gg_chr01k_box
dev.off()

chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(mw_g_chr_tab, aes_string(x = 'subgrp', y = cn)) +
  geom_boxplot() +
  ylim(min_y_val, max_y_val) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = paste('Mean Gulf alleles in ' , cn, ' vs Genome-wide average,', 
    '\nin MW subgroups', sep = ''))
  chr_box_list[[cn]] <- tmp_gg
}

allchr_box_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Gulf_rda.pdf', 
  sep = '')

pdf(allchr_box_out, width = 5*length(chr_box_list), height = 5)
grid.arrange(grobs = chr_box_list, nrow = 1)
dev.off()

chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(mw_g_chr_tab, aes_string(x = 'sub_2', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' Mean\nGulf alleles vs\nGenome-wide\naverage,', 
      sep = '')) +
    xlab('MW groups')
  chr_box_list_2[[cn]] <- tmp_gg
}

allchr_box_out_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Gulf_rda_v2.pdf',
  sep = '')

pdf(allchr_box_out_2, width = 2*length(chr_box_list_2), height = 5)
grid.arrange(grobs = chr_box_list_2, nrow = 1)
dev.off()

chr_box_list_3 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(mw_g_chr_tab, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' Mean\nGulf alleles vs\nGenome-wide\naverage,',
      sep = '')) +
    xlab('MW groups')
  chr_box_list_3[[cn]] <- tmp_gg
}

allchr_box_out_3 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Gulf_rda_v3.pdf',
  sep = '')

pdf(allchr_box_out_3, width = 2*length(chr_box_list_3), height = 5)
grid.arrange(grobs = chr_box_list_3, nrow = 1)
dev.off()

allchr_box_out_3_1 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Gulf_rda_v3.1.pdf',
  sep = '')

pdf(allchr_box_out_3_1, width = 2*5, height = 3.5*4)
grid.arrange(grobs = chr_box_list_3, nrow = 4)
dev.off()

### show patterns by group (all chromosomes next to eachother)

test <- as.vector(mw_g_chr_mat_2)
sub2_4X_inds <- which(mw_g_chr_tab$sub_2 == '4X')

sub2_4X_geno_vec <- as.vector(mw_g_chr_mat_2[sub2_4X_inds, ])
sub2_4X_chr_vec <- rep(chr_names, each = length(sub2_4X_inds))
sub2_4X_tab <- data.table(avg_pop2_allele = sub2_4X_geno_vec,
  CHROM = sub2_4X_chr_vec)

test_pop_chr_gg <- ggplot(sub2_4X_tab, 
    aes_string(x = 'CHROM', y = 'avg_pop2_allele')) + 
  geom_boxplot() +
  ylim(min_y_val, max_y_val) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = paste('Mean Gulf alleles vs\nGenome-wide average in MW-4X', 
    sep = ''))

chr_group_test_pdf <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_MW4X_Chromosomes_Gulf_rda.pdf',
  sep = '')

pdf(chr_group_test_pdf, width = 8, height = 5)
test_pop_chr_gg
dev.off()

##

sub2_grps <- unique(mw_g_chr_tab$sub_2)

sub2_chr_gg_list <- list()
for(s2 in sub2_grps){
  tmp_inds <- which(mw_g_chr_tab$sub_2 == s2)
  tmp_geno_vec <- as.vector(mw_g_chr_mat_2[tmp_inds, ])
  tmp_chr_vec <- rep(chr_names, each = length(tmp_inds))
  tmp_tab <- data.table(avg_pop2_allele = tmp_geno_vec,
    CHROM = tmp_chr_vec)
  tmp_gg <- ggplot(tmp_tab, aes_string(x = 'CHROM', y = 'avg_pop2_allele')) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste('Mean Gulf alleles by CHROM ', 
      '\nvs Genome-wide average for ', s2, sep = ''))
  sub2_chr_gg_list[[s2]] <- tmp_gg
}

sub2_group_chr_pdf <-  paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_Sub2Groups_Chromosomes_Gulf_rda.pdf',
  sep = '')

pdf(sub2_group_chr_pdf, width = 8, height = 5*length(sub2_chr_gg_list))
grid.arrange(grobs = sub2_chr_gg_list, ncol = 1)
dev.off()

##
sub3_grps <- sort(unique(mw_g_chr_tab$sub_3))

sub3_chr_gg_list <- list()
for(s3 in sub3_grps){
  tmp_inds <- which(mw_g_chr_tab$sub_3 == s3)
  tmp_geno_vec <- as.vector(mw_g_chr_mat_2[tmp_inds, ])
  tmp_chr_vec <- rep(chr_names, each = length(tmp_inds))
  tmp_tab <- data.table(avg_pop2_allele = tmp_geno_vec,
    CHROM = tmp_chr_vec)
  tmp_gg <- ggplot(tmp_tab, aes_string(x = 'CHROM', y = 'avg_pop2_allele')) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste('Mean Gulf alleles by CHROM ',
      '\nvs Genome-wide average for ', s3, sep = ''))
  sub3_chr_gg_list[[s3]] <- tmp_gg
}

sub3_group_chr_pdf <-  paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_Sub3Groups_Chromosomes_Gulf_rda.pdf',
  sep = '')

pdf(sub3_group_chr_pdf, width = 8, height = 5*length(sub3_chr_gg_list))
grid.arrange(grobs = sub3_chr_gg_list, ncol = 1)
dev.off()

####
####
####

# MW vs ATL results
mw_a_gw_list <- list()

for(mwa in names(mw_a_res_list)){
  tmp_tab <- data.table(mw_a_pos,
    geno_state = mw_a_res_list[[mwa]]$geno_state_vec,
    geno_cat = mw_a_res_list[[mwa]]$geno_cat_vec,
    rda_geno = mw_a_res_list[[mwa]]$rda_geno_vec)
  tmp_score_vec <- rep(NA, times = nrow(tmp_tab))
  for(cn in unique(tmp_tab$geno_cat)){
    tmp_score_vec[which(tmp_tab$geno_cat == cn)] <- intro_wt_vec[[cn]]
  }
  tmp_tab[ , geno_score := tmp_score_vec]
  mw_a_gw_list[[mwa]] <- tmp_tab
}

mw_a_chr_list <- lapply(mw_a_gw_list, function(x)
  tapply(x$rda_geno, x$CHR, mean))
mw_a_chr_mat <- matrix(unlist(mw_a_chr_list), nrow = length(mw_a_chr_list),
  byrow = T)
colnames(mw_a_chr_mat) <- names(mw_a_chr_list[[1]])
rownames(mw_a_chr_mat) <- names(mw_a_chr_list)

# Find mean rda value when removing highest and lowest 3 chromosomes
atl_cent_mean_list <- list()
for(i in seq(nrow(mw_a_chr_mat))){
  tmp_samp <- rownames(mw_a_chr_mat)[i]
  tmp_ord <- colnames(mw_a_chr_mat)[order(mw_a_chr_mat[i, ])]
  keep_chrs <- tmp_ord[c(4:15)]
  tmp_grep_string <- paste(keep_chrs, collapse = '|')
  tmp_mean_inds <- grep(tmp_grep_string, mw_a_gw_list[[tmp_samp]]$CHR)
  tmp_cent_mean <- mean(mw_a_gw_list[[tmp_samp]][tmp_mean_inds, rda_geno])
  atl_cent_mean_list[[tmp_samp]] <- tmp_cent_mean
}

atl_cent_mean_vec <- unlist(atl_cent_mean_list)
summary(atl_cent_mean_vec)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.05222 0.06363 0.11357 0.10943 0.13607 0.26312 

atl_tot_mean_vec <- unlist(lapply(mw_a_gw_list, 
  function(x) mean(x[, rda_geno])))
summary(atl_tot_mean_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.05258 0.06375 0.11536 0.10997 0.13630 0.26121 

# divide chromosomal mean values by "central mean"
mw_a_chr_mat_2 <- mw_a_chr_mat
for(j in seq(nrow(mw_a_chr_mat_2))){
  mw_a_chr_mat_2[j, ] <- mw_a_chr_mat_2[j, ] / atl_cent_mean_vec[j]
}

# Generate data table with Chr values

mw_a_chr_tab <- data.table(samp_name = rownames(mw_a_chr_mat), mw_a_chr_mat_2,
  cent_mean = atl_cent_mean_vec)
mw_a_chr_tab[, subgrp := as.character(NA)]
for(i in seq(nrow(mw_a_chr_tab))){
  tmp_info_ind <- which(samp_info$samp_name == mw_a_chr_tab$samp_name[i])
  mw_a_chr_tab[i, subgrp := samp_info$sub_grp[tmp_info_ind]]
}

# split off MW_01_hi
# use Gulf values as cutoff
mean_mw_v_g_rda <- unlist(lapply(mw_g_gw_list, function(x) mean(x$rda_geno)))

mw_a_chr_tab[intersect(which(mw_a_chr_tab$subgrp == 'MW_01'),
  which(mean_mw_v_g_rda > 0.21)), subgrp := 'MW_01_hi']

mw_a_chr_tab[, sub_2 := '4X']
mw_a_chr_tab[grep('MW_01|MW_02|MW_03|MW_04|MW_05', mw_a_chr_tab$subgrp),
  sub_2 := '8X']
mw_a_chr_tab[which(mw_a_chr_tab$subgrp == 'MW_01_hi'), sub_2 := '8X_hi']

mw_a_chr_tab[, sub_3 := mw_a_chr_tab$sub_2]
mw_a_chr_tab[grep('MW_06|MW_07', mw_a_chr_tab$subgrp), sub_3 := '4X_intro']

######
min_y_val <- min(mw_a_chr_mat_2)
max_y_val <- 2.5

#####
chr_names <- colnames(mw_a_chr_mat_2)
#####

atl_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(mw_a_chr_tab, aes_string(x = 'subgrp', y = cn)) +
  geom_boxplot() +
  ylim(min_y_val, max_y_val) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = paste('Mean Atlantic alleles in ' , cn, 
    ' vs Genome-wide average\nin MW subgroups', sep = ''))
  atl_chr_box_list[[cn]] <- tmp_gg
}

atl_allchr_box_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Atlantic_rda.pdf', 
  sep = '')

pdf(atl_allchr_box_out, width = 5*length(atl_chr_box_list), height = 5)
grid.arrange(grobs = atl_chr_box_list, nrow = 1)
dev.off()

atl_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(mw_a_chr_tab, aes_string(x = 'sub_2', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' Mean\nAtlantic alleles vs\nGenome-wide\naverage,', 
      sep = '')) +
    xlab('MW groups')
  atl_chr_box_list_2[[cn]] <- tmp_gg
}

atl_allchr_box_out_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Atlantic_rda_v2.pdf',
  sep = '')

pdf(atl_allchr_box_out_2, width = 2*length(atl_chr_box_list_2), height = 5)
grid.arrange(grobs = atl_chr_box_list_2, nrow = 1)
dev.off()

atl_chr_box_list_3 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(mw_a_chr_tab, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' Mean\nAtlantic alleles vs\nGenome-wide\naverage,',
      sep = '')) +
    xlab('MW groups')
  atl_chr_box_list_3[[cn]] <- tmp_gg
}

atl_allchr_box_out_3 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_allChromosomes_Atlantic_rda_v3.pdf',
  sep = '')

pdf(atl_allchr_box_out_3, width = 2*length(atl_chr_box_list_3), height = 5)
grid.arrange(grobs = atl_chr_box_list_3, nrow = 1)
dev.off()

# CONTINUE FROM HERE
# Next: 1) ATL results - by group; 2) PP-friendly plots

### show patterns by group (all chromosomes next to eachother)

atl_sub2_4X_inds <- which(mw_a_chr_tab$sub_2 == '4X')

atl_sub2_4X_geno_vec <- as.vector(mw_a_chr_mat_2[atl_sub2_4X_inds, ])
atl_sub2_4X_chr_vec <- rep(chr_names, each = length(atl_sub2_4X_inds))
atl_sub2_4X_tab <- data.table(avg_pop2_allele = atl_sub2_4X_geno_vec,
  CHROM = atl_sub2_4X_chr_vec)

atl_test_pop_chr_gg <- ggplot(atl_sub2_4X_tab, 
    aes_string(x = 'CHROM', y = 'avg_pop2_allele')) + 
  geom_boxplot() +
  ylim(min_y_val, max_y_val) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = paste('Mean Atlantic alleles vs\nGenome-wide average in MW-4X', 
    sep = ''))

atl_chr_group_test_pdf <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MWall_MW4X_Chromosomes_Atlantic_rda.pdf',
  sep = '')

pdf(atl_chr_group_test_pdf, width = 8, height = 5)
atl_test_pop_chr_gg
dev.off()

##

sub2_grps <- unique(mw_a_chr_tab$sub_2)

atl_sub2_chr_gg_list <- list()
for(s2 in sub2_grps){
  tmp_inds <- which(mw_a_chr_tab$sub_2 == s2)
  tmp_geno_vec <- as.vector(mw_a_chr_mat_2[tmp_inds, ])
  tmp_chr_vec <- rep(chr_names, each = length(tmp_inds))
  tmp_tab <- data.table(avg_pop2_allele = tmp_geno_vec,
    CHROM = tmp_chr_vec)
  tmp_gg <- ggplot(tmp_tab, aes_string(x = 'CHROM', y = 'avg_pop2_allele')) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste('Mean Atlantic alleles by CHROM ', 
      '\nvs Genome-wide average for ', s2, sep = ''))
  atl_sub2_chr_gg_list[[s2]] <- tmp_gg
}

atl_sub2_group_chr_pdf <-  paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 
  'MWall_Sub2Groups_Chromosomes_Atlantic_rda.pdf',
  sep = '')

pdf(atl_sub2_group_chr_pdf, width = 8, height = 5*length(atl_sub2_chr_gg_list))
grid.arrange(grobs = atl_sub2_chr_gg_list, ncol = 1)
dev.off()

##
sub3_grps <- sort(unique(mw_a_chr_tab$sub_3))

atl_sub3_chr_gg_list <- list()
for(s3 in sub3_grps){
  tmp_inds <- which(mw_a_chr_tab$sub_3 == s3)
  tmp_geno_vec <- as.vector(mw_a_chr_mat_2[tmp_inds, ])
  tmp_chr_vec <- rep(chr_names, each = length(tmp_inds))
  tmp_tab <- data.table(avg_pop2_allele = tmp_geno_vec,
    CHROM = tmp_chr_vec)
  tmp_gg <- ggplot(tmp_tab, aes_string(x = 'CHROM', y = 'avg_pop2_allele')) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste('Mean Atlantic alleles by CHROM ',
      '\nvs Genome-wide average for ', s3, sep = ''))
  atl_sub3_chr_gg_list[[s3]] <- tmp_gg
}

atl_sub3_group_chr_pdf <-  paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 
  'MWall_Sub3Groups_Chromosomes_Atlantic_rda.pdf',
  sep = '')

pdf(atl_sub3_group_chr_pdf, width = 8, height = 5*length(atl_sub3_chr_gg_list))
grid.arrange(grobs = atl_sub3_chr_gg_list, ncol = 1)
dev.off()

#





# Explore results
low_inds <- grep('MW_10|MW_12|MW_13|MW_14', mw_g_chr_tab$subgrp)


apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_01'), ], 2, mean)

apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_01_hi'), ], 2, mean)

apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_03'), ], 2, mean)

apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_06'), ], 2, mean)

apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_07'), ], 2, mean)

apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_10'), ], 2, mean)
apply(mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_14'), ], 2, mean)

apply(mw_g_chr_mat_2[low_inds, ], 2, mean)

mw_g_chr_mat_2[which(mw_g_chr_tab$subgrp == 'MW_11'),]

mean(which(mw_g_chr_tab))


# Still need to look at ATL results





