# Analysis of Chr08N because of odd introgression patterns

# Plan:
# 1) tally rda genotypes by chromosome
#  - does Chr08N have higher fraction of homozygous Pop2 genotypes?
# 2) Look at overall homozygosity/heterozygosity across chromosomes
#  - is there an overall difference in heterozogysity that could explain
#     the patterns

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

# VCFtools HET/HOM results
mw_chr_HOM_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MWall_Chr_HOM_vals.txt'
mw_chr_HOM <- fread(mw_chr_HOM_file)

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

table(mw_g_gw_list[[1]]$rda_geno[mw_g_gw_list[[1]]$CHR == 'Chr01K'])

allsamp_rda_tmp_mat <- matrix(as.numeric(NA), ncol = length(mw_g_gw_list), 
  nrow = nrow(mw_g_gw_list[[1]]))
for(i in seq(length(mw_g_gw_list))){
  allsamp_rda_tmp_mat[,i] <- mw_g_gw_list[[i]]$rda_geno
}
colnames(allsamp_rda_tmp_mat) <- names(mw_g_gw_list)

allsamp_rda_tab <- data.table(
  CHR = mw_g_gw_list[[1]]$CHR, 
  POS = mw_g_gw_list[[1]]$POS,
  allsamp_rda_tmp_mat)

# Split off MW_01_hi
mean_mw_v_g_rda <- unlist(lapply(mw_g_gw_list, function(x) mean(x$rda_geno)))
mw_01_hi_names <- names(mean_mw_v_g_rda)[which(mean_mw_v_g_rda > 0.21)]

# tabulate the number of 1's and 2's and find ratio
rda_2_list <- list()
for(sn in names(mw_g_gw_list)){
  tmp_table_list <- tapply(allsamp_rda_tab[, get(sn)], allsamp_rda_tab$CHR,
    table) 
  tmp_ratio <- unlist(lapply(tmp_table_list, function(x) x['2']/x['1']))
  rda_2_list[[sn]] <- tmp_ratio
}

rda_2_mat <- matrix(unlist(rda_2_list), nrow = length(rda_2_list), byrow = T)
colnames(rda_2_mat) <- unique(mw_g_gw_list[[1]]$CHR)

rda_2_tab <- data.table(samp_name = names(rda_2_list), rda_2_mat)

# Add subgroup info
rda_2_tab[, subgrp := as.character(NA)]
for(i in seq(nrow(rda_2_tab))){
  tmp_info_ind <- which(samp_info$samp_name == rda_2_tab$samp_name[i])
  rda_2_tab[i, subgrp := samp_info$sub_grp[tmp_info_ind]]
}
rda_2_tab[which(rda_2_tab$samp_name %in% mw_01_hi_names), subgrp := 'MW_01_hi']

# group subgroups together
rda_2_tab[, sub_2 := '4X']
rda_2_tab[grep('MW_01|MW_02|MW_03|MW_04|MW_05', rda_2_tab$subgrp),
  sub_2 := '8X']
rda_2_tab[which(rda_2_tab$subgrp == 'MW_01_hi'), sub_2 := '8X_hi']

rda_2_tab[, sub_3 := rda_2_tab$sub_2]
rda_2_tab[grep('MW_06|MW_07', rda_2_tab$subgrp), sub_3 := '4X_intro']

# Plot 2 to 1 ratios
min_y_val <- 0
max_y_val <- 4.05

chr_names <- colnames(rda_2_mat)

rda_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(rda_2_tab, aes_string(x = 'subgrp', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste('Ratio of 2- to 1-Pop2 Allele genotypes', sep = ''))
  rda_chr_box_list[[cn]] <- tmp_gg
}

rda_box_out_allsubgrps <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MW_Gulf_Pop2Allele_ratio_boxplot_v1.pdf',
  sep = '')

pdf(rda_box_out_allsubgrps, width = 5*4, height = 5*5)
grid.arrange(grobs = rda_chr_box_list, nrow = 5)
dev.off()

rda_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(rda_2_tab, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' HOM:HET Pop2-Allele\nGenotype Ratios', sep = ''))
  rda_chr_box_list_2[[cn]] <- tmp_gg
}

rda_box_out_allsubgrps_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 'MW_Gulf_Pop2Allele_ratio_boxplot_v2.pdf',
  sep = '')

pdf(rda_box_out_allsubgrps_2, width = 3*4, height = 4*5)
grid.arrange(grobs = rda_chr_box_list_2, nrow = 5)
dev.off()

cent_ratio_list <- list()
for(i in seq(nrow(rda_2_mat))){
  tmp_samp <- names(mw_g_gw_list)[i]
  tmp_ord <- colnames(rda_2_mat)[order(rda_2_mat[i, ])]
  keep_chrs <- tmp_ord[c(4:15)]
  tmp_grep_string <- paste(keep_chrs, collapse = '|')
  tmp_cent_inds <- grep(tmp_grep_string, mw_g_gw_list[[tmp_samp]]$CHR)
  tmp_cent_table <- table(mw_g_gw_list[[tmp_samp]][tmp_cent_inds, rda_geno])
  tmp_cent_ratio <- tmp_cent_table['2']/tmp_cent_table['1']
  cent_ratio_list[[tmp_samp]] <- tmp_cent_ratio
}

cent_ratio_vec <- unlist(cent_ratio_list)

rda_2_mat_2 <- rda_2_mat
for(j in seq(nrow(rda_2_mat_2))){
  rda_2_mat_2[j, ] <- rda_2_mat_2[j, ] / cent_ratio_vec[j]
}

rda_2_tab_2 <- data.table(
  samp_name = rda_2_tab$samp_name,
  subgrp = rda_2_tab$subgrp,
  sub_2 = rda_2_tab$sub_2,
  sub_3 = rda_2_tab$sub_3,
  rda_2_mat_2)

####
min_y_val <- 0
max_y_val <- 5

chr_names <- colnames(rda_2_mat)
###

trim_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(rda_2_tab_2, aes_string(x = 'subgrp', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste('HOM:HET Pop2-allele\nRatios Divided by', 
      '\nTrimmed GW Mean', sep = ''))
  trim_chr_box_list[[cn]] <- tmp_gg
}

trim_box_out_allsubgrps <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 
  'MW_Gulf_Pop2Allele_ratio_trimMean_boxplot.pdf',
  sep = '')

pdf(trim_box_out_allsubgrps, width = 5*4, height = 4*5)
grid.arrange(grobs = trim_chr_box_list, nrow = 5)
dev.off()

trim_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(rda_2_tab_2, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' HOM:HET Pop2-\nallele Ratios Divided by', 
      '\nTrimmed GW Mean', sep = ''))
  trim_chr_box_list_2[[cn]] <- tmp_gg
}

trim_box_out_allsubgrps_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/', 
  'MW_Gulf_Pop2Allele_ratio_trimMean_boxplot_v2.pdf',
  sep = '')

pdf(trim_box_out_allsubgrps_2, width = 2.5*4, height = 4*5)
grid.arrange(grobs = trim_chr_box_list_2, nrow = 5)
dev.off()

# Calculate overall heterozygosity at each chromosome and look for patterns

hom_mat <- as.matrix(mw_chr_HOM[, c(2:19)])

cent_hom_list <- list()
for(i in seq(nrow(hom_mat))){
  tmp_samp <- mw_chr_HOM$samp_name[i]
  tmp_ord <- colnames(hom_mat)[order(hom_mat[i, ])]
  keep_chrs <- tmp_ord[c(4:15)]
  tmp_cent_mean <- mean(hom_mat[i, keep_chrs])
  cent_hom_list[[tmp_samp]] <- tmp_cent_mean
}

cent_hom_vec <- unlist(cent_hom_list)

hom_mat_2 <- hom_mat
for(j in seq(nrow(hom_mat_2))){
  hom_mat_2[j, ] <- hom_mat_2[j, ] / cent_hom_vec[j]
}

hom_tab <- data.table(
  samp_name = mw_chr_HOM$samp_name,
  subgrp = as.character(NA),
  sub_2 = as.character(NA),
  sub_3 = as.character(NA), 
  hom_mat_2
)

for(j in seq(nrow(hom_tab))){
  tmp_ind <- which(rda_2_tab_2$samp_name == hom_tab$samp_name[j])
  hom_tab[j, subgrp := rda_2_tab_2$subgrp[tmp_ind]]
  hom_tab[j, sub_2 := rda_2_tab_2$sub_2[tmp_ind]]
  hom_tab[j, sub_3 := rda_2_tab_2$sub_3[tmp_ind]]
}

####
min_y_val <- 0.9
max_y_val <- 1.1

chr_names <- colnames(hom_mat)
###

hom_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(hom_tab, aes_string(x = 'subgrp', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' Homozygosity divided by',
      '\nTrimmed GW Mean', sep = ''))
  hom_chr_box_list[[cn]] <- tmp_gg
}

hom_box_out_allsubgrps <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_HOM_trimMean_boxplot_v1.pdf',
  sep = '')

pdf(hom_box_out_allsubgrps, width = 5*4, height = 4*5)
grid.arrange(grobs = hom_chr_box_list, nrow = 5)
dev.off()

hom_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(hom_tab, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' Homozygosity divided by',
      '\nTrimmed GW Mean', sep = ''))
  hom_chr_box_list_2[[cn]] <- tmp_gg
}

hom_box_out_allsubgrps_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_HOM_trimMean_boxplot_v2.pdf',
  sep = '')

pdf(hom_box_out_allsubgrps_2, width = 2.5*4, height = 4*5)
grid.arrange(grobs = hom_chr_box_list_2, nrow = 5)
dev.off()

# CONTINUE WITH HET Boxplots
het_mat <- 1-hom_mat
cent_het_vec <- 1-cent_hom_vec

et_mat_2 <- het_mat
for(j in seq(nrow(het_mat_2))){
  het_mat_2[j, ] <- het_mat_2[j, ] / cent_het_vec[j]
}

het_tab <- data.table(
  samp_name = mw_chr_HOM$samp_name,
  subgrp = hom_tab$subgrp,
  sub_2 = hom_tab$sub_2,
  sub_3 = hom_tab$sub_3,
  het_mat_2
)

####
min_y_val <- 0.4
max_y_val <- 1.4

chr_names <- colnames(hom_mat)

het_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(het_tab, aes_string(x = 'subgrp', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' HET divided by',
      '\nTrimmed GW Mean', sep = ''))
  het_chr_box_list[[cn]] <- tmp_gg
}

het_box_out_allsubgrps <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_HET_trimMean_boxplot_v1.pdf',
  sep = '')

pdf(het_box_out_allsubgrps, width = 5*4, height = 4*5)
grid.arrange(grobs = het_chr_box_list, nrow = 5)
dev.off()

het_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(het_tab, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' HET\ndivided by',
      '\nTrimmed GW Mean', sep = ''))
  het_chr_box_list_2[[cn]] <- tmp_gg
}

het_box_out_allsubgrps_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_HET_trimMean_boxplot_v2.pdf',
  sep = '')

pdf(het_box_out_allsubgrps_2, width = 2.5*4, height = 4*5)
grid.arrange(grobs = het_chr_box_list_2, nrow = 5)
dev.off()

########
# tabulate the number of Pop2 to non-Pop2 SNPs
pop2_list <- list()
for(sn in names(mw_g_gw_list)){
  tmp_table_list <- tapply(allsamp_rda_tab[, get(sn)], allsamp_rda_tab$CHR,
    table)
  tmp_ratio <- unlist(lapply(tmp_table_list, function(x) 
    (x['2']+x['1'])/sum(x)))
  pop2_list[[sn]] <- tmp_ratio
}

pop2_mat <- matrix(unlist(pop2_list), nrow = length(pop2_list), byrow = T)
colnames(pop2_mat) <- unique(mw_g_gw_list[[1]]$CHR)

pop2_tab <- data.table(
  samp_name = rda_2_tab$samp_name,
  subgrp = rda_2_tab$subgrp,
  sub_2 = rda_2_tab$sub_2,
  sub_3 = rda_2_tab$sub_3,
  pop2_mat)


####
min_y_val <- 0
max_y_val <- 0.4

chr_names <- colnames(pop2_mat)

pop2_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(pop2_tab, aes_string(x = 'subgrp', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' % SNPs\nwith Pop2 Alleles', sep = ''))
  pop2_chr_box_list[[cn]] <- tmp_gg
}

pop2_box_out_allsubgrps <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_per_Pop2_chr_boxplot_v1.pdf',
  sep = '')

pdf(pop2_box_out_allsubgrps, width = 5*4, height = 4*5)
grid.arrange(grobs = pop2_chr_box_list, nrow = 5)
dev.off()

pop2_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(pop2_tab, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' % SNPs\nwith Pop2 Alleles', sep = ''))
  pop2_chr_box_list_2[[cn]] <- tmp_gg
}

pop2_box_out_allsubgrps_2 <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_per_Pop2_chr_boxplot_v2.pdf',
  sep = '')

pdf(pop2_box_out_allsubgrps_2, width = 2.5*4, height = 4*5)
grid.arrange(grobs = pop2_chr_box_list_2, nrow = 5)
dev.off()

cent_pop2_list <- list()
for(i in seq(nrow(pop2_mat))){
  tmp_samp <- names(mw_g_gw_list)[i]
  tmp_ord <- colnames(pop2_mat)[order(pop2_mat[i, ])]
  keep_chrs <- tmp_ord[c(4:15)]
  tmp_grep_string <- paste(keep_chrs, collapse = '|')
  tmp_cent_inds <- grep(tmp_grep_string, mw_g_gw_list[[tmp_samp]]$CHR)
  tmp_cent_table <- table(mw_g_gw_list[[tmp_samp]][tmp_cent_inds, rda_geno])
  tmp_cent_per <- ((tmp_cent_table['2'] + tmp_cent_table['1'])/
    sum(tmp_cent_table))
  cent_pop2_list[[tmp_samp]] <- tmp_cent_per
}

cent_pop2_vec <- unlist(cent_pop2_list)

pop2_mat_2 <- pop2_mat
for(j in seq(nrow(pop2_mat_2))){
  pop2_mat_2[j, ] <- pop2_mat_2[j, ] / cent_pop2_vec[j]
}
colnames(pop2_mat_2) <- unique(mw_g_gw_list[[1]]$CHR)

pop2_tab_2 <- data.table(
  samp_name = pop2_tab$samp_name,
  subgrp = pop2_tab$subgrp,
  sub_2 = pop2_tab$sub_2,
  sub_3 = pop2_tab$sub_3,
  pop2_mat_2)

####
min_y_val <- 0
max_y_val <- 3

chr_names <- colnames(pop2_mat)

pop2trim_chr_box_list <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(pop2_tab_2, aes_string(x = 'subgrp', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' % SNPs\nwith Pop2 Alleles', 
      '\ndivided by GW mean', sep = ''))
  pop2trim_chr_box_list[[cn]] <- tmp_gg
}

pop2trim_box_out_allsubgrps <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_per_Pop2_trimMean_chr_boxplot_v1.pdf',
  sep = '')

pdf(pop2trim_box_out_allsubgrps, width = 5*4, height = 4*5)
grid.arrange(grobs = pop2trim_chr_box_list, nrow = 5)
dev.off()

pop2trim_chr_box_list_2 <- list()
for(cn in chr_names){
  tmp_gg <- ggplot(pop2_tab_2, aes_string(x = 'sub_3', y = cn)) +
    geom_boxplot() +
    ylim(min_y_val, max_y_val) +
    theme(axis.text.x = element_text(angle = 90)) +
    labs(title = paste(cn, ' % SNPs\nwith Pop2 Alleles',
      '\ndivided by GW mean', sep = ''))
  pop2trim_chr_box_list_2[[cn]] <- tmp_gg
}

pop2trim_box_out_allsubgrps_2 <- paste(
  '/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_per_Pop2_trimMean_chr_boxplot_v2.pdf',
  sep = '')

pdf(pop2trim_box_out_allsubgrps_2, width = 2.5*4, height = 4*5)
grid.arrange(grobs = pop2trim_chr_box_list_2, nrow = 5)
dev.off()


