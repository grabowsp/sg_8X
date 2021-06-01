# Assign genotypes and plot introgression statistics  using windows

### LOAD MODULES ###
# bash
# source activate R_analysis

args = commandArgs(trailingOnly=T)

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)
library(gridExtra)

func_dir <- paste('/home/grabowsky/tools/workflows/sg_8X/introgression/', 
  'comp_introgression_files/', sep = '')
introg_indiv_func_file <- paste(func_dir, 'comp_introg_indiv_functions.r', 
  sep = '')
source(introg_indiv_func_file)

### INPUT DATA ###
vcf_file <- args[1]
#vcf_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
#  'introgression_v3/', 'grp2_MWvGULF.vcf', sep = '')
vcf_in <- read.table(vcf_file, header = F, stringsAsFactors = F)

vcf_head_tmp <- system(paste('grep CHR ', vcf_file, sep = ''), intern = T)
vcf_head_2 <- sub('#CHROM', 'CHROM', vcf_head_tmp)
vcf_head_3 <- unlist(strsplit(vcf_head_2, split = '\t'))

colnames(vcf_in) <- vcf_head_3
vcf_in <- data.table(vcf_in)

allele_state_file <- args[2]
#allele_state_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
#  'introgression_v3/', 'MW_v_GULF_allele_states.txt', sep = '')
allele_states <- fread(allele_state_file)

### SET VARIABLES ###
testgrp_name <- args[3]
#testgrp_name <- 'grp2'

pop1_name <- args[4]
# pop1 = main background group
#pop1_name <- 'MW'

pop2_name <- args[5]
# pop2 = cadidate introgression group
#pop2_name <- 'GULF'

test_window_size <- as.numeric(args[6])
# test_window_size = the number of SNPs in each window
#test_window_size <- 10

### SET OUPUTS ###

out_dir <- args[7]
#out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'
if(rev(unlist(strsplit(out_dir, split = '')))[1] != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

# results from window-based analysis
window_result_file_out <- paste(out_dir, testgrp_name, '_', pop2_name,
  '_into_', pop1_name, '_', test_window_size, 'SNP_window', '_results.rds',
  sep = '')

#window_result_file_out <-  paste('/home/f2p1/work/grabowsk/data/switchgrass/',
#  'introgression_v3/', 'grp2_MW_v_GULF_10SNP_window_results.rds', sep = '')

# table and RDS file for genotypes to be used for environmental association
#   and Redundancy
rda_geno_file_out <- paste(out_dir, testgrp_name, '_', pop2_name,
  '_into_', pop1_name, '.NPOP2alleles.txt', sep = '')
rda_geno_file_rds <- paste(out_dir, testgrp_name, '_', pop2_name,
  '_into_', pop1_name, '.NPOP2alleles.rds', sep = '')

# plot files
tot_adj_n_pop2_plot_file <- paste(out_dir, testgrp_name, '_', pop2_name, 
  '_into_', pop1_name, '_', test_window_size, 'SNP_window', '.adj_n_pop2.pdf',
  sep = '')

tot_score_plot_file <- paste(out_dir, testgrp_name, '_', pop2_name,
  '_into_', pop1_name, '_', test_window_size, 'SNP_window', '.pop2_score.pdf',
  sep = '')

###
# Genotype Categories

# cat_a = 1:1 = 3:3 = dark blue (HOM pop1_ONLY)
# cat_b = 2:2 = 4:4 = medium blue (HOM pop1_MAINLY)
# cat_c = 1:Z = 3:Y = dark blue, less saturation (HET pop1_ONLY:NO_INFO)
# cat_d = 2:Z = 4:Y = medium blue, less saturation (HET pop1_MAINLY:NO_INFO)

# cat_e = A:A = C:C = dark red (HOM pop2_ONLY)
# cat_f = B:B = D:D = medium red (HOM pop2_MAINLY)
# cat_g = A:Z = C:Y = dark red, less saturaion (HET pop2_ONLY:NO_INFO)
# cat_h = B:Z = D:Y = medium red, less saturation (HET pop2_MAINLY:NO_INFO)

# cat_i = 1:C = 3:A =  magenta (HET pop1_ONLY:pop2_ONLY)
# cat_j = 1:D = 3:B =  blue-er magenta (more blue) (HET pop1_ONLY:pop2_MAINLY)
# cat_k = 2:C = 4:A =  red-er magenta (more red) (HET pop1_MAINLY:pop2_ONLY)
# cat_l = 2:D = 4:B =  light magenta (HET pop1_MAINLY:pop2_MAINLY)

# cat_m = Z:Z = Y:Y = grey (HOM NO_INFO)
# cat_n = missing data

### Introgression weights for calculating "score"
pop1_only_weight <- -1
pop1_mainly_weight <- -0.5
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

### RDS genos
# generate Gulf/pop2-presence table for RDA analysis
#  0 = No pop2 alleles; 1 = 1 pop2 allele, 2 = 2 pop2 alleles
# cat_e and cat_f = 2; cat_g, cat_h cat_i, cat_j, cat_k, cat_l = 1;
#  everything else 0

rda_geno_list <- list()
rda_geno_list[['cat_e']] <- 2
rda_geno_list[['cat_f']] <- 2
rda_geno_list[['cat_g']] <- 1
rda_geno_list[['cat_h']] <- 1
rda_geno_list[['cat_i']] <- 1
rda_geno_list[['cat_j']] <- 1
rda_geno_list[['cat_k']] <- 1
rda_geno_list[['cat_l']] <- 1
rda_geno_list[c('cat_a', 'cat_b', 'cat_c', 'cat_d', 'cat_m', 'cat_n')] <- 0

#####################
# adjust vcf so SNPs are ordered and there is a cumulative position for
#   plotting
# vector of sample names
samp_vec <- colnames(vcf_in)[c(10:ncol(vcf_in))]
# find max position for each chromosome
max_chrom_pos <- c()
for(chrm in unique(vcf_in$CHROM)){
  max_pos <- max(vcf_in$POS[vcf_in$CHROM == chrm])
  max_chrom_pos[[chrm]] <- max_pos
}
max_chrom_pos <- unlist(max_chrom_pos)
# figure out the cumulative amount to add to each chromosome
cum_chrom_amount <- rep(0, times = length(max_chrom_pos))
for(i in c(2:length(cum_chrom_amount))){
#  cum_chrom_amount[i] <- cum_chrom_amount[(i-1)] + max_chrom_pos[(i-1)]
  cum_chrom_amount[i] <- sum(max_chrom_pos[1:(i-1)])+(100000*(i-1))
}
names(cum_chrom_amount) <- names(max_chrom_pos)

vcf_in[, POS_CUM := as.numeric(NA)]
for(chrm in unique(vcf_in$CHROM)){
  tmp_inds <- which(vcf_in$CHROM == chrm)
  tmp_new_pos <- vcf_in$POS[tmp_inds] + cum_chrom_amount[chrm]
  vcf_in[tmp_inds, POS_CUM := tmp_new_pos]
}

vcf_in <- vcf_in[order(vcf_in$POS_CUM)]

if(sum(vcf_in$POS != allele_states$POS) > 0){
  printe('Need to re-order ALLELE STATES table SNPs')
}

#####################

# tot_pop2_list <- readRDS(window_result_file_out)

tot_pop2_list <- list()
for(sv in samp_vec){
  tot_pop2_list[[sv]] <- list()
  # genotype state based on allele states
  tmp_state_vec <- gen_geno_state_vec(vcf = vcf_in, samp_name = sv,
    allele_states = allele_states)
  tot_pop2_list[[sv]][['geno_state_vec']] <- tmp_state_vec
  #
  # genotype categories
  geno_cat_vec <- gen_cat_state_vec(geno_state_vec = tmp_state_vec,
    pop1_name = pop1_name, pop2_name = pop2_name)
  tot_pop2_list[[sv]][['geno_cat_vec']] <- geno_cat_vec
  #
  # calculate genotype scores
  geno_score_vec <- rep(NA, times = length(geno_cat_vec))
  for(cn in unique(geno_cat_vec)){
    geno_score_vec[which(geno_cat_vec == cn)] <- intro_wt_vec[[cn]]
  }
  # generate RDS genotype
  rda_geno_vec <- rep(NA, times = length(geno_cat_vec))
    for(cn in unique(geno_cat_vec)){
    rda_geno_vec[which(geno_cat_vec == cn)] <- rda_geno_list[[cn]]
  }
  tot_pop2_list[[sv]][['rda_geno_vec']] <- rda_geno_vec
  # Generate non-overlapping window table 
  tmp_nonoverlap_ind <- list()
  for(chrm in unique(vcf_in$CHROM)){
    tmp_tab <- gen_pop2_presence_nonoverlap_window_tab(
      samp_state_vec = tmp_state_vec, vcf = vcf_in, chrom = chrm, 
      window_size = test_window_size)
    tmp_score_tab <- gen_pop2_score_nonoverlap_window_tab(
      samp_score_vec = geno_score_vec, vcf = vcf_in, chrom = chrm,
    window_size = test_window_size)  
    tmp_tab[, POP2_SCORE := tmp_score_tab$POP2_SCORE]
    tmp_nonoverlap_ind[[chrm]] <- tmp_tab
  }
  tmp_nonoverlap_1 <- rbindlist(tmp_nonoverlap_ind)
  tot_pop2_list[[sv]][['non_overlap_wind']] <- tmp_nonoverlap_1
  #
  # Generate overlapping window table
  tmp_overlap_ind <- list()
  for(chrm in unique(vcf_in$CHROM)){
    tmp_tab <- gen_pop2_presence_overlap_window_tab(
      samp_state_vec = tmp_state_vec, vcf = vcf_in, chrom = chrm, 
      window_size = test_window_size)
    tmp_score_tab <- gen_pop2_score_overlap_window_tab(
      samp_score_vec = geno_score_vec, vcf = vcf_in, chrom = chrm,
    window_size = test_window_size)  
    tmp_tab[, POP2_SCORE := tmp_score_tab$POP2_SCORE]
    tmp_overlap_ind[[chrm]] <- tmp_tab
  }
  tmp_overlap_1 <- rbindlist(tmp_overlap_ind)
  tot_pop2_list[[sv]][['overlap_wind']] <- tmp_overlap_1
}

saveRDS(tot_pop2_list, file = window_result_file_out)

### Generate table of "genotypes", showing the number of Pop2 alleles, that
###   can be used for RDA analysis
rda_mat <- matrix(
  unlist(lapply(tot_pop2_list, function(x) x[['rda_geno_vec']])), byrow = F,
  ncol = length(tot_pop2_list))

colnames(rda_mat) <- names(tot_pop2_list)

rda_tab_1 <- data.table(rda_mat)

tmp_sub_vcf <- vcf_in[, c(1:5)]

tmp_sub_allele_states <- allele_states[, c(3,4)]

rda_tab_2 <- cbind(cbind(tmp_sub_vcf, tmp_sub_allele_states), rda_mat)

fwrite(rda_tab_2, file = rda_geno_file_out, sep = '\t')
saveRDS(rda_tab_2, file = rda_geno_file_rds)

### Plot results

# set palettes for plotting
black_chrom_palette <- rep(c('black', 'gray30'), times = 9)
names(black_chrom_palette) <- unique(tot_pop2_list[[1]][[4]]$CHROM)
black_chr_pal <- scale_colour_manual(name = 'CHROM', 
  values = black_chrom_palette)

blue_chrom_palette <- rep(c('blue3', 'gray30'), times = 9)
names(blue_chrom_palette) <- unique(tot_pop2_list[[1]][[4]]$CHROM)
blue_chr_pal <- scale_colour_manual(name = 'CHROM', 
  values = blue_chrom_palette)

# set y-axis limits for plotting
max_tot_adj_n_pop2 <- 50
#max_tot_adj_n_pop2 <- max(unlist(lapply(tot_pop2_list, 
#  function(x) max(x[[4]]$ADJ_N_POP2))))

max_tot_score <- 12.5
#max_tot_score <- max(unlist(lapply(tot_pop2_list,
#  function(x) max(x[[4]]$POP2_SCORE))))

# plot of "Adjusted Score" for all samples
tot_adj_gglist <- list()
for(i in seq(length(tot_pop2_list))){
  test_tab <- tot_pop2_list[[i]][[4]]
  tot_adj_gglist[[i]] <- ggplot(test_tab, 
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    black_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
      color = 'gray50', linetype = 'dashed') +
    xlab('') + 
    ggtitle(paste(names(tot_pop2_list)[i], '; ADJ_N_POP2 score; ', pop2_name, 
      ' into ', pop1_name, ' background; ', test_window_size, ' SNP windows', 
       sep = ''))
}

pdf(tot_adj_n_pop2_plot_file,  width = 18, height = 2*length(tot_adj_gglist))
grid.arrange(grobs = tot_adj_gglist, ncol = 1)
dev.off()

# plot of "total score" for all samples

tot_score_gglist <- list()
for(i in seq(length(tot_pop2_list))){
  test_tab <- tot_pop2_list[[i]][[4]]
  tot_score_gglist[[i]] <- ggplot(test_tab,   
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    blue_chr_pal + 
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(tot_pop2_list)[i], '; POP2_SCORE; ', pop2_name,
      ' into ', pop1_name, ' background; ', test_window_size, ' SNP windows',
       sep = ''))
}

pdf(tot_score_plot_file, width = 18, height = 2*length(tot_adj_gglist))
grid.arrange(grobs = tot_score_gglist, ncol = 1)
dev.off()

quit(save = 'no')


