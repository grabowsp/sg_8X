# Look at allele frequencies in test group and and look at movement towards 
#  pop2 freqs, the population suspected of being a source of introgression

# Steps
# calculate F_value
# Generate window values
# Select top windows using sd-based cutoff
# Find SNPs that are found X+ times (ex 5) in top windows
# Those SNPs/loci are candidate regions

# bash
# source activate R_analysis

args <- commandArgs(trailingOnly = T)

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)
library(gridExtra)

func_dir <- paste('/home/grabowsky/tools/workflows/sg_8X/introgression/',
  'comp_introgression_files/', sep = '')
introg_pop_func_file <- paste(func_dir, 'comp_introg_pop_functions.r',
  sep = '')
source(introg_pop_func_file)

### INPUT DATA ###
data_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'

subgrp_hw_file <- paste(data_dir, 'MW8X_all_MWvGULF.hwe', sep = '')
#subgrp_hw_file <- args[1]
subgrp_hw <- fread(subgrp_hw_file)

train_freq_file <- paste(data_dir, 'MW_v_GULF_ref_freq.txt', sep = '')
#train_freq_file <- args[2]
train_freq <- fread(train_freq_file)

allele_state_file <- paste(data_dir, 'MW_v_GULF_allele_states.txt', sep = '')
#allele_state_file <- args[3]
allele_states <- fread(allele_state_file)

### SET VARIABLES ###

# minimum percentage of subgroup samples have genotypes at a SNP
miss_cut <- 0.8

#window_size <- as.numeric(args[4])
window_size <- 10

# the number of SD above mean to set the cutoff for selecting "hi windows"
sd_cut <- 3

# the number of "hi windows" a SNP needs to be in to be selected
#snp_wind_cut <- as.numeric(args[5])
snp_wind_cut <- 5

### SET OUTPUTS ###

#out_dir <- args[6]
out_dir <- paste(data_dir, 'MW_analysis/', sep = '')

#subgrp_name <- args[7]
subgrp_name <- 'MW8X_all'

#train_pop1 <- args[8]
train_pop1 <- 'MW'

#train_pop2 <- args[9]
train_pop2 <- 'GULF'

## Table to save
# F_val_tab_all_SNPs
fval_all_snps_file <- paste(out_dir, subgrp_name, '_Fval_', train_pop2, 
  '_into_', train_pop1, '_allsnps.txt', sep = '')

# F_val_tab_outlier_SNPs
fval_outlier_snp_file <- paste(out_dir, subgrp_name, '_Fval_', train_pop2,
  '_into_', train_pop1, '_outliersnps.txt', sep = '')

# F_val_tab_hi_window_SNPs
fval_hi_snp_file <- paste(out_dir, subgrp_name, '_Fval_', train_pop2,
  '_into_', train_pop1, '_topwindowsnps.txt', sep = '')

# window_tab
fval_window_full_file <- paste(out_dir, subgrp_name, '_Fval_', train_pop2,
  '_into_', train_pop1, '_', window_size, 'SNP', '_allwindows.txt', sep = '')

# hi_window
## include all positions in window
fval_hi_window_file <-  paste(out_dir, subgrp_name, '_Fval_', train_pop2,
  '_into_', train_pop1, '_', window_size, 'SNP', '_topwindows.txt', sep = '')

## Plot of window values
fval_window_plot_file <- paste(out_dir, subgrp_name, '_Fval_', train_pop2,
  '_into_', train_pop1, '_', window_size, 'SNP', '_windows.pdf', sep = '')

###############################

subgrp_geno_count <- lapply(strsplit(
  unlist(subgrp_hw[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

subgrp_n_genos <- unlist(lapply(subgrp_geno_count, function(x) sum(x)))
subgrp_low_inds <- which(subgrp_n_genos < (max(subgrp_n_genos)*miss_cut))
# [1] 433

subgrp_ref_freq <- unlist(lapply(subgrp_geno_count, function(x)
  (x[1] * 2 + x[2])/(sum(x)*2)))
subgrp_alt_freq <- unlist(lapply(subgrp_geno_count, function(x)
  (x[3] * 2 + x[2])/(sum(x)*2)))

###

# calculate the REF freq difference between pop2 and pop1 and convert to 
#  their range
ref_freq_dif <- unlist(train_freq[,4] - train_freq[,3])
ref_freq_range <- abs(ref_freq_dif)

# find SNPs where pop1 has higher REF freq
pop1_ref_hi <- which(ref_freq_dif < 0)

# calc diff in REF freq between subgrp and pop1
sub_v_pop1_ref <- subgrp_ref_freq - unlist(train_freq[,3])

# adjust diff for SNPs where pop1 has higher REF freq than pop2 because
#  want this to represent difference in direction of pop2
sub_v_pop1_ref_2 <- sub_v_pop1_ref
sub_v_pop1_ref_2[pop1_ref_hi] <- sub_v_pop1_ref[pop1_ref_hi]*-1

F_sub_v_pop1 <- sub_v_pop1_ref_2 / ref_freq_range
F_sub_v_pop1[subgrp_low_inds] <- NA

train_freq <- add_cumulative_pos(train_freq)

train_freq[, REF_STATE := allele_states[, list(REF_STATE)]]
train_freq[, ALT_STATE := allele_states[, list(ALT_STATE)]]

train_freq[, subgrp_ref_freq := subgrp_ref_freq]
train_freq[subgrp_low_inds, subgrp_ref_freq := NA]

train_freq[, F_subgrp_v_pop1 := F_sub_v_pop1] 

# find "outlier" SNPs
outlier_cut <- mean(train_freq$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(train_freq$F_subgrp_v_pop1, na.rm = T))
outlier_inds <- which(train_freq$F_subgrp_v_pop1 > outlier_cut)
outlier_tab <- train_freq[outlier_inds, ]

## Look at windows
# non-overlapping windows
subgrp_Fval_tab <- gen_Fval_nonoverlap_window_tab(
  ref_freq_tab = train_freq, window_size = window_size)

# overlapping windows
subgrp_Fval_tab_2 <- gen_Fval_overlap_window_tab(
  ref_freq_tab = train_freq, window_size = window_size)

### select top windows and SNPs
# for now, I am using overlapping windows

# set cutoff value
hi_cut_1 <- mean(subgrp_Fval_tab_2$F_SUB_V_POP1) + (sd_cut*
  sd(subgrp_Fval_tab_2$F_SUB_V_POP1))
# [1] 4.597645

hi_wind_inds <- which(subgrp_Fval_tab_2$F_SUB_V_POP1 > hi_cut_1)
# length(hi_wind_inds)
# 713

# setdiff(order(subgrp_Fval_tab_2$F_SUB_V_POP1, decreasing = T)[1:20], 
#   hi_wind_inds)

# connect windows back to SNPs
hi_snp_list <- list()
for(i in hi_wind_inds){
  tmp_full_ind <- which(train_freq$POS_CUM == subgrp_Fval_tab_2$POS_CUM[i])
  tmp_snp_inds <- c(tmp_full_ind:(tmp_full_ind + (window_size-1)))
  hi_snp_list[[i]] <- tmp_snp_inds
}

hi_snp_vec <- unlist(hi_snp_list)
# length(table(hi_snp_vec))
# [1] 7130

# table(table(hi_snp_vec))
#   1   2   3   4   5   6   7   8   9  10 
# 400 308 189 182 131 121  90  64  64 109
# maybe want to use the number of times a SNP shows up in a window as a cutoff?

test_snps <- as.numeric(names(table(hi_snp_vec))[
  which(table(hi_snp_vec) >= snp_wind_cut)])

hi_snp_tab <- train_freq[test_snps, ]

low_F_cut <- 0.5
low_F_snps <- c(which(hi_snp_tab$F_subgrp_v_pop1 < low_F_cut), 
  which(is.na(hi_snp_tab$F_subgrp_v_pop1)))

hi_snp_tab_2 <- hi_snp_tab[-low_F_snps]

### Save output files

# full F_val_table with all SNPs
fwrite(train_freq, file = fval_all_snps_file, sep = '\t')

# F_val_table with "outlier" SNPs using SD-based cutoff
fwrite(outlier_tab, file = fval_outlier_snp_file, sep = '\t')

# F_val_tabel with "hi SNP" repeatedly found in high windows and with 
#   F_val > low_F_cut
fwrite(hi_snp_tab_2, file = fval_hi_snp_file, sep = '\t')

# Full window result table
fwrite(subgrp_Fval_tab_2, file = fval_window_full_file, sep = '\t')

# Table of top windows
hi_window_tab <- subgrp_Fval_tab_2[hi_wind_inds]

hi_wind_full_pos_list <- lapply(hi_snp_list[hi_wind_inds], function(x)
  train_freq$POS[x])

hi_wind_full_pos_mat <- matrix(unlist(hi_wind_full_pos_list), 
  ncol = window_size, byrow = T)
colnames(hi_wind_full_pos_mat) <- paste('POS_WIND_SNP_', seq(window_size), 
  sep = '')

hi_wind_full_pos_tab <- as.data.table(hi_wind_full_pos_mat)

hi_window_tab_2 <- cbind(hi_window_tab, hi_wind_full_pos_tab)

fwrite(hi_window_tab_2, file = fval_hi_window_file, sep = '\t')

### Plotting

# Set Chromosome palette
chrom_palette <- rep(c('black', 'blue3'), times = 9)
names(chrom_palette) <- unique(subgrp_Fval_tab_2$CHROM)
chr_pal <- scale_colour_manual(name = 'CHROM', values = chrom_palette)

chr_max_pos <- c()
for(chrm in unique(subgrp_Fval_tab_2$CHROM)){
  chr_max_pos[[chrm]] <- max(subgrp_Fval_tab_2$POS_CUM[
    subgrp_Fval_tab_2$CHROM == chrm])
}
chr_max_pos_tmp <- unlist(chr_max_pos)
chr_max_pos <- chr_max_pos_tmp[c(1:17)]

chr_names <- unique(subgrp_Fval_tab_2$CHROM)

chr_mid_pos <- rep(as.numeric(NA), times = 18)
chr_mid_pos[1] <- (((chr_max_pos_tmp[1] - min(subgrp_Fval_tab_2$POS_CUM))/2) + 
  min(subgrp_Fval_tab_2$POS_CUM))
for(i in c(2:18)){
  chr_mid_pos[i] <- (((chr_max_pos_tmp[i] - chr_max_pos_tmp[i-1])/2) + 
    chr_max_pos_tmp[i-1])
}


gg_fval <- ggplot(subgrp_Fval_tab_2, aes(x = POS_CUM, y = F_SUB_V_POP1)) +
  geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
  chr_pal +
  ggtitle(paste(subgrp_name, 'Fval in', window_size, 'SNP windows for', 
    train_pop2, 'into', train_pop1, sep = ' ')) +
  ylab('CUMMULATIVE\nF_VAL') +
  scale_x_continuous(name = '', breaks = chr_mid_pos, 
    labels = chr_names) + 
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed') +
  geom_hline(yintercept = hi_cut_1, color = 'red', linetype = 'dotted')

pdf(fval_window_plot_file, width = 18, height = 4)
gg_fval
dev.off()


#test_hist_out <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_03_mw_v_gulf_Fval_histogram.pdf'

#gg_hist <- ggplot(mw_03_Fval_tab_2, aes(x = F_SUB_V_POP1)) +
#  geom_density()

#pdf(test_hist_out, width = 5, height = 5)
#gg_hist
#dev.off()

dev.off()


