# Look at MW_03 allele frequencies and look movement towards GULF freqs

# Steps
# calculate F_value
# Generate window values
# Select top windows using sd-based cutoff
# Find SNPs that are found X+ times (ex 5) in top windows
# Those SNPs/loci are candidate regions

# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)
library(gridExtra)

pop_intro_func_file <- paste('/home/grabowsky/tools/workflows/sg_8X/', 
  'introgression/chrom_painting_files/', 'chrom_paint_pop_functions.r', 
  sep = '')
source(pop_intro_func_file)

### INPUT DATA ###
data_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'

subgrp_hw_file <- paste(data_dir, 'MW_03_MWvGULF.hwe', sep = '')
subgrp_hw <- fread(subgrp_hw_file)

train_freq_file <- paste(data_dir, 'MW_v_GULF_ref_freq.txt', sep = '')
train_freq <- fread(train_freq_file)

allele_state_file <- paste(data_dir, 'MW_v_GULF_allele_states.txt', sep = '')
allele_states <- fread(allele_state_file)

##########

subgrp_geno_count <- lapply(strsplit(
  unlist(subgrp_hw[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

subgrp_n_genos <- unlist(lapply(subgrp_geno_count, function(x) sum(x)))
subgrp_low_inds <- which(subgrp_n_genos < (max(subgrp_n_genos)*0.8))
# [1] 667

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

train_freq[, subgrp_ref_freq := subgrp_ref_freq]
train_freq[subgrp_low_inds, subgrp_ref_freq := NA]

train_freq[, F_subgrp_v_pop1 := F_sub_v_pop1] 

test <- gen_Fval_nonoverlap_window_tab(
  ref_freq_tab = train_freq,
#  chrom = 'Chr01K',
  window_size = 10
  ) 

mw_03_Fval_tab <- gen_Fval_nonoverlap_window_tab(
  ref_freq_tab = train_freq, window_size = 10)

mw_03_Fval_tab_2 <- gen_Fval_overlap_window_tab(
  ref_freq_tab = train_freq, window_size = 10)

# Set Chromosome palette
chrom_palette <- rep(c('black', 'blue3'), times = 9)
names(chrom_palette) <- unique(mw_03_Fval_tab_2$CHROM)
chr_pal <- scale_colour_manual(name = 'CHROM', values = chrom_palette)

chr_max_pos <- c()
for(chrm in unique(mw_03_Fval_tab_2$CHROM)){
  chr_max_pos[[chrm]] <- max(mw_03_Fval_tab_2$POS_CUM[
    mw_03_Fval_tab_2$CHROM == chrm])
}
chr_max_pos <- unlist(chr_max_pos)
chr_max_pos <- chr_max_pos[c(1:17)]

gg_fval <- ggplot(mw_03_Fval_tab_2, aes(x = POS_CUM, y = F_SUB_V_POP1)) +
  geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
  chr_pal +
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed')

test_out <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_03_mw_v_gulf_Fval.pdf'

pdf(test_out, width = 18, height = 4)
gg_fval
dev.off()

test_hist_out <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_03_mw_v_gulf_Fval_histogram.pdf'

gg_hist <- ggplot(mw_03_Fval_tab_2, aes(x = F_SUB_V_POP1)) +
  geom_density()

pdf(test_hist_out, width = 5, height = 5)
gg_hist
dev.off()



ref_freq_tab <- train_freq
chrom = 'Chr01K'
window_size = 10

func_file <- paste('/home/grabowsky/tools/workflows/sg_8X/introgression/', 
  'chrom_painting_files/', 'chrom_paint_pop_functions.r', sep = '')
source(func_file)



# CONTINUE FROM HERE WITH WINDOW ANALYSIS
### decide on window cutoff of F_sub_v_pop1


hi_cut_1 <- mean(mw_03_Fval_tab_2$F_SUB_V_POP1) + (3*
  sd(mw_03_Fval_tab_2$F_SUB_V_POP1))
# [1] 4.597645

hi_wind_inds <- which(mw_03_Fval_tab_2$F_SUB_V_POP1 > hi_cut_1)
length(hi_wind_inds)
# 650

# need to connect windows back to SNPs
hi_snp_list <- list()

for(i in hi_wind_inds){
  tmp_full_ind <- which(train_freq$POS_CUM == mw_03_Fval_tab_2$POS_CUM[i])
  tmp_snp_inds <- c(tmp_full_ind:(tmp_full_ind + (window_size-1)))
  hi_snp_list[[i]] <- tmp_snp_inds
}

hi_snp_vec <- unlist(hi_snp_list)
length(table(hi_snp_vec))
# [1] 1658

table(table(hi_snp_vec))
#   1   2   3   4   5   6   7   8   9  10 
# 400 308 189 182 131 121  90  64  64 109

# maybe want to use the number of times a SNP shows up in a window as a cutoff?


test_snps <- as.numeric(names(table(hi_snp_vec))[
  which(table(hi_snp_vec) >= 5)])

tmp_1 <- order(train_freq$F_subgrp_v_pop1, decreasing = T)[1:100]




tmp_full_ind <- which(train_freq$POS_CUM == mw_03_Fval_tab_2$POS_CUM[
  hi_wind_inds[455]])
tmp_snp_inds <- c(tmp_full_ind:(tmp_full_ind + (window_size-1)))
hi_snp_list[[hi_wind_inds[455]]] <- tmp_snp_inds


which(train_freq$CHR == mw_03_Fval_tab_2$CHROM[hi_wind_inds[1]] &
  train_freq$POS)

