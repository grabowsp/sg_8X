# Make manhattan plot of Fval across chromosome colored by RDA significance

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD LIBRARIES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###
fval_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/', 
  'MW8X_all_Fval_GULF_into_MW_allsnps.txt', sep = '')
fvals <- fread(fval_file)

atl_fval_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/',
  'MW8X_all_Fval_ATL_into_MW_allsnps.txt', sep = '')
atl_fvals <- fread(atl_fval_file)

rda_res_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/',
  'Results_Gulf_Midwest_8x_4x_rda.csv', sep = '')
rda_res <- fread(rda_res_file)
# need to transfer from Downloads to HA, then load

atl_rda_res_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/',
  'Results_Atlantic_Midwest_8x_4x_rda.csv', sep = '')
atl_rda_res <- fread(atl_rda_res_file)

#allele_state_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/',
#  'MW_v_GULF_allele_states.txt', sep = '')
#allele_states <- fread(allele_state_file)

#train_freq_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/',
#  'MW_v_GULF_ref_freq.txt', sep = '')
#train_freq <- fread(train_freq_file)

### SET OUTPUT ###

out_dir <- paste('/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/',
  'MW_analysis/', sep = '')

mw_fval_manhat_out <- paste(out_dir, 'MW8X_all_Gulf_MW_Fval_manhattan_v1.png',
  sep = '')

mw_atl_fval_manhat_out <- paste(out_dir,
  'MW8X_all_Atlantic_MW_Fval_manhattan_v1.png', sep = '')


### SET VARIABLES ###
sd_cut <- 3

##########

## Determine cutoff and rda significant SNPs
# Gulf
outlier_cut <- mean(fvals$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(fvals$F_subgrp_v_pop1, na.rm = T))
outlier_inds <- which(fvals$F_subgrp_v_pop1 > outlier_cut)

fvals[, SNP_NAME := paste(fvals$CHR, fvals$POS, sep = '_')]
rda_overlap_inds <- which(fvals$SNP_NAME %in% rda_res$snp)

fvals[, PLOT_CLASS := fvals$CHR]
fvals[rda_overlap_inds, PLOT_CLASS := 'RDA_SIG']
fvals[intersect(outlier_inds, rda_overlap_inds), PLOT_CLASS := 'OVERLAP']

fvals[, EXP_FVAL := exp(fvals$F_subgrp_v_pop1)]

# Atlantic
atl_outlier_cut <- mean(atl_fvals$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(atl_fvals$F_subgrp_v_pop1, na.rm = T))
atl_outlier_inds <- which(atl_fvals$F_subgrp_v_pop1 > atl_outlier_cut)

atl_fvals[, SNP_NAME := paste(atl_fvals$CHR, atl_fvals$POS, sep = '_')]
atl_rda_overlap_inds <- which(atl_fvals$SNP_NAME %in% atl_rda_res$snp)

atl_fvals[, PLOT_CLASS := atl_fvals$CHR]
atl_fvals[atl_rda_overlap_inds, PLOT_CLASS := 'RDA_SIG']
atl_fvals[intersect(atl_outlier_inds, atl_rda_overlap_inds), 
  PLOT_CLASS := 'OVERLAP']

atl_fvals[, EXP_FVAL := exp(atl_fvals$F_subgrp_v_pop1)]


## Plotting
# Gulf
chrom_palette <- rep(c('black', 'blue3'), times = 9)
names(chrom_palette) <- unique(fvals$CHR)
chr_pal <- scale_colour_manual(name = 'CHR', values = chrom_palette)

sub_class_palette <- c('orange3', 'red2')
names(sub_class_palette) <- c('RDA_SIG', 'OVERLAP')
plot_class_palette <- c(chrom_palette, sub_class_palette)
plot_class_pal <- scale_colour_manual(name = 'PLOT_CLASS', 
  values = plot_class_palette)

chr_max_pos <- c()
for(chrm in unique(fvals$CHR)){
  chr_max_pos[[chrm]] <- max(fvals$POS_CUM[
    fvals$CHR == chrm])
}
chr_max_pos_tmp <- unlist(chr_max_pos)
chr_max_pos <- chr_max_pos_tmp[c(1:17)]

chr_names <- unique(fvals$CHR)

chr_mid_pos <- rep(as.numeric(NA), times = 18)
chr_mid_pos[1] <- (((chr_max_pos_tmp[1] - min(fvals$POS_CUM))/2) +
  min(fvals$POS_CUM))
for(i in c(2:18)){
  chr_mid_pos[i] <- (((chr_max_pos_tmp[i] - chr_max_pos_tmp[i-1])/2) +
    chr_max_pos_tmp[i-1])
}

gg_fval <- ggplot(fvals, aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = PLOT_CLASS), show.legend = F) +
  plot_class_pal +
  ggtitle(paste('MW8X_all exp(Fvals) for Gulf into Midwest SNPs,', 
    '\ncolored by RDA significance', sep = ' ')) +
  ylab('EXP(F_VAL)') +
  scale_x_continuous(name = '', breaks = chr_mid_pos, labels = chr_names) +
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed') +
  geom_hline(yintercept = exp(outlier_cut), color = 'red', linetype = 'dotted')

png(mw_fval_manhat_out, width = 480*4, height = 480)
gg_fval
dev.off()

# Atlantic
chrom_palette <- rep(c('black', 'blue3'), times = 9)
names(chrom_palette) <- unique(atl_fvals$CHR)
chr_pal <- scale_colour_manual(name = 'CHR', values = chrom_palette)

sub_class_palette <- c('orange3', 'red2')
names(sub_class_palette) <- c('RDA_SIG', 'OVERLAP')
plot_class_palette <- c(chrom_palette, sub_class_palette)
plot_class_pal <- scale_colour_manual(name = 'PLOT_CLASS',
  values = plot_class_palette)

chr_max_pos <- c()
for(chrm in unique(atl_fvals$CHR)){
  chr_max_pos[[chrm]] <- max(atl_fvals$POS_CUM[
    atl_fvals$CHR == chrm])
}
chr_max_pos_tmp <- unlist(chr_max_pos)
chr_max_pos <- chr_max_pos_tmp[c(1:17)]

chr_names <- unique(atl_fvals$CHR)

chr_mid_pos <- rep(as.numeric(NA), times = 18)
chr_mid_pos[1] <- (((chr_max_pos_tmp[1] - min(atl_fvals$POS_CUM))/2) +
  min(atl_fvals$POS_CUM))
for(i in c(2:18)){
  chr_mid_pos[i] <- (((chr_max_pos_tmp[i] - chr_max_pos_tmp[i-1])/2) +
    chr_max_pos_tmp[i-1])
}

atl_gg_fval <- ggplot(atl_fvals, aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = PLOT_CLASS), show.legend = F) +
  plot_class_pal +
  ggtitle(paste('MW8X_all exp(Fvals) for Atlantic into Midwest SNPs,',
    '\ncolored by RDA significance', sep = ' ')) +
  ylab('EXP(F_VAL)') +
  scale_x_continuous(name = '', breaks = chr_mid_pos, labels = chr_names) +
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed') +
  geom_hline(yintercept = exp(atl_outlier_cut), color = 'red', 
    linetype = 'dotted')

png(mw_atl_fval_manhat_out, width = 480*4, height = 480)
atl_gg_fval
dev.off()

quit(save = 'no')



