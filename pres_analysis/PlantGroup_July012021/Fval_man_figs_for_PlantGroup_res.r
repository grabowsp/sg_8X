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

### SET OUTPUT ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/'

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

### Plotting stuff
chrom_palette <- rep(c('black', 'blue3'), times = 9)
names(chrom_palette) <- unique(fvals$CHR)
chr_pal <- scale_colour_manual(name = 'CHR', values = chrom_palette)

sub_class_palette <- c('orange3', 'red2')
names(sub_class_palette) <- c('RDA_SIG', 'OVERLAP')
plot_class_palette <- c(chrom_palette, sub_class_palette)
plot_class_pal <- scale_colour_manual(name = 'PLOT_CLASS',
  values = plot_class_palette)

### Atlantic Full Manhattan, no RDA results

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

atl_gg_noRDA_fval <- ggplot(atl_fvals, aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = CHR), show.legend = F) +
  plot_class_pal +
  ggtitle('MW8X_all exp(Fvals) for Atlantic into Midwest SNPs') +
  ylab('EXP(F_VAL)') +
  scale_x_continuous(name = '', breaks = chr_mid_pos, labels = chr_names) +
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed') +
  geom_hline(yintercept = exp(atl_outlier_cut), color = 'red',
    linetype = 'dotted')

atl_noRDA_man_out <- paste(out_dir, 'ATL_into_MW_manhattan_noRDA.png', sep = '')

png(atl_noRDA_man_out, width = 480*4, height = 480)
atl_gg_noRDA_fval
dev.off()

### Atlantic Chr08N manhattan, no RDA results

atl_gg_Chr08N_noRDA_fval <- ggplot(atl_fvals[CHR == 'Chr08N'], 
  aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = CHR), show.legend = F) +
  plot_class_pal +
  ggtitle('MW8X_all exp(Fvals) for Atlantic\ninto Midwest SNPs\nChr08N') +
  ylab('EXP(F_VAL)') +
#  scale_x_continuous(labels = '') +
  geom_hline(yintercept = exp(atl_outlier_cut), color = 'red', 
    linetype = 'dotted')

atl_Chr08N_noRDA_man_out <- paste(out_dir, 
  'ATL_into_MW_manhattan_noRDA_Chr08N.png', sep = '')

png(atl_Chr08N_noRDA_man_out, width = 480/2, height = 480)
atl_gg_Chr08N_noRDA_fval
dev.off()

### Atlantic Chr03K manhattan, no RDA results - as contrast to Chr08N

atl_gg_Chr03K_noRDA_fval <- ggplot(atl_fvals[CHR == 'Chr03K'],
  aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = CHR), show.legend = F) +
  plot_class_pal +
  ggtitle('MW8X_all exp(Fvals) for Atlantic\ninto Midwest SNPs\nChr03K') +
  ylab('EXP(F_VAL)') +
#  scale_x_continuous(labels = '') +
  geom_hline(yintercept = exp(atl_outlier_cut), color = 'red',
    linetype = 'dotted')

atl_Chr03K_noRDA_man_out <- paste(out_dir,
  'ATL_into_MW_manhattan_noRDA_Chr03K.png', sep = '')

png(atl_Chr03K_noRDA_man_out, width = 480/2, height = 480)
atl_gg_Chr03K_noRDA_fval
dev.off()


#### Test for significance for enrichment of RDA outliers in top Fval SNPs

outlier_tab <- data.frame(rda_sig = c(105, 6135), not_sig = c(492, 52567))
rownames(outlier_tab) <- c('Low_shift', 'All_SNPs')

chisq.test(outlier_tab)
# X-squared = 31.217, df = 1, p-value = 2.308e-08

### Atlantic Chr08N Manhattan, with RDA results

atl_gg_Chr08N_withRDA_fval <- ggplot(atl_fvals[CHR == 'Chr08N'],
  aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = PLOT_CLASS), show.legend = F) +
  plot_class_pal +
  ggtitle('MW8X_all exp(Fvals) for Atlantic\ninto Midwest SNPs\nChr08N') +
  ylab('EXP(F_VAL)') +
#  scale_x_continuous(labels = '') +
  geom_hline(yintercept = exp(atl_outlier_cut), color = 'red',
    linetype = 'dotted')

atl_Chr08N_withRDA_man_out <- paste(out_dir,
  'ATL_into_MW_manhattan_withRDA_Chr08N.png', sep = '')

png(atl_Chr08N_withRDA_man_out, width = 480/2, height = 480)
atl_gg_Chr08N_withRDA_fval
dev.off()






