# Manhattan plot of Fval colored by RDA significance

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD LIBRARIES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###
file_dir <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_fval_files/'

fval_file <- paste(file_dir, 'MW8X_Fval_ATL_into_MW_allsnps_Fst0.5.txt', 
  sep = '')
fvals <- fread(fval_file)

rda_res_file <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_rda_results/Updated_Results_Atlantic_Midwest_8x_4x_rda_JN.csv'
rda_res <- fread(rda_res_file)

rda_snp_file <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_ATL_into_MW_Fst_0.5.RDA_SNP_info.txt'
rda_snps <- fread(rda_snp_file)

### SET OUTPUT ###

manhattan_out <- paste(file_dir, 'MW8X_MWvATL_Fst0.5_Fval_RDA_color.png',
  sep = '')

manhattan_2_out <- paste(file_dir, 
  'MW8X_MWvATL_Fst0.5_Fval_RDA_color_RDAonly.png',
  sep = '')

rda_fval_overlap_out <- paste(file_dir,
  'MW8X_MWvATL_Fst0.5_RDA_Fval_overlap.txt', sep = '')

rda_fval_2_overlap_out <- paste(file_dir,
  'MW8X_MWvATL_Fst0.5_RDA_Fval_overlap_RDASNPsonly.txt', sep = '')

### SET VARIABLES ###
sd_cut <- 3

## Determine cutoff and rda significant SNPs
outlier_cut <- mean(fvals$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(fvals$F_subgrp_v_pop1, na.rm = T))
outlier_inds <- which(fvals$F_subgrp_v_pop1 > outlier_cut)
# 607 total

fvals[, SNP_NAME := paste(fvals$CHR, fvals$POS, sep = '_')]
rda_sig_inds <- which(fvals$SNP_NAME %in% rda_res$snp)

rda_test_inds <- setdiff(which(fvals$SNP_NAME %in% rda_snps$ID), rda_sig_inds)

rda_overlap_inds <- intersect(outlier_inds, rda_sig_inds)
# 60

fvals[, PLOT_CLASS := fvals$CHR]
fvals[rda_sig_inds, PLOT_CLASS := 'RDA_SIG']
fvals[rda_test_inds, PLOT_CLASS := 'RDA_TEST']
fvals[rda_overlap_inds, PLOT_CLASS := 'OVERLAP']

fvals[, EXP_FVAL := exp(fvals$F_subgrp_v_pop1)]

## Plotting
chrom_palette <- rep(c('black', 'blue3'), times = 9)
names(chrom_palette) <- unique(fvals$CHR)
chr_pal <- scale_colour_manual(name = 'CHR', values = chrom_palette)

sub_class_palette <- c('yellow4', 'orange3', 'red2')
names(sub_class_palette) <- c('RDA_TEST', 'RDA_SIG', 'OVERLAP')

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
  ggtitle(paste('MW8X exp(Fvals) for Atlantic into Midwest SNPs,',
    '\ncolored by RDA significance', sep = ' ')) +
  ylab('EXP(F_VAL)') +
  scale_x_continuous(name = '', breaks = chr_mid_pos, labels = chr_names) +
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed') +
  geom_hline(yintercept = exp(outlier_cut), color = 'red', linetype = 'dotted')


png(manhattan_out, width = 480*4, height = 480)
gg_fval
dev.off()

fwrite(fvals, rda_fval_overlap_out)

# Now just with SNPs used in the RDA

fvals_2 <- fvals[which(fvals$SNP_NAME %in% rda_snps$ID), ]

fvals_2[, PLOT_CLASS := fvals_2$CHR]

outlier_2_cut <- mean(fvals_2$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(fvals_2$F_subgrp_v_pop1, na.rm = T))

outlier_2_inds <- which(fvals_2$F_subgrp_v_pop1 > outlier_2_cut)
# 262 total

rda_sig_inds_2 <- which(fvals_2$SNP_NAME %in% rda_res$snp)
rda_overlap_inds_2 <- intersect(outlier_2_inds, rda_sig_inds_2)
#125

fvals_2[rda_sig_inds_2, PLOT_CLASS := 'RDA_SIG']
fvals_2[rda_overlap_inds_2, PLOT_CLASS := 'OVERLAP']

gg_fval_2 <- ggplot(fvals_2, aes(x = POS_CUM, y = EXP_FVAL)) +
  geom_point(aes(color = PLOT_CLASS), show.legend = F) +
  plot_class_pal +
  ggtitle(paste('MW8X exp(Fvals) for Atlantic into Midwest SNPs,',
    '\nSNPs used in RDA, colored by RDA significance', sep = ' ')) +
  ylab('EXP(F_VAL)') +
  scale_x_continuous(name = '', breaks = chr_mid_pos, labels = chr_names) +
  geom_vline(xintercept = chr_max_pos, color = 'gray50', linetype = 'dashed') +
  geom_hline(yintercept = exp(outlier_2_cut), color = 'red', linetype = 'dotted')

png(manhattan_2_out, width = 480*4, height = 480)
gg_fval_2
dev.off()

fwrite(fvals_2, rda_fval_2_overlap_out)


