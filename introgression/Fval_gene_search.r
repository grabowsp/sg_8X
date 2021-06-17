# Make manhattan plot of Fval across chromosome colored by RDA significance

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD LIBRARIES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###
# I will need the annotation file to look at the genes, eventually

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

rda_outlier_inds <- intersect(outlier_inds, rda_overlap_inds)

fval_gulf_4K_peak <- fvals[intersect(
  rda_outlier_inds, which(fvals$CHR == 'Chr04K')),]
max(fval_gulf_4K_peak$POS) - min(fval_gulf_4K_peak$POS)
# 2,047,292 = 2MB region; a lot of genes; some of it is gene poor.
#  it looks like there are also a bunch of RDA-sig SNPs just below the
#  F_val outlier cutoff

fval_gulf_3N_peak <- fvals[intersect(
  rda_outlier_inds, which(fvals$CHR == 'Chr03N')),]
# 56452006-56280071 = 171935; 172kb region
# Includes 3 aspartyl proteases

fval_gulf_7N_hits <- fvals[intersect(
  rda_outlier_inds, which(fvals$CHR == 'Chr07N')),]

fval_gulf_7N_peak_1 <- fval_gulf_7N_hits[c(2:11),]
# 4640625 - 4058205 = 582420; 580kb region
# huge region, so a bunch of genes

fval_gulf_7N_peak_2 <- fval_gulf_7N_hits[c(13:15)]
# 9339584 - 9263785 = 75799; 76 kb
# mainly a gene desert

fval_gulf_7N_peak_3 <- fval_gulf_7N_hits[c(18:34)]
# 27351270 - 16953170 = 10398100; 10 GB region 

## Look at the F_val peaks on Chr07K and Chr07N
fval_gulf_7K_all <- fvals[intersect(outlier_inds,
  which(fvals$CHR == 'Chr07K'))]

gulf_7K_1 <- fval_gulf_7K_all[c(8:14), ]
# 3382320-3057495 = 324825; 325kb

fval_gulf_7N_all <- fvals[intersect(outlier_inds,
  which(fvals$CHR == 'Chr07N'))]

gulf_7N_1 <- fval_gulf_7N_all[c(2:15), ]
# 4436027 - 4023356 = 412671; 412kb
# this covers a homologous region to gulf_7K_1...

###
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

###
atl_rda_outlier_inds <- intersect(atl_outlier_inds, atl_rda_overlap_inds)

fval_atl_3K_peak <- atl_fvals[intersect(
  atl_rda_outlier_inds, which(atl_fvals$CHR == 'Chr03K')),]
# 52026630 - 51871367 = 155263; 155kb
# 7 genes, including 3 Aspartyl proteases
# this looks like the homeologous interval as the Chr03N peak...

fval_atl_3N_peak <- atl_fvals[intersect(
  atl_rda_outlier_inds, which(atl_fvals$CHR == 'Chr03N')),]
# 56452006 - 56281436 = 170570, 170kb
# homeologous interval as the Chr03K peak...

# makes me wonder if there's some sort of 


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



