# Analysis of PCoA reproducibility - on laptop

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
# install.packages('patchwork')
library(patchwork)

### INPUT DATA ###
pcoa_list_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2.100k.pcoa_mats.list.rds'
pcoa_list <- readRDS(pcoa_list_file)

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

### SET OUTPUT ###
out_dir <- '/Users/grabowsk/Documents/Switchgrass_8X/pcoa_reproducibility/'

### SET VARIABLES ###
tet_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '4X')]
oct_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '8X' |
  samp_meta$NQUIRE_PLOIDY == '6X')]

### ANALYSIS FUNCTIONS ###
standardize_pco <- function(pco_mat, pco_num, control_name){
  # Standardize a PCo vector from 0 to 1; includes setting a specific sample
  #   as always having a positive value, so that PCo's from replicate runs
  #   (should) have the same sign
  # INPUTS
  # pco_mat = results from cmdscale or some other PCoA function
  # pco_num = the PCo number to analyze
  # control_name = the "control" sample that will always be positive
  # OUTPUT
  # vector of standardized PCo$pco_num values
  ###########
  control_ind <- which(rownames(pco_mat) == control_name)
  tmp_vec <- pco_mat[, pco_num]
  if(tmp_vec[control_ind] < 0){
    tmp_vec <- tmp_vec * -1
  }
  tmp_vec <- tmp_vec - min(tmp_vec)
  tmp_vec <- tmp_vec / max(tmp_vec)
  return(tmp_vec)
}

gen_pco_comp_mat <- function(pco_res_list, pco_num, control_name){
  tmp_pco_list <- list()
  for(dpl in seq(length(pco_res_list))){
    tmp_pco_list[[dpl]] <- standardize_pco(pco_mat = pco_res_list[[dpl]],
      pco_num = pco_num, control_name = control_name)
  }
  tmp_pco_mat <- matrix(unlist(tmp_pco_list), nrow = length(tmp_pco_list),
    byrow = T)
  colnames(tmp_pco_mat) <- names(tmp_pco_list[[1]])
  return(tmp_pco_mat)
}

################
### Make PCoA result lists

dip_pcoa_list <- pcoa_list[['dip_pcoa_list']]
tet_pcoa_list <- pcoa_list[['tet_pcoa_list']]
poly_pcoa_list <- pcoa_list[['poly_pcoa_list']] 

comp_res_list <- data.table(samp_name = rownames(dip_pcoa_list[[1]]))
comp_res_list[, ploidy := as.character(NA)]
comp_res_list[which(comp_res_list$samp_name %in% tet_names), ploidy := '4X']
comp_res_list[which(comp_res_list$samp_name %in% oct_names), ploidy := '8X']

cr_tet_inds <- which(comp_res_list$ploidy == '4X')
cr_oct_inds <- which(comp_res_list$ploidy == '8X')

### PCo1 ###

dip_pco_1_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list, 
  pco_num = 1, control_name = 'J430.A')
comp_res_list[, dip_pco_1_var := apply(dip_pco_1_stand_mat, 2, var)]
comp_res_list[, dip_pco_1_mean := apply(dip_pco_1_stand_mat, 2, mean)]

tet_pco_1_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
  pco_num = 1, control_name = 'J430.A')
comp_res_list[, tet_pco_1_var := apply(tet_pco_1_stand_mat, 2, var)]
comp_res_list[, tet_pco_1_mean := apply(tet_pco_1_stand_mat, 2, mean)]

poly_pco_1_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
  pco_num = 1, control_name = 'J430.A')
comp_res_list[, poly_pco_1_var := apply(poly_pco_1_stand_mat, 2, var)]
comp_res_list[, poly_pco_1_mean := apply(poly_pco_1_stand_mat, 2, mean)]

## Mean ##
comp_res_list[ , tet_dif_dip_pco_1_mean := 
  comp_res_list$tet_pco_1_mean - comp_res_list$dip_pco_1_mean]
summary(abs(comp_res_list$tet_dif_dip_pco_1_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.240e-06 1.049e-03 2.482e-03 6.331e-03 7.004e-03 7.689e-02
summary(abs((comp_res_list$tet_dif_dip_pco_1_mean[cr_tet_inds])))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.243e-06 8.232e-04 1.604e-03 2.685e-03 4.037e-03 2.216e-02
summary(abs((comp_res_list$tet_dif_dip_pco_1_mean[cr_oct_inds])))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3.339e-05 1.296e-02 1.885e-02 1.857e-02 2.305e-02 7.689e-02
# Generally small difference between means in 4X individuals when using
#   tetrasomic vs disomc genotypes; difference an order of magnitude 
#   bigger for 8X
t.test(x = abs((comp_res_list$tet_dif_dip_pco_1_mean[cr_tet_inds])),
  y = abs((comp_res_list$tet_dif_dip_pco_1_mean[cr_oct_inds])))
# t = -21.051, df = 168.3, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.01737998 -0.01439975
# sample estimates:
#   mean of x   mean of y 
# 0.002684845 0.018574709 

tet_dip_mean_1_melt <- data.table(tet_dip_dif = c(
    abs(comp_res_list$tet_dif_dip_pco_1_mean[cr_tet_inds]),
    abs(comp_res_list$tet_dif_dip_pco_1_mean[cr_oct_inds])
  ),
  comp = c(
    rep('4X', times = length(cr_tet_inds)),
    rep('8X', times = length(cr_oct_inds))
  )
)

gg_pco1_tet_v_dip_mean_density <- ggplot(tet_dip_mean_1_melt) +
  geom_density(aes(x = tet_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.09)) +
  ggtitle(paste('Difference of mean PCo1 values for\ndisomic vs ', 
                'tetrasomic genotypes', sep = ''))

gg_pco1_tet_v_dip_mean_box <- ggplot(tet_dip_mean_1_melt) +
  geom_boxplot(aes(x = comp, y = tet_dip_dif, fill = comp)) +
  ggtitle(paste('Difference of mean PCo1 values for\ndisomic vs ', 
                'tetrasomic genotypes', sep = '')) +
  coord_cartesian(ylim = c(0, 0.08))

comp_res_list[ , poly_dif_dip_pco_1_mean := 
                 comp_res_list$poly_pco_1_mean - comp_res_list$dip_pco_1_mean]
summary(abs(comp_res_list$poly_dif_dip_pco_1_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.240e-06 5.772e-04 1.254e-03 6.067e-03 4.159e-03 7.264e-02
summary(abs(comp_res_list$poly_dif_dip_pco_1_mean[cr_tet_inds]))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.236e-06 4.585e-04 9.293e-04 1.217e-03 1.590e-03 5.697e-03
summary(abs(comp_res_list$poly_dif_dip_pco_1_mean[cr_oct_inds]))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003526 0.016690 0.022552 0.022357 0.028462 0.072642

# compare tet_v_dip and poly_v_dip for 4X
t.test( x = abs(comp_res_list$tet_dif_dip_pco_1_mean[cr_tet_inds]), 
        y = abs(comp_res_list$poly_dif_dip_pco_1_mean[cr_tet_inds]))
# t = 12.096, df = 727.48, p-value < 2.2e-16
# sample estimates:
#   mean of x   mean of y 
# 0.002684845 0.001216612
# polyploid genotype produce smaller difference from dip than do tetrasomic
#   genotypes

comp_res_list[ , poly_dif_tet_pco_1_mean := 
                 comp_res_list$poly_pco_1_mean - comp_res_list$tet_pco_1_mean]
summary(abs(comp_res_list$poly_dif_tet_pco_1_mean))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.793e-06 9.371e-04 2.149e-03 3.858e-03 7.094e-03 2.171e-02
summary(abs(comp_res_list$poly_dif_tet_pco_1_mean[cr_tet_inds]))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.793e-06 7.737e-04 1.601e-03 3.151e-03 4.962e-03 2.171e-02
summary(abs(comp_res_list$poly_dif_tet_pco_1_mean[cr_oct_inds]))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3.261e-05 4.391e-03 7.149e-03 6.230e-03 7.740e-03 1.186e-02

poly_tet_mean_1_melt <- data.table(poly_tet_dif = c(
    abs(comp_res_list$poly_dif_tet_pco_1_mean[cr_tet_inds]),
    abs(comp_res_list$poly_dif_tet_pco_1_mean[cr_oct_inds])
  ),
  comp = c(
    rep('4X', times = length(cr_tet_inds)),
    rep('8X', times = length(cr_oct_inds))
  )
)

gg_pco1_poly_v_tet_mean_density <- ggplot(poly_tet_mean_1_melt) +
  geom_density(aes(x = poly_tet_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.09)) +
  ggtitle(paste('Difference of mean PCo1 values for\npolyploid vs ', 
                'tetrasomic genotypes', sep = ''))

gg_pco1_poly_v_tet_mean_box <- ggplot(poly_tet_mean_1_melt) +
  geom_boxplot(aes(x = comp, y = poly_tet_dif, fill = comp)) +
  ggtitle(paste('Difference of mean PCo1 values for\npolyploid vs ', 
                'tetrasomic genotypes', sep = '')) +
  coord_cartesian(ylim = c(0, 0.08))

poly_ploidy_mean_1_melt <- data.table(poly_ploidy_dif = c(
  abs(comp_res_list$poly_dif_dip_pco_1_mean[cr_tet_inds]),
  abs(comp_res_list$poly_dif_tet_pco_1_mean[cr_oct_inds])
  ),
  comp = c(
   rep('4X', times = length(cr_tet_inds)), 
   rep('8X', times = length(cr_oct_inds))
     )
)

gg_pco1_poly_v_ploidy_mean_density <- ggplot(poly_ploidy_mean_1_melt) +
  geom_density(aes(x = poly_ploidy_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.09)) + 
  ggtitle(paste('Difference of mean PCo1 values for\npolyploid vs ', 
                'ploidy-appropriate genotypes', sep = ''))

gg_pco1_poly_v_ploidy_mean_box <- ggplot(poly_ploidy_mean_1_melt) +
  geom_boxplot(aes(x = comp, y = poly_ploidy_dif, fill = comp)) +
  ggtitle(paste('Difference of mean PCo1 values for\npolyploid vs ', 
                'ploidy-appropriate genotypes', sep = '')) +
  coord_cartesian(ylim = c(0, 0.08))

# Do polyploid genotypes result in same differences from ploidy-appropriate
#   results in each ploidy level
t.test(x = abs(comp_res_list$poly_dif_dip_pco_1_mean[cr_tet_inds]),
       y = abs(comp_res_list$poly_dif_tet_pco_1_mean[cr_oct_inds]))
# t = -22.209, df = 175.72, p-value < 2.2e-16
# sample estimates:
# mean of x   mean of y 
# 0.001216612 0.006230198
# Not quite the same, but much closer than before

pco_1_mean_pdf_name <- paste(out_dir, 'natv2_pco_1_mean_reproducibility.pdf',
                             sep = '') 
pdf(pco_1_mean_pdf_name, height = 10, width = 10)
(gg_pco1_tet_v_dip_mean_density + gg_pco1_tet_v_dip_mean_box) /
  (gg_pco1_poly_v_tet_mean_density + gg_pco1_poly_v_tet_mean_box) /
  (gg_pco1_poly_v_ploidy_mean_density + gg_pco1_poly_v_ploidy_mean_box)
dev.off()

## Variance ##
comp_res_list[ , tet_dif_dip_pco_1_var := 
                 comp_res_list$tet_pco_1_var - comp_res_list$dip_pco_1_var]
summary(comp_res_list$tet_dif_dip_pco_1_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.417e-05  3.281e-06  9.031e-06  8.858e-06  1.475e-05  4.068e-05
# variance in PCo1 value generally greater for tet than dip genotypes
summary(comp_res_list$tet_dif_dip_pco_1_var[cr_tet_inds])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -7.303e-06  5.623e-06  1.017e-05  1.078e-05  1.544e-05  4.068e-05
summary(comp_res_list$tet_dif_dip_pco_1_var[cr_oct_inds])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.417e-05 -6.629e-06  1.162e-06  2.395e-06  1.102e-05  3.480e-05
# variance in PCo1 generally greater for all samples, but bigger difference
#   for 4X than for 8X samples

tet_dip_var_1_melt <- data.table(tet_dip_dif = c(
  comp_res_list$tet_dif_dip_pco_1_var[cr_tet_inds],
  comp_res_list$tet_dif_dip_pco_1_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco1_tet_v_dip_var_density <- ggplot(tet_dip_var_1_melt) +
  geom_density(aes(x = tet_dip_dif, fill = comp), alpha = 0.5) +
#  coord_cartesian(xlim = c(0, 0.09)) +
  ggtitle(paste('Difference of variance in PCo1 values for\ntetrasomic vs ', 
                'disomic genotypes', sep = ''))

gg_pco1_tet_v_dip_var_box <- ggplot(tet_dip_var_1_melt) +
  geom_boxplot(aes(x = comp, y = tet_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(-3e-05, 4e-05)) +
    ggtitle(paste('Difference of variance in PCo1 values for\ntetrasomic vs ', 
                'disomic genotypes', sep = ''))

comp_res_list[ , poly_dif_dip_pco_1_var := 
                 comp_res_list$poly_pco_1_var - comp_res_list$dip_pco_1_var]
summary(comp_res_list$poly_dif_dip_pco_1_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.899e-05 -1.427e-06  3.135e-06  1.785e-06  5.598e-06  2.979e-05 
# Still see higher variance for polyploid genotypes, but to be expected since
#   see higher variance in 8X samples using tetrasomic genotypes
summary(comp_res_list$poly_dif_dip_pco_1_var[cr_tet_inds])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.609e-06  1.529e-06  3.838e-06  3.495e-06  5.888e-06  1.347e-05
summary(comp_res_list$poly_dif_dip_pco_1_var[cr_oct_inds])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.899e-05 -1.192e-05 -5.922e-06 -3.957e-06  3.344e-06  2.979e-05

poly_dip_var_1_melt <- data.table(poly_dip_dif = c(
  comp_res_list$poly_dif_dip_pco_1_var[cr_tet_inds],
  comp_res_list$poly_dif_dip_pco_1_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco1_poly_v_dip_var_density <- ggplot(poly_dip_var_1_melt) +
  geom_density(aes(x = poly_dip_dif, fill = comp), alpha = 0.5) +
  #  coord_cartesian(xlim = c(0, 0.09)) +
  ggtitle(paste('Difference of variance in PCo1 values for\npolyploid vs ', 
                'disomic genotypes', sep = ''))

gg_pco1_poly_v_dip_var_box <- ggplot(poly_dip_var_1_melt) +
  geom_boxplot(aes(x = comp, y = poly_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(-3e-05, 4e-05)) +
  ggtitle(paste('Difference of variance in PCo1 values for\npolyploid vs ', 
                'disomic genotypes', sep = ''))

pco_1_var_pdf_name <- paste(out_dir, 'natv2_pco_1_var_reproducibility.pdf',
                             sep = '') 
pdf(pco_1_var_pdf_name, height = 7, width = 10)
(gg_pco1_tet_v_dip_var_density + gg_pco1_tet_v_dip_var_box) /
  (gg_pco1_poly_v_dip_var_density + gg_pco1_poly_v_dip_var_box)
dev.off()


### PCo2 ###
dip_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
  pco_num = 2, control_name = 'J660.B')
comp_res_list[, dip_pco_2_mean := apply(dip_pco_2_stand_mat, 2, mean)]
comp_res_list[, dip_pco_2_var := apply(dip_pco_2_stand_mat, 2, var)]

tet_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
  pco_num = 2, control_name = 'J660.B')
comp_res_list[, tet_pco_2_mean := apply(tet_pco_2_stand_mat, 2, mean)]
comp_res_list[, tet_pco_2_var := apply(tet_pco_2_stand_mat, 2, var)]

poly_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
  pco_num = 2, control_name = 'J660.B')
comp_res_list[, poly_pco_2_mean := apply(poly_pco_2_stand_mat, 2, mean)]
comp_res_list[, poly_pco_2_var := apply(poly_pco_2_stand_mat, 2, var)]

## MEAN ##
comp_res_list[ , tet_dif_dip_pco_2_mean := 
                 comp_res_list$tet_pco_2_mean - comp_res_list$dip_pco_2_mean]
summary(abs(comp_res_list$tet_dif_dip_pco_2_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000112 0.0018358 0.0034843 0.0049691 0.0058001 0.0506270
summary(abs(comp_res_list$tet_dif_dip_pco_2_mean)[cr_tet_inds])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000208 0.0015074 0.0028658 0.0031238 0.0044047 0.0143409
summary(abs(comp_res_list$tet_dif_dip_pco_2_mean)[cr_oct_inds])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000112 0.0063457 0.0107519 0.0111656 0.0144338 0.0506270

tet_dip_mean_2_melt <- data.table(tet_dip_dif = c(
  abs(comp_res_list$tet_dif_dip_pco_2_mean[cr_tet_inds]),
  abs(comp_res_list$tet_dif_dip_pco_2_mean[cr_oct_inds])
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco2_tet_v_dip_mean_density <- ggplot(tet_dip_mean_2_melt) +
  geom_density(aes(x = tet_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.05)) +
  ggtitle(paste('Difference of mean PCo2 values for\ndisomic vs ', 
                'tetrasomic genotypes', sep = ''))

gg_pco2_tet_v_dip_mean_box <- ggplot(tet_dip_mean_2_melt) +
  geom_boxplot(aes(x = comp, y = tet_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(0, 0.05)) +
  ggtitle(paste('Difference of mean PCo2 values for\ndisomic vs ', 
                'tetrasomic genotypes', sep = ''))

comp_res_list[ , poly_dif_dip_pco_2_mean := 
                 comp_res_list$poly_pco_2_mean - comp_res_list$dip_pco_2_mean]
summary(abs(comp_res_list$poly_dif_dip_pco_2_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000005 0.0004565 0.0016477 0.0050767 0.0081374 0.0494012
summary(abs(comp_res_list$poly_dif_dip_pco_2_mean)[cr_tet_inds])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 5.030e-07 3.333e-04 9.169e-04 2.448e-03 2.905e-03 8.448e-03
summary(abs(comp_res_list$poly_dif_dip_pco_2_mean)[cr_oct_inds])
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000562 0.009457 0.013541 0.013904 0.017317 0.049401

poly_dip_mean_2_melt <- data.table(poly_dip_dif = c(
  abs(comp_res_list$poly_dif_dip_pco_2_mean[cr_tet_inds]),
  abs(comp_res_list$poly_dif_dip_pco_2_mean[cr_oct_inds])
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco2_poly_v_dip_mean_density <- ggplot(poly_dip_mean_2_melt) +
  geom_density(aes(x = poly_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.05)) +
  ggtitle(paste('Difference of mean PCo2 values for\ndisomic vs ', 
                'polyploid genotypes', sep = ''))

gg_pco2_poly_v_dip_mean_box <- ggplot(poly_dip_mean_2_melt) +
  geom_boxplot(aes(x = comp, y = poly_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(0, 0.05)) +
  ggtitle(paste('Difference of mean PCo2 values for\ndisomic vs ', 
                'polyploid genotypes', sep = ''))

comp_res_list[ , poly_dif_tet_pco_2_mean := 
                 comp_res_list$poly_pco_2_mean - comp_res_list$tet_pco_2_mean]
summary(abs(comp_res_list$poly_dif_tet_pco_2_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.116e-06 1.734e-03 3.404e-03 4.822e-03 6.013e-03 2.006e-02

poly_ploidy_mean_2_melt <- data.table(poly_ploidy_dif = c(
  abs(comp_res_list$poly_dif_dip_pco_2_mean[cr_tet_inds]),
  abs(comp_res_list$poly_dif_tet_pco_2_mean[cr_oct_inds])
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco2_poly_v_ploidy_mean_density <- ggplot(poly_ploidy_mean_2_melt) +
  geom_density(aes(x = poly_ploidy_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.05)) + 
  ggtitle(paste('Difference of mean PCo2 values for\npolyploid vs ', 
                'ploidy-appropriate genotypes', sep = ''))

gg_pco2_poly_v_ploidy_mean_box <- ggplot(poly_ploidy_mean_2_melt) +
  geom_boxplot(aes(x = comp, y = poly_ploidy_dif, fill = comp)) +
  coord_cartesian(ylim = c(0, 0.05)) +
  ggtitle(paste('Difference of mean PCo1 values for\npolyploid vs ', 
                'ploidy-appropriate genotypes', sep = ''))

pco_2_mean_pdf_name <- paste(out_dir, 'natv2_pco_2_mean_reproducibility.pdf',
                             sep = '') 
pdf(pco_2_mean_pdf_name, height = 10, width = 10)
(gg_pco2_tet_v_dip_mean_density + gg_pco2_tet_v_dip_mean_box) /
  (gg_pco2_poly_v_dip_mean_density + gg_pco2_poly_v_dip_mean_box) /
  (gg_pco2_poly_v_ploidy_mean_density + gg_pco2_poly_v_ploidy_mean_box)
dev.off()

## VARIANCE ##
comp_res_list[ , tet_dif_dip_pco_2_var := 
                 comp_res_list$tet_pco_2_var - comp_res_list$dip_pco_2_var]
summary(comp_res_list$tet_dif_dip_pco_2_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.073e-05 -6.910e-07  1.012e-05  8.203e-06  1.931e-05  5.832e-05
summary(comp_res_list$tet_dif_dip_pco_2_var[cr_tet_inds])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -3.048e-05  5.756e-06  1.329e-05  1.305e-05  2.058e-05  5.193e-05
summary(comp_res_list$tet_dif_dip_pco_2_var[cr_oct_inds])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.073e-05 -1.869e-05 -8.370e-06 -8.064e-06  9.468e-07  5.832e-05
# On average, variance is higher in 4X and lower in 8X when using 
#  tetrasomic genotypes

tet_dip_var_2_melt <- data.table(tet_dip_dif = c(
  comp_res_list$tet_dif_dip_pco_2_var[cr_tet_inds],
  comp_res_list$tet_dif_dip_pco_2_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco2_tet_v_dip_var_density <- ggplot(tet_dip_var_2_melt) +
  geom_density(aes(x = tet_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(-6e-05, 6e-05)) +
  ggtitle(paste('Difference of variance in PCo2 values for\ntetrasomic vs ', 
                'disomic genotypes', sep = ''))

gg_pco2_tet_v_dip_var_box <- ggplot(tet_dip_var_2_melt) +
  geom_boxplot(aes(x = comp, y = tet_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(-6e-05, 6e-05)) +
  ggtitle(paste('Difference of variance in PCo2 values for\ntetrasomic vs ', 
                'disomic genotypes', sep = ''))

comp_res_list[ , poly_dif_dip_pco_2_var := 
                 comp_res_list$poly_pco_2_var - comp_res_list$dip_pco_2_var]
summary(comp_res_list$poly_dif_dip_pco_2_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.921e-05 -2.505e-06  1.577e-06 -3.664e-07  3.552e-06  3.317e-05 
# Variance is comparable to disomic when  using polyploid genotypes

poly_dip_var_2_melt <- data.table(poly_dip_dif = c(
  comp_res_list$poly_dif_dip_pco_2_var[cr_tet_inds],
  comp_res_list$poly_dif_dip_pco_2_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco2_poly_v_dip_var_density <- ggplot(poly_dip_var_2_melt) +
  geom_density(aes(x = poly_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(-6e-05, 6e-05)) +
  ggtitle(paste('Difference of variance in PCo2 values for\npolyploid vs ', 
                'disomic genotypes', sep = ''))

gg_pco2_poly_v_dip_var_box <- ggplot(poly_dip_var_2_melt) +
  geom_boxplot(aes(x = comp, y = poly_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(-6e-05, 6e-05)) +
  ggtitle(paste('Difference of variance in PCo2 values for\npolyploid vs ', 
                'disomic genotypes', sep = ''))

comp_res_list[ , poly_dif_tet_pco_2_var := 
                 comp_res_list$poly_pco_2_var - comp_res_list$tet_pco_2_var]
summary(comp_res_list$poly_dif_tet_pco_2_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -4.942e-05 -1.514e-05 -7.494e-06 -8.569e-06 -1.762e-06  2.816e-05

poly_ploidy_var_2_melt <- data.table(poly_ploidy_dif = c(
  comp_res_list$poly_dif_dip_pco_2_var[cr_tet_inds],
  comp_res_list$poly_dif_tet_pco_2_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco2_poly_v_ploidy_var_density <- ggplot(poly_ploidy_var_2_melt) +
  geom_density(aes(x = poly_ploidy_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(-6e-05, 6e-05)) +
  ggtitle(paste('Difference of variance in PCo2 values for\npolyploid vs ', 
                'ploidy-specific genotypes', sep = ''))

gg_pco2_poly_v_ploidy_var_box <- ggplot(poly_ploidy_var_2_melt) +
  geom_boxplot(aes(x = comp, y = poly_ploidy_dif, fill = comp)) +
  coord_cartesian(ylim = c(-6e-05, 6e-05)) +
  ggtitle(paste('Difference of variance in PCo2 values for\npolyploid vs ', 
                'ploidy-specific genotypes', sep = ''))

pco_2_var_pdf_name <- paste(out_dir, 'natv2_pco_2_var_reproducibility.pdf',
                            sep = '') 
pdf(pco_2_var_pdf_name, height = 7, width = 10)
(gg_pco2_tet_v_dip_var_density + gg_pco2_tet_v_dip_var_box) /
  (gg_pco2_poly_v_dip_var_density + gg_pco2_poly_v_dip_var_box) /
  (gg_pco2_poly_v_ploidy_var_density + gg_pco2_poly_v_ploidy_var_box)
dev.off()

## PCo3
dip_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
  pco_num = 3, control_name = 'J456.B')
comp_res_list[, dip_pco_3_mean := apply(dip_pco_3_stand_mat, 2, mean)]
comp_res_list[, dip_pco_3_var := apply(dip_pco_3_stand_mat, 2, var)]

tet_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
  pco_num = 3, control_name = 'J456.B')
comp_res_list[, tet_pco_3_mean := apply(tet_pco_3_stand_mat, 2, mean)]
comp_res_list[, tet_pco_3_var := apply(tet_pco_3_stand_mat, 2, var)]

poly_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
  pco_num = 3, control_name = 'J456.B')
comp_res_list[, poly_pco_3_mean := apply(poly_pco_3_stand_mat, 2, mean)]
comp_res_list[, poly_pco_3_var := apply(poly_pco_3_stand_mat, 2, var)]

comp_res_list[ , tet_dif_dip_pco_3_mean := 
                 comp_res_list$dip_pco_3_mean - comp_res_list$tet_pco_3_mean]
comp_res_list[ , poly_dif_dip_pco_3_mean := 
                 comp_res_list$poly_pco_3_mean - comp_res_list$dip_pco_3_mean]
comp_res_list[ , poly_dif_tet_pco_3_mean := 
                 comp_res_list$poly_pco_3_mean - comp_res_list$tet_pco_3_mean]

tet_dip_mean_3_melt <- data.table(tet_dip_dif = c(
  abs(comp_res_list$tet_dif_dip_pco_3_mean[cr_tet_inds]),
  abs(comp_res_list$tet_dif_dip_pco_3_mean[cr_oct_inds])
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco3_tet_v_dip_mean_density <- ggplot(tet_dip_mean_3_melt) +
  geom_density(aes(x = tet_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.06)) +
  ggtitle(paste('Difference of mean PCo3 values for\ndisomic vs ', 
                'tetrasomic genotypes', sep = ''))

gg_pco3_tet_v_dip_mean_box <- ggplot(tet_dip_mean_3_melt) +
  geom_boxplot(aes(x = comp, y = tet_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(0, 0.06)) +
  ggtitle(paste('Difference of mean PCo3 values for\ndisomic vs ', 
                'tetrasomic genotypes', sep = ''))

poly_tet_mean_3_melt <- data.table(poly_tet_dif = c(
  abs(comp_res_list$poly_dif_tet_pco_3_mean[cr_tet_inds]),
  abs(comp_res_list$poly_dif_tet_pco_3_mean[cr_oct_inds])
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco3_poly_v_tet_mean_density <- ggplot(poly_tet_mean_3_melt) +
  geom_density(aes(x = poly_tet_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.06)) +
  ggtitle(paste('Difference of mean PCo3 values for\npolyploid vs ', 
                'tetrasomic genotypes', sep = ''))

gg_pco3_poly_v_tet_mean_box <- ggplot(poly_tet_mean_3_melt) +
  geom_boxplot(aes(x = comp, y = poly_tet_dif, fill = comp)) +
  coord_cartesian(ylim = c(0, 0.06)) +
  ggtitle(paste('Difference of mean PCo3 values for\npolyploid vs ', 
                'tetrasomic genotypes', sep = ''))

poly_ploidy_mean_3_melt <- data.table(poly_ploidy_dif = c(
  abs(comp_res_list$poly_dif_dip_pco_3_mean[cr_tet_inds]),
  abs(comp_res_list$poly_dif_tet_pco_3_mean[cr_oct_inds])
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco3_poly_v_ploidy_mean_density <- ggplot(poly_ploidy_mean_3_melt) +
  geom_density(aes(x = poly_ploidy_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.06)) +
  ggtitle(paste('Difference of mean PCo3 values for\npolyploid vs ', 
                'ploidy-appropriate genotypes', sep = ''))

gg_pco3_poly_v_ploidy_mean_box <- ggplot(poly_ploidy_mean_3_melt) +
  geom_boxplot(aes(x = comp, y = poly_ploidy_dif, fill = comp)) +
  coord_cartesian(ylim = c(0, 0.06)) +
  ggtitle(paste('Difference of mean PCo3 values for\npolyploid vs ', 
                'ploidy-appropriate genotypes', sep = ''))

pco_3_mean_pdf_name <- paste(out_dir, 'natv2_pco_3_mean_reproducibility.pdf',
                             sep = '') 
pdf(pco_3_mean_pdf_name, height = 10, width = 10)
(gg_pco3_tet_v_dip_mean_density + gg_pco3_tet_v_dip_mean_box) /
  (gg_pco3_poly_v_tet_mean_density + gg_pco3_poly_v_tet_mean_box) /
  (gg_pco3_poly_v_ploidy_mean_density + gg_pco3_poly_v_ploidy_mean_box)
dev.off()

### VARIANCE ###
comp_res_list[ , tet_dif_dip_pco_3_var := 
                 comp_res_list$tet_pco_3_var - comp_res_list$dip_pco_3_var]

tet_dip_var_3_melt <- data.table(tet_dip_dif = c(
  comp_res_list$tet_dif_dip_pco_3_var[cr_tet_inds],
  comp_res_list$tet_dif_dip_pco_3_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco3_tet_v_dip_var_density <- ggplot(tet_dip_var_3_melt) +
  geom_density(aes(x = tet_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(-4e-04, 4e-04)) +
  ggtitle(paste('Difference of variance in PCo3 values for\ntetrasomic vs ', 
                'disomic genotypes', sep = ''))

gg_pco3_tet_v_dip_var_box <- ggplot(tet_dip_var_3_melt) +
  geom_boxplot(aes(x = comp, y = tet_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(-4e-04, 4e-04)) +
  ggtitle(paste('Difference of variance in PCo3 values for\ntetrasomic vs ', 
                'disomic genotypes', sep = ''))

comp_res_list[ , poly_dif_dip_pco_3_var := 
                 comp_res_list$poly_pco_3_var - comp_res_list$dip_pco_3_var]

poly_dip_var_3_melt <- data.table(poly_dip_dif = c(
  comp_res_list$poly_dif_dip_pco_3_var[cr_tet_inds],
  comp_res_list$poly_dif_dip_pco_3_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco3_poly_v_dip_var_density <- ggplot(poly_dip_var_3_melt) +
  geom_density(aes(x = poly_dip_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(-4e-04, 4e-04)) +
  ggtitle(paste('Difference of variance in PCo3 values for\npolyploid vs ', 
                'disomic genotypes', sep = ''))

gg_pco3_poly_v_dip_var_box <- ggplot(poly_dip_var_3_melt) +
  geom_boxplot(aes(x = comp, y = poly_dip_dif, fill = comp)) +
  coord_cartesian(ylim = c(-4e-04, 4e-04)) +
  ggtitle(paste('Difference of variance in PCo3 values for\npolyploid vs ', 
                'disomic genotypes', sep = ''))

comp_res_list[ , poly_dif_tet_pco_3_var := 
                 comp_res_list$poly_pco_3_var - comp_res_list$tet_pco_3_var]

poly_ploidy_var_3_melt <- data.table(poly_ploidy_dif = c(
  comp_res_list$poly_dif_dip_pco_3_var[cr_tet_inds],
  comp_res_list$poly_dif_tet_pco_3_var[cr_oct_inds]
),
comp = c(
  rep('4X', times = length(cr_tet_inds)),
  rep('8X', times = length(cr_oct_inds))
)
)

gg_pco3_poly_v_ploidy_var_density <- ggplot(poly_ploidy_var_3_melt) +
  geom_density(aes(x = poly_ploidy_dif, fill = comp), alpha = 0.5) +
  coord_cartesian(xlim = c(-4e-04, 4e-04)) +
  ggtitle(paste('Difference of variance in PCo3 values for\npolyploid vs ', 
                'ploidy-specific genotypes', sep = ''))

gg_pco3_poly_v_ploidy_var_box <- ggplot(poly_ploidy_var_3_melt) +
  geom_boxplot(aes(x = comp, y = poly_ploidy_dif, fill = comp)) +
  coord_cartesian(ylim = c(-4e-04, 4e-04)) +
  ggtitle(paste('Difference of variance in PCo3 values for\npolyploid vs ', 
                'ploidy-specific genotypes', sep = ''))

pco_3_var_pdf_name <- paste(out_dir, 'natv2_pco_3_var_reproducibility.pdf',
                            sep = '') 
pdf(pco_3_var_pdf_name, height = 10, width = 10)
(gg_pco3_tet_v_dip_var_density + gg_pco3_tet_v_dip_var_box) /
  (gg_pco3_poly_v_dip_var_density + gg_pco3_poly_v_dip_var_box) /
  (gg_pco3_poly_v_ploidy_var_density + gg_pco3_poly_v_ploidy_var_box)
dev.off()



quit(save = 'no')

