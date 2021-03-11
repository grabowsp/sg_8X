# Analysis of tree reproducibility for different genotypes

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
#library(ape, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
#library(phangorn, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
library(ggplot2)

### INPUT DATA ###
dist_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist/'

dip_files <- system(paste('ls ', dist_dir, '*diploid_DistMat.rds', sep = ''), 
  intern = T)
tet_files <- system(paste('ls ', dist_dir, '*tetraploid_DistMat.rds', 
  sep = ''), intern = T)
poly_files <- system(paste('ls ', dist_dir, '*polyploid_DistMat.rds', 
  sep = ''), intern = T)

dip_list <- list()
for(dip in seq(length(dip_files))){
  tmp_data <- readRDS(dip_files[dip])
  tmp_mat <- as.matrix(tmp_data[['euclidean_dist']])
  rm_inds <- grep('TET', rownames(tmp_mat))
  tmp_dist <- as.dist(tmp_mat[-rm_inds, -rm_inds])
  dip_list[[dip]] <- tmp_dist
}

tet_list <- list()
for(tet in seq(length(tet_files))){
  tmp_data <- readRDS(tet_files[tet])
  tmp_mat <- as.matrix(tmp_data[['euclidean_dist']])
  rm_inds <- grep('TET', rownames(tmp_mat))
  tmp_dist <- as.dist(tmp_mat[-rm_inds, -rm_inds])
  tet_list[[tet]] <- tmp_dist
}

poly_list <- list()
for(pf in seq(length(poly_files))){
  tmp_data <- readRDS(poly_files[pf])
  tmp_mat <- as.matrix(tmp_data[['euclidean_dist']])
  rm_inds <- grep('TET', rownames(tmp_mat))
  tmp_dist <- as.dist(tmp_mat[-rm_inds, -rm_inds])
  poly_list[[pf]] <- tmp_dist
}

samp_meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

### SET OUTPUT ###
res_file <- paste(dist_dir, 'natv2.100k.NJtree.reproducibility.results.rds',
  sep = '')

#dip_tree_file <- paste(dist_dir, 'natv2.100k.dip_NJtree.list.v1.rds', sep = '')
#tet_tree_file <- paste(dist_dir, 'natv2.100k.tet_NJtree.list.v1.rds', sep = '')
#poly_tree_file <- paste(dist_dir, 'natv2.100k.poly_NJtree.list.v1.rds', 
#  sep = '')

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

dip_pcoa_list <- lapply(dip_list, function(x) cmdscale(x, k = 20))
tet_pcoa_list <- lapply(tet_list, function(x) cmdscale(x, k = 20))
poly_pcoa_list <- lapply(poly_list, function(x) cmdscale(x, k = 20))

comp_res_list <- data.table(samp_name = rownames(dip_pcoa_list[[1]]))
comp_res_list[, ploidy := as.character(NA)]
comp_res_list[which(comp_res_list$samp_name %in% tet_names), ploidy := '4X']
comp_res_list[which(comp_res_list$samp_name %in% oct_names), ploidy := '8X']

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
summary(abs((tet_pco_1_mean - dip_pco_1_mean))[which(names(dip_pco_1_mean)
  %in% tet_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.243e-06 8.232e-04 1.604e-03 2.685e-03 4.037e-03 2.216e-02
summary(abs((tet_pco_1_mean - dip_pco_1_mean))[which(names(dip_pco_1_mean)
  %in% oct_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3.339e-05 1.296e-02 1.885e-02 1.857e-02 2.305e-02 7.689e-02
# Generally small difference between means in 4X individuals when using
#   tetrasomic vs disomc genotypes; difference an order of magnitude 
#   bigger for 8X
gg_pco1_tet_v_dip_mean <- ggplot(comp_res_list) + 
  geom_density(aes(x = abs(tet_dif_dip_pco_1_mean)), fill = 'grey50') +
  geom_density(aes(x = 
      abs(tet_dif_dip_pco_1_mean[which(comp_res_list$ploidy == '4X')])), 
    fill = 'red') +
  geom_density(aes(x = 
      abs(tet_dif_dip_pco_1_mean[which(comp_res_list$ploidy == '8X')])),
    fill = 'blue')

# CONTINUE ADAPTING FROM HERE #


summary(abs(poly_pco_1_mean - dip_pco_1_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.240e-06 5.772e-04 1.254e-03 6.067e-03 4.159e-03 7.264e-02
summary(abs(poly_pco_1_mean - dip_pco_1_mean)[which(names(dip_pco_1_mean)
  %in% tet_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.236e-06 4.585e-04 9.293e-04 1.217e-03 1.590e-03 5.697e-03
summary(abs(poly_pco_1_mean - tet_pco_1_mean)[which(names(dip_pco_1_mean)
  %in% oct_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 3.261e-05 4.391e-03 7.149e-03 6.230e-03 7.740e-03 1.186e-02
# Mean values within ploidy level similar to for polyploid results to  when 
#   using ploidy-specific genotypes   

## Variance ##
summary(tet_pco_1_var - dip_pco_1_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.417e-05  3.281e-06  9.031e-06  8.858e-06  1.475e-05  4.068e-05
# variance in PCo1 value generally greater for tet than dip genotypes
summary((tet_pco_1_var - dip_pco_1_var)[which(names(dip_pco_1_var) 
  %in% tet_names)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -7.303e-06  5.623e-06  1.017e-05  1.078e-05  1.544e-05  4.068e-05
summary((tet_pco_1_var - dip_pco_1_var)[which(names(dip_pco_1_var) 
  %in% oct_names)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.417e-05 -6.629e-06  1.162e-06  2.395e-06  1.102e-05  3.480e-05
# variance in PCo1 generally greater for all samples, but bigger difference
#   for 4X than for 8X samples

summary(poly_pco_1_var - dip_pco_1_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.899e-05 -1.427e-06  3.135e-06  1.785e-06  5.598e-06  2.979e-05 
# Still see higher variance for polyploid genotypes, but to be expected since
#   see higher variance in 8X samples using tetrasomic genotypes

### PCo2 ###

dip_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
  pco_num = 2, control_name = 'J660.B')
dip_pco_2_var <- apply(dip_pco_2_stand_mat, 2, var)
dip_pco_2_mean <- apply(dip_pco_2_stand_mat, 2, mean)

tet_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
  pco_num = 2, control_name = 'J660.B')
tet_pco_2_var <- apply(tet_pco_2_stand_mat, 2, var)
tet_pco_2_mean <- apply(tet_pco_2_stand_mat, 2, mean)

poly_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
  pco_num = 2, control_name = 'J660.B')
poly_pco_2_var <- apply(poly_pco_2_stand_mat, 2, var)
poly_pco_2_mean <- apply(poly_pco_2_stand_mat, 2, mean)

## MEAN ##
summary(abs(tet_pco_2_mean - dip_pco_2_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000112 0.0018358 0.0034843 0.0049691 0.0058001 0.0506270
summary(abs((tet_pco_2_mean - dip_pco_2_mean))[which(names(dip_pco_1_mean)
  %in% tet_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000208 0.0015074 0.0028658 0.0031238 0.0044047 0.0143409
summary(abs((tet_pco_2_mean - dip_pco_2_mean))[which(names(dip_pco_1_mean)
  %in% oct_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000112 0.0063457 0.0107519 0.0111656 0.0144338 0.0506270

summary(abs(poly_pco_2_mean - dip_pco_2_mean))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000005 0.0004565 0.0016477 0.0050767 0.0081374 0.0494012
summary(abs(poly_pco_2_mean - dip_pco_2_mean)[which(names(dip_pco_1_mean)
  %in% tet_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 5.030e-07 3.333e-04 9.169e-04 2.448e-03 2.905e-03 8.448e-03
summary(abs(poly_pco_2_mean - tet_pco_2_mean)[which(names(dip_pco_1_mean)
  %in% oct_names)])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 7.579e-05 3.351e-03 4.049e-03 5.283e-03 4.943e-03 2.006e-02

## VARIANCE ##
summary(tet_pco_2_var - dip_pco_2_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.073e-05 -6.910e-07  1.012e-05  8.203e-06  1.931e-05  5.832e-05
summary((tet_pco_2_var - dip_pco_2_var)[which(names(dip_pco_1_var)
  %in% tet_names)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -3.048e-05  5.756e-06  1.329e-05  1.305e-05  2.058e-05  5.193e-05
summary((tet_pco_2_var - dip_pco_2_var)[which(names(dip_pco_1_var)
  %in% oct_names)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.073e-05 -1.869e-05 -8.370e-06 -8.064e-06  9.468e-07  5.832e-05
# On average, variance is higher in 4X and lower in 8X when using 
#  tetrasomic genotypes

summary(poly_pco_2_var - dip_pco_2_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -5.921e-05 -2.505e-06  1.577e-06 -3.664e-07  3.552e-06  3.317e-05 
# Variance is comparable to disomic when  using polyploid genotypes

## PCo3

dip_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
  pco_num = 3, control_name = 'J456.B')
dip_pco_3_var <- apply(dip_pco_3_stand_mat, 2, var)

tet_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
  pco_num = 3, control_name = 'J456.B')
tet_pco_3_var <- apply(tet_pco_3_stand_mat, 2, var)

poly_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
  pco_num = 3, control_name = 'J456.B')
poly_pco_3_var <- apply(poly_pco_3_stand_mat, 2, var)

summary(tet_pco_3_var - dip_pco_3_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -3.528e-04 -5.267e-05  2.728e-05  1.849e-05  9.109e-05  3.680e-04
summary((tet_pco_3_var - dip_pco_3_var)[which(names(dip_pco_1_var)
  %in% tet_names)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -1.397e-04  8.623e-06  4.249e-05  5.902e-05  1.199e-04  3.680e-04
summary((tet_pco_3_var - dip_pco_3_var)[which(names(dip_pco_1_var)
  %in% oct_names)])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -3.528e-04 -2.219e-04 -1.014e-04 -1.176e-04 -3.183e-05  2.803e-04
# On average, variance is higher in 4X and lower in 8X when using tetrasomic
#  genotypes

summary(poly_pco_3_var - dip_pco_3_var)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -2.934e-04 -4.891e-05  8.540e-07 -8.854e-06  2.797e-05  1.942e-04
# Variance is comparable to disomic when using polyploid genotypes







quit(save = 'no')

