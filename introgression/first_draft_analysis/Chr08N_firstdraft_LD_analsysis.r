# Analysis of LD results for first draft of manuscript

# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/Chr08N_ATLvsMW_Fst0.5_results/'

# LD result files
atl_rda_ld_file <- 'ATL_Chr08N_ATLvsMW_Fst0.5_RDAonly.ld'
gulf_rda_ld_file <- 'GULF_Chr08N_ATLvsMW_Fst0.5_RDAonly.ld'

atl_rda_ld_mat <- as.matrix(fread(paste(res_dir, atl_rda_ld_file,
  sep = ''), header = F))
gulf_rda_ld_mat <- as.matrix(fread(paste(res_dir, gulf_rda_ld_file,
  sep = ''), header = F))

mwall_rda_ld_file <- 'MWall_Chr08N_ATLvsMW_Fst0.5_RDAonly.ld'
mwall_control_ld_file <- 'MWall_Chr08N_ATLvsMW_ControlSNPs.ld'

mw4X_rda_ld_file <- 'MW4X_NoIntrogression_Chr08N_ATLvsMW_Fst0.5_RDAonly.ld'
mw8X_rda_ld_file <- 'MW8X_Chr08N_ATLvsMW_Fst0.5_RDAonly.ld'

mwall_rda_ld_mat <- as.matrix(fread(paste(res_dir, mwall_rda_ld_file,
  sep = ''), header = F))

# SNP files
rda_snp_file <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_ATL_into_MW_Fst_0.5.RDA_SNP_info.txt'

rda_snps <- fread(rda_snp_file)
control_snp_file <- '/home/f1p1/tmp/switchgrass_8X/MWall_control_snps/ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.txt'
control_snps <- fread(control_snp_file)

# RDA result file
rda_res_file <- paste('/home/f1p1/tmp/switchgrass_8X/firstdraft_fval_files/',
  'MW8X_MWvATL_Fst0.5_RDA_Fval_overlap_RDASNPsonly.txt', sep = '')
rda_res <- fread(rda_res_file)

### SET VARIABLES ###

dist_inds <- c(0, 1e4, 1e5, 1e6, 1e7, 1e8)

#####

chr08N_rda_snps <- rda_snps[CHROM == 'Chr08N', ]

# generate distance matrix between snps
chr08N_rda_dist <- as.matrix(dist(chr08N_rda_snps$POS, method = 'manhattan', 
  upper = T, diag = T))

# the "overlap hit" SNPs
overlap_inds <- intersect(which(rda_res$PLOT_CLASS == 'OVERLAP'),
  grep('Chr08N', rda_res$SNP_NAME))
overlap_names <- rda_res[overlap_inds, SNP_NAME]
hit_inds <- which(chr08N_rda_snps$ID %in% overlap_names)

# chromosome before the hits
start_lims <- c(1, 21300000)
start_inds <- which(chr08N_rda_snps$POS > start_lims[1] &
  chr08N_rda_snps$POS < start_lims[2])

# regions with hits
big_locus_lims <- c(21300000, 37300000)
big_locus_inds <- which(chr08N_rda_snps$POS > big_locus_lims[1] &
  chr08N_rda_snps$POS < big_locus_lims[2])

# chromosome after the hits
end_lims <- c(37300000, 1e8)
end_inds <- which(chr08N_rda_snps$POS > end_lims[1] &
  chr08N_rda_snps$POS < end_lims[2])

###

# MWall values

test_ld_mat <- mwall_rda_ld_mat
test_dist_mat <- chr08N_rda_dist

ld_decay_list <- list()
for(j in seq(nrow(test_ld_mat))){
  test_samp <- j
  tmp_res_vec <- c()

for(i in seq(length(dist_inds)-1)){
 tmp_res <- mean(test_ld_mat[which(test_dist_mat[, test_samp] > dist_inds[i] &
  test_dist_mat[, test_samp] <= dist_inds[(i+1)]) , test_samp], na.rm = T)
 tmp_res_vec <- c(tmp_res_vec, tmp_res)
}
  ld_decay_list[[j]] <- tmp_res_vec
}

mw_decay_mat <- matrix(data = unlist(ld_decay_list),
  nrow = length(ld_decay_list), byrow = T)

mw_start_decay <- apply(mw_decay_mat[start_inds, ], 2, mean, na.rm = T)
mw_locus_decay <- apply(mw_decay_mat[big_locus_inds, ], 2, mean, na.rm = T)
mw_end_decay <- apply(mw_decay_mat[end_inds, ], 2, mean, na.rm = T)

# Gulf values

test_ld_mat <- gulf_rda_ld_mat
test_dist_mat <- chr08N_rda_dist

ld_decay_list <- list()
for(j in seq(nrow(test_ld_mat))){
  test_samp <- j
  tmp_res_vec <- c()

for(i in seq(length(dist_inds)-1)){
 tmp_res <- mean(test_ld_mat[which(test_dist_mat[, test_samp] > dist_inds[i] &
  test_dist_mat[, test_samp] <= dist_inds[(i+1)]) , test_samp], na.rm = T)
 tmp_res_vec <- c(tmp_res_vec, tmp_res)
}
  ld_decay_list[[j]] <- tmp_res_vec
}

gulf_decay_mat <- matrix(data = unlist(ld_decay_list),
  nrow = length(ld_decay_list), byrow = T)

gulf_start_decay <- apply(gulf_decay_mat[start_inds, ], 2, mean, na.rm = T)
gulf_locus_decay <- apply(gulf_decay_mat[big_locus_inds, ], 2, mean, na.rm = T)
gulf_end_decay <- apply(gulf_decay_mat[end_inds, ], 2, mean, na.rm = T)


# Atl results
test_ld_mat <- atl_rda_ld_mat
test_dist_mat <- chr08N_rda_dist

ld_decay_list <- list()
for(j in seq(nrow(test_ld_mat))){
  test_samp <- j
  tmp_res_vec <- c()

for(i in seq(length(dist_inds)-1)){
 tmp_res <- mean(test_ld_mat[which(test_dist_mat[, test_samp] > dist_inds[i] &
  test_dist_mat[, test_samp] <= dist_inds[(i+1)]) , test_samp], na.rm = T)
 tmp_res_vec <- c(tmp_res_vec, tmp_res)
}
  ld_decay_list[[j]] <- tmp_res_vec
}

atl_decay_mat <- matrix(data = unlist(ld_decay_list),
  nrow = length(ld_decay_list), byrow = T)

atl_start_decay <- apply(atl_decay_mat[start_inds, ], 2, mean, na.rm = T)
atl_locus_decay <- apply(atl_decay_mat[big_locus_inds, ], 2, mean, na.rm = T)
atl_end_decay <- apply(atl_decay_mat[end_inds, ], 2, mean, na.rm = T)

# Comparisons

mw_locus_decay[2]/(mean(c(mw_start_decay[2], mw_end_decay[2])))
# 3.74512
atl_locus_decay[2]/(mean(c(atl_start_decay[2], atl_end_decay[2])))
# 2.51928
gulf_locus_decay[2]/(mean(c(gulf_start_decay[2], gulf_end_decay[2])))
# 2.608377

mw_locus_decay[4]/(mean(c(mw_start_decay[4], mw_end_decay[4])))
# 5.275808
atl_locus_decay[4]/(mean(c(atl_start_decay[4], atl_end_decay[4])))
# 2.171906
gulf_locus_decay[4]/(mean(c(gulf_start_decay[4], gulf_end_decay[4])))
# 1.724465


############






start_inds <- which(chr08N_rda_snps$POS > start_lims[1] &
  chr08N_rda_snps$POS < start_lims[2])

# 

locus_1_lims <- c(21300000,22000000)
locus_1_inds <- which(chr08N_rda_snps$POS > locus_1_lims[1] &
  chr08N_rda_snps$POS < locus_1_lims[2])

inter_1_lims <- c(22000000, 23300000)
inter_1_inds <- which(chr08N_rda_snps$POS > inter_1_lims[1] &
  chr08N_rda_snps$POS < inter_1_lims[2])

locus_2_lims <- c(23300000, 25100000)
locus_2_inds <- which(chr08N_rda_snps$POS > locus_2_lims[1] &
  chr08N_rda_snps$POS < locus_2_lims[2])

inter_2_lims <- c(25100000, 27000000)
inter_2_inds <- which(chr08N_rda_snps$POS > inter_2_lims[1] &
  chr08N_rda_snps$POS < inter_2_lims[2])

locus_3_lims <- c(27000000, 32000000)
locus_3_inds <- which(chr08N_rda_snps$POS > locus_3_lims[1] &
  chr08N_rda_snps$POS < locus_3_lims[2])

inter_3_lims <- c(32000000, 34000000)
inter_3_inds <- which(chr08N_rda_snps$POS > inter_3_lims[1] &
  chr08N_rda_snps$POS < inter_3_lims[2])

locus_4_lims <- c(34000000, 37300000)
locus_4_inds <- which(chr08N_rda_snps$POS > locus_4_lims[1] &
  chr08N_rda_snps$POS < locus_4_lims[2])

end_lims <- c(37300000, 1e8)
end_inds <- which(chr08N_rda_snps$POS > end_lims[1] &
  chr08N_rda_snps$POS < end_lims[2])

big_locus_lims <- c(21300000, 37300000)
big_locus_inds <- which(chr08N_rda_snps$POS > big_locus_lims[1] &
  chr08N_rda_snps$POS < big_locus_lims[2])

test_ld_mat <- mwall_rda_ld_mat
test_dist_mat <- chr08N_rda_dist
test_samp <- 1

ld_decay_list <- list()
for(j in seq(nrow(test_ld_mat))){
  test_samp <- j
  tmp_res_vec <- c()

for(i in seq(length(dist_inds)-1)){
 tmp_res <- mean(test_ld_mat[which(test_dist_mat[, test_samp] > dist_inds[i] &
  test_dist_mat[, test_samp] <= dist_inds[(i+1)]) , test_samp])
 tmp_res_vec <- c(tmp_res_vec, tmp_res)
}
  ld_decay_list[[j]] <- tmp_res_vec
}

ld_decay_mat <- matrix(data = unlist(ld_decay_list), 
  nrow = length(ld_decay_list), byrow = T)

apply(ld_decay_mat[-hit_inds, ], 2, mean, na.rm = T)

ld_decay_mat[inter_1_inds, ]

start_ld_decay <- apply(ld_decay_mat[start_inds, ], 2, mean, na.rm = T)

locus_1_ld_decay <- apply(ld_decay_mat[locus_1_inds, ], 2, mean, na.rm = T)

inter_1_ld_decay <- apply(ld_decay_mat[inter_1_inds, ], 2, mean, na.rm = T)

locus_2_ld_decay <- apply(ld_decay_mat[locus_2_inds, ], 2, mean, na.rm = T)

inter_2_ld_decay <- apply(ld_decay_mat[inter_2_inds, ], 2, mean, na.rm = T)

locus_3_ld_decay <- apply(ld_decay_mat[locus_3_inds, ], 2, mean, na.rm = T)

inter_3_ld_decay <- apply(ld_decay_mat[inter_3_inds, ], 2, mean, na.rm = T)

locus_4_ld_decay <- apply(ld_decay_mat[locus_4_inds, ], 2, mean, na.rm = T)

end_ld_decay <- apply(ld_decay_mat[end_inds, ], 2, mean, na.rm = T)

big_locus_ld_decay <- apply(ld_decay_mat[big_locus_inds, ], 2, mean, na.rm = T)


