# Generate trees of subsets of samples to better select comparison samples
#   for 4X/8X kmer analysis

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
library(ape, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
library(phangorn, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')

### INPUT DATA ###
meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

dist_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/GW_all_samps_ld0.3_symmetric.mdist'
dist_res <- fread(dist_res_file, header = F)
dist_ids <- fread(paste(dist_res_file, 'id', sep = '.'), header = F)
dist_mat <- as.matrix(dist_res)
colnames(dist_mat) <- rownames(dist_mat) <- dist_ids$V1

admix_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt'
admix_res <- fread(admix_res_file)
# 1 = Atlantic, 2 = Gulf, 3 = Midwest

ni_admix_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_NorthInland_admix/GW_50k_geobigNorthInland.3.results.txt'
ni_res_0 <- fread(ni_admix_file)
# 1 = 8X West, 2 = 8X East, 3 = 4X

natv2_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_file, header = F)$V1

#

test_samp_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/'

atlantic_file <- paste(test_samp_dir, 'Atlantic_ref_samps_40.txt', sep = '')
atl_samps <- fread(atlantic_file, header = T)[, ID]

gulf_file <- paste(test_samp_dir, 'Gulf_ref_samps_40.txt', sep = '')
gulf_samps <- fread(gulf_file, header = T)[, ID]

mw_file <- paste(test_samp_dir, 'Midwest_ref_samps_40.txt', sep = '')
mw_samps <- fread(mw_file, header = T)[, ID]

candidate_samp_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/Introgression_target_individuals.txt'
candidate_samps <- fread(candidate_samp_file)

grp3_file <- paste(test_samp_dir, 'Grp4_subgroups.csv', sep = '')
grp3_samps <- fread(grp3_file)

### SET OUTPUTS ###


### SET VARIABLES ###
# colors used for ancestry groups in trees
gp_col_vec <- c()
gp_col_vec['ATLANTIC'] <- 'yellow2'
gp_col_vec['GULF'] <- 'red2'
gp_col_vec['MIDWEST'] <- 'blue2'
gp_col_vec['ATLANTIC+GULF'] <-'goldenrod3'
gp_col_vec['ATLANTIC+MIDWEST'] <- 'green2'
gp_col_vec['GULF+ATLANTIC'] <- 'orangered1'
gp_col_vec['GULF+MIDWEST'] <- 'orchid2'
gp_col_vec['MIDWEST+ATLANTIC'] <- 'forestgreen'
gp_col_vec['MIDWEST+GULF'] <- 'purple3'
gp_col_vec['NOT_AVAILABLE'] <- 'black'
gp_col_vec['MW_8X_WEST'] <- 'cyan2'
gp_col_vec['MW_8X_EAST'] <- 'seagreen3'
gp_col_vec['MW_4X'] <- 'dodgerblue4'
gp_col_vec['MW_8X_WEST+MW_8X_EAST'] <- 'cyan4'
gp_col_vec['MW_8X_WEST+MW_4X'] <- 'deepskyblue3'
gp_col_vec['MW_8X_EAST+MW_8X_WEST'] <- 'darkseagreen3'
gp_col_vec['MW_8X_EAST+MW_4X'] <- 'seagreen4'
gp_col_vec['MW_4X+MW_8X_WEST'] <- 'dodgerblue1'
gp_col_vec['MW_4X+MW_8X_EAST'] <- 'seagreen4'

### FUNCTIONS ###
make_adm_ancestry_vector <- function(admix_res, mw_adm_res = NA, pure_cut = 0.9,
  use_mw_adm = F){
  # Generate vector assigning samples to ancestry groups based on ADMIXTURE
  #    results
  #  Hardcoded using $V2 = Atlantic, $V3 = Gulf, $V4 = Midwest; so need to 
  #    format admix_res to have same order of columns
  #  Note: Good idea for admix_tab to have same sample order as distance matrix
  #          that will be associated with 
  # INPUTS
  # admix_res = admixture table
  # pure_cut = the ADMIXTURE coefficient used for calling a sample "pure"
  # OUTPUT
  # Vector with the index of each sample assigned an ancestry group
  ##############3
  am_atlantic_inds <- which(admix_res$V2 > 0.9)
  am_gulf_inds <- which(admix_res$V3 > 0.9)
  am_mw_inds <- which(admix_res$V4 > 0.9)
  #
#  order(admix_res[1, c(2:4)], decreasing = T)
  #
  ancestry_order <- apply(admix_res[, c(2:4)], 1, order, decreasing = T)
  ancestry_order[c(1:3), which(is.na(admix_res$V2))] <- NA
  #
  am_a_g_inds <- setdiff(
    which(ancestry_order[1,] == 1 & ancestry_order[2,] == 2), am_atlantic_inds)
  am_a_m_inds <- setdiff(
    which(ancestry_order[1,] == 1 & ancestry_order[2,] == 3), am_atlantic_inds)
  am_g_a_inds <- setdiff(
    which(ancestry_order[1,] == 2 & ancestry_order[2,] == 1), am_gulf_inds)
  am_g_m_inds <- setdiff(
    which(ancestry_order[1,] == 2 & ancestry_order[2,] == 3), am_gulf_inds)
  am_m_a_inds <- setdiff(
    which(ancestry_order[1,] == 3 & ancestry_order[2,] == 1), am_mw_inds)
  am_m_g_inds <- setdiff(
    which(ancestry_order[1,] == 3 & ancestry_order[2,] == 2), am_mw_inds)
  #
  all_gp_vec <- rep(NA, times = nrow(admix_res))
  all_gp_vec[am_atlantic_inds] <- 'ATLANTIC'
  all_gp_vec[am_gulf_inds] <- 'GULF'
  all_gp_vec[am_mw_inds] <- 'MIDWEST'
  all_gp_vec[am_a_g_inds] <- 'ATLANTIC+GULF'
  all_gp_vec[am_a_m_inds] <- 'ATLANTIC+MIDWEST'
  all_gp_vec[am_g_a_inds] <- 'GULF+ATLANTIC'
  all_gp_vec[am_g_m_inds] <- 'GULF+MIDWEST'
  all_gp_vec[am_m_a_inds] <- 'MIDWEST+ATLANTIC'
  all_gp_vec[am_m_g_inds] <- 'MIDWEST+GULF'
  all_gp_vec[which(is.na(admix_res$V2))] <- 'NOT_AVAILABLE'
  #
  if(use_mw_adm == T){
    mw_res <- admix_res
    mw_res$V2 <- mw_res$V3 <- mw_res$V4 <- NA
    mw_adm_ord <- c()
    for(i in mw_adm_res$V1){
      tmp_ind <- which(admix_res$V1 == i)
      mw_adm_ord <- c(mw_adm_ord, tmp_ind)
    }
    mw_res$V2[mw_adm_ord] <- mw_adm_res$V2
    mw_res$V3[mw_adm_ord] <- mw_adm_res$V3
    mw_res$V4[mw_adm_ord] <- mw_adm_res$V4
    #
    mw_8W <- which(mw_res$V2 > 0.9)
    mw_8E <- which(mw_res$V3 > 0.9)
    mw_4 <- which(mw_res$V4 > 0.9)
    #
    mw_anc_order <- apply(mw_res[, c(2:4)], 1, order, decreasing = T)
    mw_anc_order[c(1:3), which(is.na(mw_res$V2))] <- NA
    #
    mw_8W_8E_inds <- setdiff(
      which(mw_anc_order[1,] == 1 & mw_anc_order[2,] == 2), mw_8W)
    mw_8W_4_inds <- setdiff(
      which(mw_anc_order[1,] == 1 & mw_anc_order[2,] == 3), mw_8W)
    mw_8E_8W_inds <- setdiff(
      which(mw_anc_order[1,] == 2 & mw_anc_order[2,] == 1), mw_8E)
    mw_8E_4_inds <- setdiff(
      which(mw_anc_order[1,] == 2 & mw_anc_order[2,] == 3), mw_8E)
    mw_4_8W <- setdiff(
      which(mw_anc_order[1,] == 3 & mw_anc_order[2,] == 1), mw_4)
    mw_4_8E <- setdiff(
      which(mw_anc_order[1,] == 3 & mw_anc_order[2,] == 2), mw_4)
    #
    all_gp_vec[mw_8W] <- 'MW_8X_WEST'
    all_gp_vec[mw_8E] <- 'MW_8X_EAST'
    all_gp_vec[mw_4] <- 'MW_4X'
    all_gp_vec[mw_8W_8E_inds] <- 'MW_8X_WEST+MW_8X_EAST'
    all_gp_vec[mw_8W_4_inds] <- 'MW_8X_WEST+MW_4X'
    all_gp_vec[mw_8E_8W_inds] <- 'MW_8X_EAST+MW_8X_WEST'
    all_gp_vec[mw_8E_4_inds] <- 'MW_8X_EAST+MW_4X'
    all_gp_vec[mw_4_8W] <- 'MW_4X+MW_8X_WEST'
    all_gp_vec[mw_4_8E] <- 'MW_4X+MW_8X_EAST'
  }
  return(all_gp_vec)
}

###################

# Select filtered 8X groups
grp_1_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp1', ID],
  natv2_samps)

grp_2_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp2', ID],
  natv2_samps)

grp_4_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp4', ID],
  natv2_samps)

grp_5_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp5', ID],
  natv2_samps)

grp_3_1_filt <- intersect(grp3_samps[Sub_Group == 'subgrp1', ID],
  natv2_samps)

grp_3_2_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp2', ID],
  natv2_samps)

grp_3_2_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp2', ID],
  natv2_samps)

grp_3_3_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp3', ID],
  natv2_samps)

grp_3_4_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp4', ID],
  natv2_samps)

# Full natv2 tree
adm_nat <- admix_res[admix_res$V1 %in% natv2_samps]

dist_nat_ord <- c()
for(j in adm_nat$V1){
  tmp_ind <- which(rownames(dist_mat) == j)
  dist_nat_ord <- c(dist_nat_ord, tmp_ind)
}

nat_dist_mat <- dist_mat[dist_nat_ord, dist_nat_ord]

nat_gp_vec <- make_adm_ancestry_vector(admix_res = adm_nat)

nat_gp_colors <- nat_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(nat_gp_colors == gpname)
  nat_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

nat_dist <- as.dist(nat_dist_mat)
nat_nj <- NJ(nat_dist)

nat_nj_tree_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/natv2_nj_tree.newick'
write.tree(nat_nj, file = nat_nj_tree_file)

nat_nj_tree_pdf_out <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/natv2_nj_tree.pdf'
pdf(file = nat_nj_tree_pdf_out, height = 50, width = 10)
plot(nat_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = nat_gp_colors)
dev.off()



# Atlantic tree

atl_samps <- intersect(natv2_samps, admix_res$V1[admix_res$V2 > 0.5])

adm_atl_ord <- c()
for(i in atl_samps){
  tmp_ind <- which(admix_res$V1 == i)
  adm_atl_ord <- c(adm_atl_ord, tmp_ind)
}

adm_atl_res <- admix_res[adm_atl_ord,]

dist_atl_ord <- c()
for(j in atl_samps){
  tmp_ind <- which(rownames(dist_mat) == j)
  dist_atl_ord <- c(dist_atl_ord, tmp_ind)
}

atl_dist_mat <- dist_mat[dist_atl_ord, dist_atl_ord]

atl_gp_vec <- make_adm_ancestry_vector(admix_res = adm_atl_res)

atl_gp_colors <- atl_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(atl_gp_colors == gpname)
  atl_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

atl_dist <- as.dist(atl_dist_mat)
atl_nj <- NJ(atl_dist)

atl_nj_tree_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/atlantic_nj_tree.newick'
write.tree(atl_nj, file = atl_nj_tree_file)

atl_nj_tree_pdf_out <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/atlantic_nj_tree.pdf'
pdf(file = atl_nj_tree_pdf_out, height = 30, width = 10)
plot(atl_nj, type = 'phylogram', use.edge.length = F, no.margin = T, 
  cex = 0.3, tip.color = atl_gp_colors)
dev.off()

# Midwest tree

mw_samps <- intersect(natv2_samps, admix_res$V1[admix_res$V4 > 0.5])

adm_mw_ord <- c()
for(i in mw_samps){
  tmp_ind <- which(admix_res$V1 == i)
  adm_mw_ord <- c(adm_mw_ord, tmp_ind)
}

adm_mw_res <- admix_res[adm_mw_ord,]

dist_mw_ord <- c()
for(j in mw_samps){
  tmp_ind <- which(rownames(dist_mat) == j)
  dist_mw_ord <- c(dist_mw_ord, tmp_ind)
}

mw_dist_mat <- dist_mat[dist_mw_ord, dist_mw_ord]

ni_res_1 <- ni_res_0[ni_res_0$V1 %in% adm_mw_res$V1, ]

mw_gp_vec <- make_adm_ancestry_vector(admix_res = adm_mw_res, 
  mw_adm_res = ni_res_1, use_mw_adm = T)

mw_gp_colors <- mw_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(mw_gp_colors == gpname)
  mw_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

mw_dist <- as.dist(mw_dist_mat)
mw_nj <- NJ(mw_dist)

mw_nj_tree_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/midwest_nj_tree.newick'
write.tree(mw_nj, file = mw_nj_tree_file)

mw_nj_tree_pdf_out <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/midwest_nj_tree.pdf'
pdf(file = mw_nj_tree_pdf_out, height = 30, width = 10)
plot(mw_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = mw_gp_colors)
dev.off()

# some notes: almost all the 8X samples form a single clade, so can't rule
#  out single origin of major 8X clades
# grp3_1, grp3_2, and grep3_3 look to be most closely related to grp4 (8X West)

# Gulf tree 

g_samps <- intersect(natv2_samps, admix_res$V1[admix_res$V3 > 0.5])

adm_g_ord <- c()
for(i in g_samps){
  tmp_ind <- which(admix_res$V1 == i)
  adm_g_ord <- c(adm_g_ord, tmp_ind)
}

adm_g_res <- admix_res[adm_g_ord,]

dist_g_ord <- c()
for(j in g_samps){
  tmp_ind <- which(rownames(dist_mat) == j)
  dist_g_ord <- c(dist_g_ord, tmp_ind)
}

g_dist_mat <- dist_mat[dist_g_ord, dist_g_ord]

g_gp_vec <- make_adm_ancestry_vector(admix_res = adm_g_res, use_mw_adm = F)

g_gp_colors <- g_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(g_gp_colors == gpname)
  g_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

g_dist <- as.dist(g_dist_mat)
g_nj <- NJ(g_dist)

g_nj_tree_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/gulf_nj_tree.newick'
write.tree(g_nj, file = g_nj_tree_file)

g_nj_tree_pdf_out <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/gulf_nj_tree.pdf'
pdf(file = g_nj_tree_pdf_out, height = 30, width = 10)
plot(g_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = g_gp_colors)
dev.off()

# Gulf tree including Grp_3 samples
g_2_samps <- setdiff(
  intersect(c(admix_res$V1[admix_res$V3 > 0.5], grp3_samps$ID), natv2_samps), 
  'J530.B')

adm_g2_ord <- c()
for(i in g_2_samps){
  tmp_ind <- which(admix_res$V1 == i)
  adm_g2_ord <- c(adm_g2_ord, tmp_ind)
}

adm_g2_res <- admix_res[adm_g2_ord, ]

dist_g2_ord <- c()
for(j in g_2_samps){
  tmp_ind <- which(rownames(dist_mat) == j)
  dist_g2_ord <- c(dist_g2_ord, tmp_ind)
}

g2_dist_mat <- dist_mat[dist_g2_ord, dist_g2_ord]

g2_gp_vec <- make_adm_ancestry_vector(admix_res = adm_g2_res, use_mw_adm = F)

g2_gp_colors <- g2_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(g2_gp_colors == gpname)
  g2_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

g2_dist <- as.dist(g2_dist_mat)
g2_nj <- NJ(g2_dist)

g2_nj_tree_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/gulf_and_Grp3_nj_tree.newick'
write.tree(g2_nj, file = g2_nj_tree_file)

g2_nj_tree_pdf_out <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/gulf_and_Grp3_nj_tree.pdf'
pdf(file = g2_nj_tree_pdf_out, height = 30, width = 10)
plot(g2_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = g2_gp_colors)
dev.off()

quit(save = 'no')


