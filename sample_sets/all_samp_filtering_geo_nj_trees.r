# NJ Tree using only Natural Samples

### LOAD MODULES ###

#module load python/3.7-anaconda-2019.10
#module swap PrgEnv-intel PrgEnv-gnu
#source activate R_tidy

# install.packages('ape', lib = '/global/homes/g/grabowsp/tools/r_libs')
# install.packages('phangorn', lib = '/global/homes/g/grabowsp/tools/r_libs')

### LOAD PACKAGES ###
library(data.table)
library(ape, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
library(phangorn, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')

### INPUT DATA ###
meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v2.0.csv'
samp_meta <- fread(meta_file)

dist_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/GW_all_samps_ld0.3_symmetric.mdist'
dist_res <- fread(dist_res_file, header = F)
dist_ids <- fread(paste(dist_res_file, 'id', sep = '.'), header = F)
dist_mat <- as.matrix(dist_res)
colnames(dist_mat) <- rownames(dist_mat) <- dist_ids$V1

geobig_name_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_big_names.txt'
geobig_names <- fread(geobig_name_file, header=F)$V1

admix_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt'
admix_res_0 <- fread(admix_res_file)

### SET OUTPUT ###

geobig_ploidy_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/geobig_NJ_tree_ploidy.pdf'

geobig_gp_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/geobig_NJ_tree_genepools.pdf'

##################

geobig_mat <- dist_mat[geobig_names, geobig_names]
geobig_dist <- as.dist(geobig_mat)
geobig_nj <- NJ(geobig_dist)

meta_ord <- c()
for(j in seq(nrow(geobig_mat))){
  tmp_ind <- which(samp_meta$VCF_NAME == rownames(geobig_mat)[j])
  meta_ord <- c(meta_ord, tmp_ind)
}

ploid_col_vec <- c()
ploid_col_vec['8X'] <- 'turquoise3'
ploid_col_vec['6X'] <- 'hotpink2'
ploid_col_vec['4X'] <- 'black'

geobig_ploidy_cols <- samp_meta$NQUIRE_PLOIDY[meta_ord]
for(SP in names(ploid_col_vec)){
  tmp_inds <- which(geobig_ploidy_cols == SP)
  geobig_ploidy_cols[tmp_inds] <- ploid_col_vec[SP]
}


pdf(file = geobig_ploidy_outfile)
plot(geobig_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = geobig_ploidy_cols)
dev.off()

# color by gene pool

admix_ord <- c()
for(j in seq(nrow(geobig_mat))){
  tmp_ind <- which(admix_res_0 == rownames(geobig_mat)[j])
  admix_ord <- c(admix_ord, tmp_ind)
}

admix_res <- admix_res_0[admix_ord, ]
#admix_res$V1[admix_res$V4 > 0.9]
# 1 (V2) = Atlantic
#2 (V3) = Gulf
# 3 (V4) = Midwest

am_atlantic_inds <- which(admix_res$V2 > 0.9)
am_gulf_inds <- which(admix_res$V3 > 0.9) 
am_mw_inds <- which(admix_res$V4 > 0.9)

order(admix_res[1, c(2:4)], decreasing = T)

ancestry_order <- apply(admix_res[, c(2:4)], 1, order, decreasing = T)
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

geobig_gp_vec <- rep(NA, times = nrow(geobig_mat))
geobig_gp_vec[am_atlantic_inds] <- 'ATLANTIC'
geobig_gp_vec[am_gulf_inds] <- 'GULF'
geobig_gp_vec[am_mw_inds] <- 'MIDWEST'
geobig_gp_vec[am_a_g_inds] <- 'ATLANTIC+GULF'
geobig_gp_vec[am_a_m_inds] <- 'ATLANTIC+MIDWEST'
geobig_gp_vec[am_g_a_inds] <- 'GULF+ATLANTIC'
geobig_gp_vec[am_g_m_inds] <- 'GULF+MIDWEST'
geobig_gp_vec[am_m_a_inds] <- 'MIDWEST+ATLANTIC'
geobig_gp_vec[am_m_g_inds] <- 'MIDWEST+GULF'

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

geobig_gp_colors <- geobig_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(geobig_gp_colors == gpname)
  geobig_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

pdf(file = geobig_gp_outfile)
plot(geobig_nj, type = 'unrooted', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = geobig_gp_colors)
dev.off()

quit(save = 'no')


