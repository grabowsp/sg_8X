# Big NJ Tree using all samples and UNI_ACC as labels

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

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

admix_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt'
admix_res_0 <- fread(admix_res_file)

ni_admix_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_NorthInland_admix/GW_50k_geobigNorthInland.3.results.txt'
ni_res_0 <- fread(ni_admix_file)

#####################

admix_ord <- c()
for(j in seq(nrow(admix_res_0))){
  tmp_ind <- which(rownames(dist_mat) == admix_res_0$V1[j])
  admix_ord <- c(admix_ord, tmp_ind)
}

admix_res <- data.table(V1 = rep(as.character(NA), times = nrow(dist_mat)),
  V2 = rep(as.numeric(NA), times = nrow(dist_mat)),
  V3 = rep(as.numeric(NA), times = nrow(dist_mat)),
  V4 = rep(as.numeric(NA), times = nrow(dist_mat))
)

admix_res[admix_ord, V1 := admix_res_0$V1]
admix_res[admix_ord, V2 := admix_res_0$V2]
admix_res[admix_ord, V3 := admix_res_0$V3]
admix_res[admix_ord, V4 := admix_res_0$V4]

#admix_res$V1[admix_res$V4 > 0.9]
# 1 (V2) = Atlantic
#2 (V3) = Gulf
# 3 (V4) = Midwest

am_atlantic_inds <- which(admix_res$V2 > 0.9)
am_gulf_inds <- which(admix_res$V3 > 0.9)
am_mw_inds <- which(admix_res$V4 > 0.9)

order(admix_res[1, c(2:4)], decreasing = T)

ancestry_order <- apply(admix_res[, c(2:4)], 1, order, decreasing = T)
ancestry_order[c(1:3), which(is.na(admix_res$V2))] <- NA

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

all_gp_vec <- rep(NA, times = nrow(dist_mat))
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

allsamps_gp_colors <- all_gp_vec
for(gpname in names(gp_col_vec)){
  tmp_inds <- which(allsamps_gp_colors == gpname)
  allsamps_gp_colors[tmp_inds] <- gp_col_vec[gpname]
}

allsamp_dist <- as.dist(dist_mat)
allsamp_nj <- NJ(allsamp_dist)

allsamp_sampname_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_sampnames.pdf'

pdf(file = allsamp_sampname_outfile, height = 50, width = 10)
plot(allsamp_nj, type = 'phylogram', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

allsamp_sampname_outfile_nolength <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_sampnames_nolength.pdf'

pdf(file = allsamp_sampname_outfile_nolength, height = 50, width = 10)
plot(allsamp_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

meta_ord <- c()
for(i in seq(nrow(dist_mat))){
  tmp_ind <- which(samp_meta$VCF_NAME == rownames(dist_mat)[i])
  meta_ord <- c(meta_ord, tmp_ind)
}

UA_mat <- dist_mat
rownames(UA_mat) <- colnames(UA_mat) <- samp_meta$UNI_ACC[meta_ord]

UA_dist <- as.dist(UA_mat)
UA_nj <- NJ(UA_dist)

allsamp_UA_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_unified_accesions.pdf'

pdf(file = allsamp_UA_outfile, height = 50, width = 10)
plot(UA_nj, type = 'phylogram', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

allsamp_UA_outfile_nolength <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_unified_accesions_nolength.pdf'

pdf(file = allsamp_UA_outfile_nolength, height = 50, width = 10)
plot(UA_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

samp_meta$STATE[meta_ord]
samp_meta[which(is.na(samp_meta$STATE)), VCF_NAME]

samp_meta[which(is.na(samp_meta$STATE))[1], ]

samp_meta[intersect(which(is.na(samp_meta$STATE)), grep('MX_', samp_meta$VCF_NAME)), STATE := 'MX']
samp_meta[VCF_NAME == 'J007_C', STATE := "Argentina"]
samp_meta[VCF_NAME == 'J043_A', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J448_B', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J285_A', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J285_B', STATE := 'Hybrid']

samp_meta[VCF_NAME == 'J608_B', STATE := 'Unknown']
samp_meta[VCF_NAME == 'J652_A', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J419.B', STATE := 'Hybrid']

samp_meta[VCF_NAME == 'J301.A', STATE := 'Unknown_NAM']
samp_meta[intersect(which(is.na(samp_meta$STATE)), 
  grep('NAM_', samp_meta$UNI_ACC)), STATE := 'Unknown_NAM']

samp_meta[VCF_NAME == 'J419.A', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J448.A', STATE := 'Hybrid']

samp_meta[VCF_NAME == 'J608.A', STATE := 'Unknown']
samp_meta[VCF_NAME == 'J651.A', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J651.B', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J651.C', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J652.B', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J652.C', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'J681.A', STATE := 'Hybrid']

samp_meta[VCF_NAME == 'NFGA34_08', STATE := 'TX']
samp_meta[VCF_NAME == 'NFGA16_12', STATE := 'NC']
samp_meta[VCF_NAME == 'NFGA34_10', STATE := 'TX']
samp_meta[VCF_NAME == 'NFGA37_05', STATE := 'Unknown']
samp_meta[VCF_NAME == 'NFGA15_11', STATE := 'MS']

samp_meta[VCF_NAME == 'NFGA09_02', STATE := 'OK']
samp_meta[VCF_NAME == 'NFGA09_05', STATE := 'OK']
samp_meta[VCF_NAME == 'NFGA32_06', STATE := 'TX']

samp_meta[VCF_NAME == 'NFGA16_05', STATE := 'NC']
samp_meta[VCF_NAME == 'NFGA02_06', STATE := 'AR']
samp_meta[VCF_NAME == 'NFGA37_03', STATE := 'Unknown']

samp_meta[VCF_NAME == 'NL94', STATE := 'OK']
samp_meta[VCF_NAME == 'SL93', STATE := 'TX']
samp_meta[VCF_NAME == 'Pv346_1', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'GAasmb6', STATE := 'FL']
samp_meta[VCF_NAME == 'Pv458_1', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'Pv304_1', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'Pv317_1', STATE := 'Hybrid']

samp_meta[VCF_NAME == 'NFGA16_02', STATE := 'NC']
samp_meta[VCF_NAME == '9001-3_BN389-69S', STATE := 'NC']
samp_meta[VCF_NAME == '9020-5', STATE := 'NC']

samp_meta[VCF_NAME == 'SW786_WO1', STATE := 'Contaminant']
samp_meta[VCF_NAME == '2302', STATE := 'OK']
samp_meta[VCF_NAME == '2785b', STATE := 'Hybrid']
samp_meta[VCF_NAME == 'NL94_85_1', STATE := 'OK']
samp_meta[VCF_NAME == 'SL93_16_1', STATE := 'TX']

samp_meta[VCF_NAME == '2732', STATE := 'OK']
samp_meta[VCF_NAME == 'NL94_85_1_5', STATE := 'OK']
samp_meta[VCF_NAME == 'SL93_16', STATE := 'TX']

state_mat <- dist_mat
rownames(state_mat) <- colnames(state_mat) <- samp_meta$STATE[meta_ord]

state_dist <- as.dist(state_mat)
state_nj <- NJ(state_dist)

allsamp_state_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_state.pdf'

pdf(file = allsamp_state_outfile, height = 50, width = 10)
plot(state_nj, type = 'phylogram', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

allsamp_state_outfile_nolength <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_state_nolength.pdf'

pdf(file = allsamp_state_outfile_nolength, height = 50, width = 10)
plot(state_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

####

combined_names <- paste(paste(samp_meta$VCF_NAME[meta_ord], 
  samp_meta$UNI_ACC[meta_ord], sep = '+'), samp_meta$STATE[meta_ord], sep = '+')

combo_mat <- dist_mat
rownames(combo_mat) <- colnames(combo_mat) <- combined_names

combo_dist <- as.dist(combo_mat)
combo_nj <- NJ(combo_dist)

allsamp_combo_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_combolabel.pdf'

pdf(file = allsamp_combo_outfile, height = 50, width = 10)
plot(combo_nj, type = 'phylogram', use.edge.length = T, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

allsamp_combo_outfile_nolength <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/allsamp_NJ_tree_combolabel_nolength.pdf'

pdf(file = allsamp_combo_outfile_nolength, height = 50, width = 10)
plot(combo_nj, type = 'phylogram', use.edge.length = F, no.margin = T,
  cex = 0.3, tip.color = allsamps_gp_colors)
dev.off()

quit(save = 'no')

