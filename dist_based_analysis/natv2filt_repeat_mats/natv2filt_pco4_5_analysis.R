### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)

### INPUT DATA ###
pcoa_dip_res_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_dist_files/GW.natv2filt.tet.100k.0001.gz_diploid_DistMat.rds'
pcoa_dip_res <- readRDS(pcoa_dip_res_file)

pcoa_tet_res_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_dist_files/GW.natv2filt.tet.100k.0001.gz_tetraploid_DistMat.rds'
pcoa_tet_res <- readRDS(pcoa_tet_res_file)

pcoa_poly_res_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_dist_files/GW.natv2filt.tet.100k.0001.gz_polyploid_DistMat.rds'
pcoa_poly_res <- readRDS(pcoa_poly_res_file)

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

natv2_samp_file <- '/Users/grabowsk/Analysis/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_samp_file, header = F)$V1

allsamp_k3_file <- '/Users/grabowsk/data/sg_8X_analysis/ADMIX_res/GW_200k_allsamps.3.results.txt'
allsamp_k3_res <- fread(allsamp_k3_file)
# $V2 = MW, $V3 = Atlantic, $V4 = Gulf

natv2_k3_file <- '/Users/grabowsk/data/sg_8X_analysis/ADMIX_res/GW_50k_natv2.3.results.txt'
natv2_k3_res <- fread(natv2_k3_file)

dip_mat_file <- '/Users/grabowsk/data/sg_8X_analysis/genotype_mats/GW.natv2.ancestral.tet.100k.0001.gz_diploid_geno_dosage_mat.rds'
dip_mat_1_raw <- readRDS(dip_mat_file)
dip_mat_1 <- dip_mat_1_raw[, -which(colnames(dip_mat_1_raw) %in% 
                                      c('ANCESTRAL_TET_GENO', 'REF_TET_GENO'))]


tet_mat_file <- '/Users/grabowsk/data/sg_8X_analysis/genotype_mats/GW.natv2.ancestral.tet.100k.0001.gz_tetraploid_geno_dosage_mat.rds'
tet_mat_1_raw <- readRDS(tet_mat_file)
tet_mat_1 <- tet_mat_1_raw[, -which(colnames(tet_mat_1_raw) %in% 
                                      c('ANCESTRAL_TET_GENO', 'REF_TET_GENO'))]
### SET VARIABLES ###
pure_cut <- 0.9
adm_cut <- 0.05

tet_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '4X')]
oct_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '8X' |
                                        samp_meta$NQUIRE_PLOIDY == '6X')]
oct_names_exact <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '8X')]
hex_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '6X')]

### SET OUTPUT ###
out_pdf <- '/Users/grabowsk/Documents/Switchgrass_8X/pcoa_reproducibility/natv2filt_PCo4vs5.pdf'

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
######

plot_tab <- data.table(samp_name = rownames(
  as.matrix(pcoa_dip_res[['euclidean_dist']])), 
                       genepool = as.character(NA), ploidy = as.character(NA))
# assign genepool and ploidy

atl_col <- 1
gulf_col <- 3
mw_col <- 2

admix_k3_res <- natv2_k3_res
ancestry_order <- apply(admix_k3_res[, c(2:4)], 1, order, decreasing = T)

mw_names <- admix_k3_res[
  which(admix_k3_res[,c(1+mw_col), with = F] > pure_cut), V1]
plot_tab[plot_tab$samp_name %in% mw_names, genepool := 'Midwest']

atl_names <- admix_k3_res[
  which(admix_k3_res[, c(1+atl_col), with = F] > pure_cut), V1]
plot_tab[plot_tab$samp_name %in% atl_names, genepool := 'Atlantic']

gulf_names <- admix_k3_res[
  which(admix_k3_res[, c(1+gulf_col), with = F] > pure_cut), V1]
plot_tab[plot_tab$samp_name %in% gulf_names, genepool := 'Gulf']

# Samples with ancestry from all 3 gene pools
multi_adm <- admix_k3_res$V1[which(admix_k3_res$V2 > adm_cut & 
                                     admix_k3_res$V3 > adm_cut & 
                                     admix_k3_res$V3 > adm_cut)]
plot_tab[plot_tab$samp_name %in% multi_adm, genepool := 'MW+ATL+GULF']

# admixed samples
mw_atl_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == mw_col & 
                          ancestry_order[2, ] == atl_col)], 
  c(mw_names, multi_adm))
mw_gulf_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == mw_col & 
                          ancestry_order[2, ] == gulf_col)], 
  c(mw_names,multi_adm))

atl_mw_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == atl_col & 
                          ancestry_order[2, ] == mw_col)], 
  c(atl_names, multi_adm))
atl_gulf_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == atl_col & 
                          ancestry_order[2, ] == gulf_col)], 
  c(atl_names, multi_adm))

gulf_mw_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == gulf_col & 
                          ancestry_order[2, ] == mw_col)], 
  c(gulf_names, multi_adm))
gulf_atl_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == gulf_col & 
                          ancestry_order[2, ] == atl_col)], 
  c(gulf_names, multi_adm))

plot_tab[plot_tab$samp_name %in% mw_atl_names, genepool := 'MW+ATL']
plot_tab[plot_tab$samp_name %in% mw_gulf_names, genepool := 'MW+GULF']
plot_tab[plot_tab$samp_name %in% atl_mw_names, genepool := 'ATL+MW']
plot_tab[plot_tab$samp_name %in% atl_gulf_names, genepool := 'ATL+GULF']
plot_tab[plot_tab$samp_name %in% gulf_mw_names, genepool := 'GULF+MW']
plot_tab[plot_tab$samp_name %in% gulf_atl_names, genepool := 'GULF+ATL']

plot_tab[plot_tab$samp_name %in% tet_names, ploidy := '4X']
plot_tab[plot_tab$samp_name %in% hex_names, ploidy := '6X']
plot_tab[plot_tab$samp_name %in% oct_names_exact, ploidy := '8X']
plot_tab[which(is.na(plot_tab$ploidy)), ploidy := '4X']

####

dip_pcoa <- cmdscale(as.matrix(pcoa_dip_res[['euclidean_dist']]), k = 40)

plot_tab[, dip_PCo4 := dip_pcoa[,4]]
plot_tab[, dip_PCo5 := dip_pcoa[,5]]

gg_dip4_5 <- ggplot(plot_tab, aes(x = dip_PCo4, y = dip_PCo5, color = genepool)) +
  geom_point() +
  labs(title = 'natv2_filtered diploid genotypes\nPCo4 vs PCo5')
  
tet_pcoa <- cmdscale(as.matrix(pcoa_tet_res[['euclidean_dist']]), k = 40)

plot_tab[, tet_PCo4 := tet_pcoa[,4]]
plot_tab[, tet_PCo5 := tet_pcoa[,5]]

gg_tet4_5 <- ggplot(plot_tab, aes(x = tet_PCo4, y = tet_PCo5, color = genepool)) +
  geom_point() +
  labs(title = 'natv2_filtered tetraploid genotypes\nPCo4 vs PCo5')

poly_pcoa <- cmdscale(as.matrix(pcoa_poly_res[['euclidean_dist']]), k = 40)

plot_tab[, poly_PCo4 := poly_pcoa[,4]]
plot_tab[, poly_PCo5 := poly_pcoa[,5]]

gg_poly4_5 <- ggplot(plot_tab, aes(x = poly_PCo4, y = poly_PCo5, color = genepool)) +
  geom_point() +
  labs(title = 'natv2_filtered polyploid genotypes\nPCo4 vs PCo5')

pdf(file = out_pdf, height = 5*3, width = 6)
gg_dip4_5 / gg_tet4_5 / gg_poly4_5
dev.off()
