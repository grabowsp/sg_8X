# Random analysis for JGI Plant Group Meeting on July 1 2021

# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library('bit', lib.loc = '/home/grabowsky/tools/r_packages')
library('bit64', lib.loc = '/home/grabowsky/tools/r_packages')
library(ggplot2)
library(ggbeeswarm, lib.loc = '/home/grabowsky/tools/r_packages')

### INPUT DATA ###
meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

natv2_samp_file <- '/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_samp_file, header = F)

prelim_ploidy_file <- '/home/f2p1/work/grabowsk/data/switchgrass/pseudohap/sg_ploidy_results_v3.0.txt'
prelim_ploidy <- fread(prelim_ploidy_file)

res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_file)

### SET OUTPUTS ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/'

##########

### Info about accession and sample numbers
length(unique(samp_meta$ACC))
# 535 Accessions

length(unique(samp_meta$UNI_ACC))
# 486 Unique Accessions

sum(samp_meta$NAT_PAPER_LIB == 'Y', na.rm = T)
# 733 samples in the Nature paper

nrow(natv2_samps)
# 706 samples in ploidy analysis

meta_natv2_inds <- which(samp_meta$VCF_NAME %in% natv2_samps$V1)

# Samples not in the Nature Paper
new_samp_inds <- setdiff(meta_natv2_inds, which(samp_meta$NAT_PAPER_LIB == 'Y'))

length(new_samp_inds)
# 178 new samples

length(setdiff(which(samp_meta$NAT_PAPER_LIB == 'Y'), meta_natv2_inds))
# 205 samples from Nature Paper excluded - cultivars and geographic questions

samp_meta[new_samp_inds, .N, by = 'NQUIRE_PLOIDY']
#   NQUIRE_PLOIDY   N
#1:            6X   2
#2:            8X 152
#3:            4X  24

### nQuire plot(s)
# nQuire plot to explain program
nSNP_cut <- 500000
prelim_ploidy[, COV_GROUP := as.character(NA)]
prelim_ploidy[nquire_nSNPS_20 < nSNP_cut, COV_GROUP := 'LOW']
prelim_ploidy[nquire_nSNPS_20 >= nSNP_cut, COV_GROUP := 'HIGH']

gg_nquire_1 <- ggplot(prelim_ploidy, 
  aes(x = nquire_d_dip_portion_20, y = nquire_d_tet_portion_20, 
    color = COV_GROUP)) +
  geom_point() +
  xlab('2N Score') +
  ylab('4N Score')

nquire_example_file <- paste(out_dir, 'nquire_example_dotplot.pdf', sep = '')

pdf(nquire_example_file, width = 5, height = 5)
gg_nquire_1
dev.off()

### PCA plots
res_tab[, genepool := as.character(NA)]
res_tab[grep('MW_', res_tab$subgrp_v2), genepool := 'Midwest']
res_tab[grep('GULF_', res_tab$subgrp_v2), genepool := 'Gulf']
res_tab[grep('ATL_', res_tab$subgrp_v2), genepool := 'Atlantic']
res_tab[grep('MF_', res_tab$subgrp_v2), genepool := 'Admixed']

ploidy_shape_vec <- c()
ploidy_shape_vec['4X'] <- 15
ploidy_shape_vec['6X'] <- 17
ploidy_shape_vec['8X'] <- 16

ploidy_shape_set <- scale_shape_manual(name = 'Ploidy', 
  values = ploidy_shape_vec)

gp_color_vec <- c()
gp_color_vec['Midwest'] <- 'royalblue1' #'blue2'
gp_color_vec['Atlantic'] <- 'yellow3'
gp_color_vec['Gulf'] <- 'red2'
gp_color_vec['Admixed'] <- 'grey50'

gp_palette <- scale_colour_manual(name = 'Genepool', values = gp_color_vec)

gg_pca_4X <- ggplot(res_tab[ploidy == '4X'], 
    aes(x = full_pc01_raw, y = full_pc02_raw,
    color = genepool, shape = ploidy)) +
  geom_point() +
  ploidy_shape_set +
  gp_palette +
  xlab('PC1 (17.9%)') +
  ylab('PC2 (6.3%)')

pca_4X_file <- paste(out_dir, 'PCA_4Xonly_dotplot.pdf', sep = '')

pdf(pca_4X_file, width = 6, height = 5)
gg_pca_4X
dev.off()

oct_pc1_vals <- res_tab[, full_pc01_raw]
oct_pc1_vals[which(res_tab$ploidy != '8X')] <- NA
oct_pc2_vals <- res_tab[, full_pc02_raw]
oct_pc2_vals[which(res_tab$ploidy != '8X')] <- NA

hex_pc1_vals <- res_tab[, full_pc01_raw]
hex_pc1_vals[which(res_tab$ploidy != '6X')] <- NA
hex_pc2_vals <- res_tab[, full_pc02_raw]
hex_pc2_vals[which(res_tab$ploidy != '6X')] <- NA

g_hex_pc1_vals <- res_tab[, full_pc01_raw]
g_hex_pc1_vals[union(which(res_tab$ploidy != '6X'), 
  which(res_tab$genepool != 'Gulf'))] <- NA
g_hex_pc2_vals <- res_tab[, full_pc02_raw]
g_hex_pc2_vals[union(which(res_tab$ploidy != '6X'), 
  which(res_tab$genepool != 'Gulf'))] <- NA

g_oct_pc1_vals <- res_tab[, full_pc01_raw]
g_oct_pc1_vals[union(which(res_tab$ploidy != '8X'),    
  which(res_tab$genepool != 'Gulf'))] <- NA
g_oct_pc2_vals <- res_tab[, full_pc02_raw]
g_oct_pc2_vals[union(which(res_tab$ploidy != '8X'),     
  which(res_tab$genepool != 'Gulf'))] <- NA

a_hex_pc1_vals <- res_tab[, full_pc01_raw]
a_hex_pc1_vals[union(which(res_tab$ploidy != '6X'),
  which(res_tab$genepool != 'Atlantic'))] <- NA
a_hex_pc2_vals <- res_tab[, full_pc02_raw]
a_hex_pc2_vals[union(which(res_tab$ploidy != '6X'),
  which(res_tab$genepool != 'Atlantic'))] <- NA

a_oct_pc1_vals <- res_tab[, full_pc01_raw]
a_oct_pc1_vals[union(which(res_tab$ploidy != '8X'),
  which(res_tab$genepool != 'Atlantic'))] <- NA
a_oct_pc2_vals <- res_tab[, full_pc02_raw]
a_oct_pc2_vals[union(which(res_tab$ploidy != '8X'),
  which(res_tab$genepool != 'Atlantic'))] <- NA

mw_hex_pc1_vals <- res_tab[, full_pc01_raw]
mw_hex_pc1_vals[union(which(res_tab$ploidy != '6X'),
  which(res_tab$genepool != 'Midwest'))] <- NA
mw_hex_pc2_vals <- res_tab[, full_pc02_raw]
mw_hex_pc2_vals[union(which(res_tab$ploidy != '6X'),
  which(res_tab$genepool != 'Midwest'))] <- NA

mw_oct_pc1_vals <- res_tab[, full_pc01_raw]
mw_oct_pc1_vals[union(which(res_tab$ploidy != '8X'),
  which(res_tab$genepool != 'Midwest'))] <- NA
mw_oct_pc2_vals <- res_tab[, full_pc02_raw]
mw_oct_pc2_vals[union(which(res_tab$ploidy != '8X'),
  which(res_tab$genepool != 'Midwest'))] <- NA

mf_hex_pc1_vals <- res_tab[, full_pc01_raw]
mf_hex_pc1_vals[union(which(res_tab$ploidy != '6X'),
  which(res_tab$genepool != 'Admixed'))] <- NA
mf_hex_pc2_vals <- res_tab[, full_pc02_raw]
mf_hex_pc2_vals[union(which(res_tab$ploidy != '6X'),
  which(res_tab$genepool != 'Admixed'))] <- NA

mf_oct_pc1_vals <- res_tab[, full_pc01_raw]
mf_oct_pc1_vals[union(which(res_tab$ploidy != '8X'),
  which(res_tab$genepool != 'Admixed'))] <- NA
mf_oct_pc2_vals <- res_tab[, full_pc02_raw]
mf_oct_pc2_vals[union(which(res_tab$ploidy != '8X'),
  which(res_tab$genepool != 'Admixed'))] <- NA

gg_pca_all <- ggplot(res_tab,
    aes(x = full_pc01_raw, y = full_pc02_raw,
    color = genepool, shape = ploidy)) +
  geom_point() + 
  ploidy_shape_set +
  gp_palette +
  xlab('PC1 (17.9%)') +
  ylab('PC2 (6.3%)') + 
#  geom_point(aes(x = oct_pc1_vals, y = oct_pc2_vals), 
#    color = 'black', shape = 1) +
#  geom_point(aes(x = hex_pc1_vals, y = hex_pc2_vals), 
#    color = 'black', shape = 2) +
  geom_point(aes(x = g_hex_pc1_vals, y = g_hex_pc2_vals),
    fill = gp_color_vec['Gulf'], shape = 24, color = 'black') +
  geom_point(aes(x = g_oct_pc1_vals, y = g_oct_pc2_vals),
    fill = gp_color_vec['Gulf'], shape = 21, color = 'black') +
  geom_point(aes(x = a_hex_pc1_vals, y = a_hex_pc2_vals),
    fill = gp_color_vec['Atlantic'], shape = 24, color = 'black') +
  geom_point(aes(x = a_oct_pc1_vals, y = a_oct_pc2_vals),
    fill = gp_color_vec['Atlantic'], shape = 21, color = 'black') +
  geom_point(aes(x = mw_hex_pc1_vals, y = mw_hex_pc2_vals),
    fill = gp_color_vec['Midwest'], shape = 24, color = 'black') +
  geom_point(aes(x = mw_oct_pc1_vals, y = mw_oct_pc2_vals),
    fill = gp_color_vec['Midwest'], shape = 21, color = 'black') +
  geom_point(aes(x = mf_hex_pc1_vals, y = mf_hex_pc2_vals),
    fill = gp_color_vec['Admixed'], shape = 24, color = 'black') +
  geom_point(aes(x = mf_oct_pc1_vals, y = mf_oct_pc2_vals),
    fill = gp_color_vec['Admixed'], shape = 21, color = 'black')


pca_all_file <- paste(out_dir, 'PCA_allsamps_dotplot.pdf', sep = '')

pdf(pca_all_file, width = 6, height = 5)
gg_pca_all
dev.off()

### Summary of 8X in each genepool
mw8X_inds <- intersect(grep('MW_', res_tab$subgrp_v2), 
  which(res_tab$ploidy == '8X'))

length(mw8X_inds)
# 105

summary(res_tab[mw8X_inds, full_admix_k3_MW])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7160  0.9730  0.9986  0.9608  1.0000  1.0000

sum(res_tab[mw8X_inds, full_admix_k3_MW] > 0.99)
# 71 MW_8X have 99%+ Midwest ancestry

g8X_inds <- intersect(grep('GULF_', res_tab$subgrp_v2), 
  which(res_tab$ploidy == '8X'))

length(g8X_inds)
# [1] 17

summary(res_tab[g8X_inds, full_admix_k3_GULF])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7567  0.8027  0.9039  0.8811  0.9695  0.9909

sum(res_tab[g8X_inds, full_admix_k3_GULF] > 0.99)
# 1
sum(res_tab[g8X_inds, full_admix_k3_GULF] > 0.97)
# 4

a8X_inds <- intersect(grep('ATL_', res_tab$subgrp_v2), 
  which(res_tab$ploidy == '8X'))

length(a8X_inds)
# [1] 11

summary(res_tab[a8X_inds, full_admix_k3_ATL])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.6977  0.7082  0.7427  0.7833  0.7985  1.0000

sum(res_tab[a8X_inds, full_admix_k3_ATL] > 0.99)
# 1 sample with 99%+ (and 95%+) Atlantic ancestry

max_top_ancestry <- apply(res_tab[, c('full_admix_k3_ATL', 
  'full_admix_k3_GULF', 'full_admix_k3_MW'), with = F], 1, max)

sum(max_top_ancestry[which(res_tab$ploidy == '4X')] > 0.97)
# 468
length(which(res_tab$ploidy == '4X'))
# 544
# 468/544 = 0.86

sum(max_top_ancestry[which(res_tab$ploidy == '4X')] <= 0.8)
# 20
# 20/544 = 0.037

sum(max_top_ancestry[which(res_tab$ploidy == '8X')] > 0.97)
# 84
length(which(res_tab$ploidy == '8X'))
# 159
# 84/159 = 0.53

sum(max_top_ancestry[which(res_tab$ploidy == '8X')] <= 0.8)
# 46/159 = 0.289

# Boxplot of top ancestry

max_top_ancestry <- apply(res_tab[, c('full_admix_k3_ATL',
  'full_admix_k3_GULF', 'full_admix_k3_MW'), with = F], 1, max)

res_tab[, per_top_ancestry := round(max_top_ancestry, digits = 2)]
res_tab[, ploidy := as.factor(res_tab$ploidy)]

gg_ancest_violin <- ggplot(res_tab, aes(y = per_top_ancestry, x = ploidy, 
    fill = ploidy)) +
  geom_violin() +
  geom_quasirandom() +
#  scale_fill_manual(values = c('deeppink', 'chartreuse3')) +
  scale_x_discrete(limits=c('4X', '8X'))


ancestry_violin_file <- paste(out_dir, 'topancestry_ploidy_violinplot.pdf', 
  sep = '')

pdf(ancestry_violin_file, width = 6, height = 5)
gg_ancest_violin
dev.off()





