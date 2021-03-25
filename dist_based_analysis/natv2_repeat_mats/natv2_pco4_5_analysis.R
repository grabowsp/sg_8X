### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
pcoa_list_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2.100k.pcoa_mats.list.rds'
pcoa_list <- readRDS(pcoa_list_file)

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

plot_tab <- data.table(samp_name = rownames(pcoa_list[[1]][[1]]), 
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
dip_pcoa_list <- pcoa_list[['dip_pcoa_list']]
tet_pcoa_list <- pcoa_list[['tet_pcoa_list']]
poly_pcoa_list <- pcoa_list[['poly_pcoa_list']]

# PCo4 mean values
dip_pco_4_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
                                         pco_num = 4, control_name = 'J645.C')
plot_tab[ , dip_pco_4_mean := apply(dip_pco_4_stand_mat, 2, mean)]
#
tet_pco_4_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
                                        pco_num = 4, control_name = 'J645.C')
plot_tab[ , tet_pco_4_mean := apply(tet_pco_4_stand_mat, 2, mean)]

# PCo5 mean values
dip_pco_5_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
                                        pco_num = 5, control_name = 'J635.C')
plot_tab[ , dip_pco_5_mean := apply(dip_pco_5_stand_mat, 2, mean)]
#
tet_pco_5_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
                                        pco_num = 5, control_name = 'J635.C')
plot_tab[ , tet_pco_5_mean := apply(tet_pco_5_stand_mat, 2, mean)]

#sum(plot_tab$samp_name == colnames(dip_mat_1[, c(3:ncol(dip_mat_1))]))
# 706

###### PCo4 Correlations ######
#### Calculate correlations ####
# Diploid genotypes vs diploid PCoA results
dipG_v_dipP_cor_m1_pco4 <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))], 1, cor, 
                         y = plot_tab$dip_pco_4_mean)
dipG_v_dipP_m1_pco4_hi_cor_snps <- order(abs(dipG_v_dipP_cor_m1_pco4), 
                                         decreasing = T)[seq(200)]

# tetraploid genotypes vs tetraploid PCoA results
tetG_v_tetP_cor_m1_pco4 <- apply(tet_mat_1[, c(3:ncol(tet_mat_1))], 1, cor, 
                                 y = plot_tab$tet_pco_4_mean)
tetG_v_tetP_m1_pco4_hi_cor_snps <- order(abs(tetG_v_tetP_cor_m1_pco4), 
                                        decreasing = T)[seq(200)]

# diploid genotypes vs tetraploid PCoA results
dipG_v_tetP_cor_m1_pco4 <- apply(dip_mat_1[, c(3:ncol(tet_mat_1))], 1, cor, 
                                 y = plot_tab$tet_pco_4_mean)
dipG_v_tetP_m1_pco4_hi_cor_snps <- order(abs(dipG_v_tetP_cor_m1_pco4), 
                                         decreasing = T)[seq(200)]

# tetraploid genotypes vs diploid PCoA results
tetG_v_dipP_cor_m1_pco4 <- apply(tet_mat_1[, c(3:ncol(tet_mat_1))], 1, cor, 
                                 y = plot_tab$dip_pco_4_mean)
tetG_v_dipP_m1_pco4_hi_cor_snps <- order(abs(tetG_v_dipP_cor_m1_pco4), 
                                         decreasing = T)[seq(200)]

#### Summarize correlations ####
# summary of top disomic correlations
summary(abs(dipG_v_dipP_cor_m1_pco4[dipG_v_dipP_m1_pco4_hi_cor_snps]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2493  0.2728  0.2978  0.3167  0.3401  0.5651

# summary of top tetrasomic correlations
summary(abs(tetG_v_tetP_cor_m1_pco4[tetG_v_tetP_m1_pco4_hi_cor_snps]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.2426  0.2600  0.2851  0.3021  0.3306  0.5323

# summary of top dip-genotype v tet PCo4 correlations
summary(abs(dipG_v_tetP_cor_m1_pco4[dipG_v_tetP_m1_pco4_hi_cor_snps]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2380  0.2593  0.2832  0.3018  0.3333  0.5366

# summary of top tet-genotype v dip-PCo4 correlations
summary(abs(tetG_v_dipP_cor_m1_pco4[tetG_v_dipP_m1_pco4_hi_cor_snps]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2390  0.2579  0.2795  0.3022  0.3242  0.5649

## Observations: 
# generally see high SNPs with same distribution of high correlations with
#  PCo4 when comparing within or between ploidy genotypes v PCo4

#### What are correlation for hi-cor SNPs in other comparisons ####

# D_G_vs_D_P hi-cor SNPs in T_G_vs_T_P correlations
summary(abs(tetG_v_tetP_cor_m1_pco4[dipG_v_dipP_m1_pco4_hi_cor_snps]))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0005452 0.0120110 0.0194873 0.0781855 0.0485603 0.5323036
# many tet-correlations are LOW for SNPs with high dip-correlations

# T_G_vs_T_P hi-cor SNPs in D_G_vs_D_P_correlations
summary(abs(dipG_v_dipP_cor_m1_pco4[tetG_v_tetP_m1_pco4_hi_cor_snps]))
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.03323 0.15587 0.18012 0.18696 0.22213 0.33986
# many dip-correlations remain HI for SNPs with high tet-correlations
#  - not as high, but still elevated

# D_G_vs_D_P hi-cor SNPs in D_G_vs_T_P correlations
summary(abs(dipG_v_tetP_cor_m1_pco4[dipG_v_dipP_m1_pco4_hi_cor_snps]))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000014 0.0141402 0.0229498 0.0798698 0.0489038 0.5365862
# The tetrasomic-PCo4 values seem to be sufficiently different so that
#   the SNPs with highest correlations to diploid-PCo4 aren't correlated
#  with the tetraploid-PCo4 values

# D_G_vs_D_P hi-cor SNPs in T_G_vs_D_P correlations
summary(abs(tetG_v_dipP_cor_m1_pco4[dipG_v_dipP_m1_pco4_hi_cor_snps]))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.1670  0.2511  0.2792  0.2967  0.3242  0.5649
# tetrasomic-genotypes show elevated correlations to diploid-PCo4

## Observations:
# It seems that PCo-4 values change enough that SNPs correlated to 
#   disomic-PCo4 using either genotype no longer show high correlations;
#   maybe it's because there is no longer the same variation in the Midwest?

#### Overlap in top-correlations between genotype classes ####
# how many SNPs are top-correlated in both genotype classes
length(intersect(dipG_v_dipP_m1_pco4_hi_cor_snps, 
                 tetG_v_tetP_m1_pco4_hi_cor_snps))
# 27
hi_cor_both_snps <- intersect(dipG_v_dipP_m1_pco4_hi_cor_snps, 
                              tetG_v_tetP_m1_pco4_hi_cor_snps)
# What populations are HET in hi_cor_both_snps
hi_cor_het_tabs <- list()
for(i in seq(length(hi_cor_both_snps))){
  hi_cor_het_tabs[[i]] <- table(plot_tab$genepool[
  which(dip_mat_1[hi_cor_both_snps[i], c(3:ncol(dip_mat_1))] == 1)])
}
# Always with Atlantic, often with some Gulf, MW and/or admixed individuals

# Out of samples with HET genotypes, what proportion are Atlantic
hi_cor_ATL_het_portion <- unlist(lapply(hi_cor_het_tabs, function(x)
  x['Atlantic']/sum(x)))

# What ploidy is HET in hi_cor_both_snps
hi_cor_het_pl_tabs <- list()
for(i in seq(length(hi_cor_both_snps))){
  hi_cor_het_pl_tabs[[i]] <- table(plot_tab$ploidy[
    which(dip_mat_1[hi_cor_both_snps[i], c(3:ncol(dip_mat_1))] == 1)])
}

# Out of samples with HET genotypes, what percentage is 4X
hi_cor_4X_het_portion <- unlist(lapply(hi_cor_het_pl_tabs, function(x)
  x['4X']/sum(x)))

#table(unlist(dip_mat_1[hi_cor_both_snps[1], c(3:ncol(dip_mat_1))]))
#table(unlist(tet_mat_1[hi_cor_both_snps[1], c(3:ncol(tet_mat_1))]))

# How does HET dosage change from dip to tet genotypes in hi-cor-both SNPs
hi_cor_dosage_tabs <- list()
for(i in seq(length(hi_cor_both_snps))){
  hi_cor_dosage_tabs[[i]] <- table(unlist(tet_mat_1[hi_cor_both_snps[i], 
                                             c(3:ncol(tet_mat_1))]))
}

samp_meta[samp_meta$VCF_NAME %in% 
            plot_tab$samp_name[which(plot_tab$tet_pco_4_mean > 0.7 & 
                                       plot_tab$genepool == 'Atlantic')], 
          STATE]
# Mainly Northern Atlantic

samp_meta[samp_meta$VCF_NAME %in% 
            plot_tab$samp_name[which(plot_tab$tet_pco_4_mean < 0.25 & 
                                       plot_tab$genepool == 'Atlantic')], 
          STATE]

samp_meta[samp_meta$VCF_NAME %in% 
            plot_tab$samp_name[which(plot_tab$tet_pco_4_mean > 0.5 & 
                                       plot_tab$tet_pco_4_mean < 0.7 &
                                       plot_tab$genepool == 'Atlantic')], 
          STATE]

samp_meta[samp_meta$VCF_NAME %in% 
            plot_tab$samp_name[which(plot_tab$tet_pco_4_mean > 0.7 & 
                                       plot_tab$genepool == 'Gulf')], 
          STATE]
# Gulf Coast (MS, LA, FL)

# Look at dip hi-cor SNPs with low correlation to tet-PCo4 
low_cor_tet_snps <- dipG_v_dipP_m1_pco4_hi_cor_snps[which(abs(
  tetG_v_tetP_cor_m1_pco4[dipG_v_dipP_m1_pco4_hi_cor_snps]) < 0.05)]
length(low_cor_tet_snps)
# 153 of 200 have correlations below 0.05

# What is distrubution of gene pools for HET genotypes in low-cor-tet SNPs
low_cor_tet_het_tabs <- list()
for(i in seq(length(low_cor_tet_snps))){
  low_cor_tet_het_tabs[[i]] <- table(plot_tab$genepool[
    which(dip_mat_1[low_cor_tet_snps[i], c(3:ncol(dip_mat_1))] == 1)])
}

#  What portion of low-cor-tet SNPs are HET in Midwest
low_cor_tet_het_MW_portion <- unlist(lapply(low_cor_tet_het_tabs, function(x)
  x['Midwest']/sum(x)))

# What ploidy is HET in low_cor_tet_snps
low_cor_tet_het_pl_tabs <- list()
for(i in seq(length(low_cor_tet_snps))){
  low_cor_tet_het_pl_tabs[[i]] <- table(plot_tab$ploidy[
    which(dip_mat_1[low_cor_tet_snps[i], c(3:ncol(dip_mat_1))] == 1)])
}
# In low-cor-tet SNPs what portion of HET genotypes are 8X
low_cor_tet_8X_het_portion <- unlist(lapply(low_cor_tet_het_pl_tabs, 
                                            function(x) x['8X']/sum(x)))

# Idea:
# Find low-correlation SNPs with < 50% Midwest and look at genotypes
test_low_cor_inds <- low_cor_tet_snps[which(low_cor_tet_het_MW_portion < 0.5)]

low_cor_tet_het_tabs[which(low_cor_tet_het_MW_portion < 0.5)[4]]
low_cor_tet_8X_het_portion[which(low_cor_tet_het_MW_portion < 0.5)[4]]

table(unlist(dip_mat_1[test_low_cor_inds[4], c(3:ncol(dip_mat_1))]))
table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))]))

table(unlist(dip_mat_1[test_low_cor_inds[4], c(3:ncol(dip_mat_1))[
  which(plot_tab$genepool == 'Atlantic')]]))
table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  which(plot_tab$genepool == 'Atlantic')]]))
table(unlist(dip_mat_1[test_low_cor_inds[4], c(3:ncol(dip_mat_1))[
  which(plot_tab$genepool == 'Midwest')]]))
table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  which(plot_tab$genepool == 'Midwest')]]))

table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  intersect(which(plot_tab$genepool == 'Atlantic'),
            which(plot_tab$ploidy == '4X'))]]))

table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  intersect(which(plot_tab$genepool == 'Gulf'),
            which(plot_tab$ploidy == '4X'))]]))
table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  intersect(which(plot_tab$genepool == 'Gulf'),
            which(plot_tab$ploidy == '8X'))]]))

table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  intersect(which(plot_tab$genepool == 'Midwest'), 
            which(plot_tab$ploidy == '4X'))]]))
table(unlist(tet_mat_1[test_low_cor_inds[4], c(3:ncol(tet_mat_1))[
  intersect(which(plot_tab$genepool == 'Midwest'), 
            which(plot_tab$ploidy == '8X'))]]))

# Observations:
# test_low_cor_inds[1] = 21060:
#   in dip genotypes, almost all Atlantic and Gulf samples are 1; most
#     MW samples are 1, though 37 (all 8X) are 2
#   in tet genotypes, most Atlantic and Gulf samples stay as 1 and most
#     4X MW samples stay as 1 BUT 38 of 52 8X (1 with dip genos) shift to 1.5
#   So, many 8X-MW are 2, and most HET 8X-MW shift to 1.5 while most 4X
#     stay as 1 with tet genotypes
# test_low_cor_inds[2] = 15865:
#   in dip genotypes, ~1/2 Atlantic are HET; 30% (63) of MW are HET, 
#     but of those, almost all are 8X
#   in tet genotypes, almost all ATL HETs shift to 1.5; only 24 of 49 8X-MW
#     shift to 1.5
# test_low_cor_inds[4] = :29086
#   in dip genotypes, only a few ATL are HETs, on het MW are 8X, and half are
#      hets
#   in tet genotypes, almost everything shifts to 1.5; not sure why this
#     disrupts the correlation...

# How many SNPs are HET in more than 50% of samples
het_tally <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))], 1, 
                   function(x) sum(x == 1, na.rm = T))

het_ATL_tally <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))[
  intersect(which(plot_tab$genepool == 'Atlantic'),
            which(plot_tab$ploidy == '4X'))]], 1, 
                       function(x) sum(x == 1, na.rm = T))
length(intersect(which(plot_tab$genepool == 'Atlantic'),
                 which(plot_tab$ploidy == '4X')))
# 264
sum(het_ATL_tally > (264/2))
# 656

het_GULF_tally <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))[
  intersect(which(plot_tab$genepool == 'Gulf'),
            which(plot_tab$ploidy == '4X'))]], 1, 
  function(x) sum(x == 1, na.rm = T))
length(intersect(which(plot_tab$genepool == 'Gulf'),
                 which(plot_tab$ploidy == '4X')))
# 126
sum(het_GULF_tally > 126/2)
# 602

het_MW_tally <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))[
  intersect(which(plot_tab$genepool == 'Midwest'),
            which(plot_tab$ploidy == '4X'))]], 1, 
  function(x) sum(x == 1, na.rm = T))
length(intersect(which(plot_tab$genepool == 'Midwest'),
                 which(plot_tab$ploidy == '4X')))
# 123
sum(het_MW_tally > 123/2)
# 821
sum(het_MW_tally > 86)
# 494

p <- 0.5
q <- 1-p
max_8X_freq <- (4*p^3*q) + (6*p^2*q^2) + (4*p*q^3)

het_MW_8X_tally <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))[
  intersect(which(plot_tab$genepool == 'Midwest'),
            which(plot_tab$ploidy == '8X'))]], 1, 
  function(x) sum(x == 1, na.rm = T))
length(intersect(which(plot_tab$genepool == 'Midwest'),
                 which(plot_tab$ploidy == '8X')))
# 89
sum(het_MW_8X_tally > (89*max_8X_freq))
# 364

excess_het_inds <- c(which(het_ATL_tally > (264/2)), 
                     which(het_GULF_tally > 126/2),
                     which(het_MW_tally > 123/2),
                     which(het_MW_8X_tally > (89*max_8X_freq)))
table(table(excess_het_inds))
#   1   2   3   4 
# 576 283 171 197 

length(unique(excess_het_inds))
# 1227
# so, more than 1% of the SNPs are probably including misassembled regions
#  in at least one of the genepools - this could DEFINITELY affect the
#  pop structure analysis

dip_mat_2 <- dip_mat_1[-excess_het_inds, ]

dip_mat_2_convert <- as.matrix(dip_mat_2[, c(3:ncol(dip_mat_2))])
dip_dist <- dist(t(dip_mat_2_convert), method = 'euclidean', 
                 diag = T, upper = T)
dip_convert_pcoa <- cmdscale(dip_dist, k = 40, eig = T)

tmp_res <- data.frame(dip_convert_pcoa$points)

ggplot(tmp_res, aes(x = X1, y = X9)) + geom_point()

which(tmp_res$X1 > 40)

# How do genotypes change at shared and not-shared hi-correlation SNPs
## shared
table(unlist(dip_mat_1[hi_cor_both_snps[9], c(3:ncol(dip_mat_1))]))
table(unlist(tet_mat_1[hi_cor_both_snps[9], c(3:ncol(tet_mat_1))]))



table(unlist(dip_mat_1[low_cor_tet_snps[4], c(3:ncol(dip_mat_1))]))
table(unlist(tet_mat_1[low_cor_tet_snps[4], c(3:ncol(tet_mat_1))]))



dipG_v_tetP_cor_m1_pco4 <- apply(dip_mat_1[, c(3:ncol(dip_mat_1))], 1, cor, 
                                 y = plot_tab$tet_pco_4_mean)
dipG_v_tetP_m1_pco4_bottom_snps <- order(dipG_v_tetP_cor_m1_pco4)[1:100]
dipG_v_tetP_m1_pco4_top_snps <- order(dipG_v_tetP_cor_m1_pco4, 
                                      decreasing = T)[1:100]



length(intersect(dipG_v_dipP_m1_pco4_bottom_snps, 
                 dipG_v_tetP_m1_pco4_bottom_snps))
# 6

length(intersect(dipG_v_dipP_m1_pco4_top_snps, 
                 dipG_v_tetP_m1_pco4_top_snps))
# 33

