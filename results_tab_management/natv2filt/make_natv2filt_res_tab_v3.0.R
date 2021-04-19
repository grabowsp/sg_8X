# Add small-scale population separations to genepool-level files
#  Plot results to make sure all groups are correct before saving
### LOAD PACKAGES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

misfit_pca_res_file <- '/Users/grabowsk/data/sg_8X_analysis/GW.100kSNPs.tetrasomic.CDS.natv2filt.MISFIT.genlight.PCAresults.rds'
misfit_res <- readRDS(misfit_pca_res_file)

# results tables
mw_res_tab_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                         'natv2filt_MW_res_tab_v1.0.txt', sep = '')
mw_res_tab <- fread(mw_res_tab_file)

gulf_res_tab_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                         'natv2filt_GULF_res_tab_v1.0.txt', sep = '')
gulf_res_tab <- fread(gulf_res_tab_file)

atl_res_tab_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                         'natv2filt_ATL_res_tab_v1.0.txt', sep = '')
atl_res_tab <- fread(atl_res_tab_file)

### SET OUTPUTS ###
mw_res_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                         'natv2filt_MW_res_tab_v2.0.txt', sep = '')

gulf_res_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                           'natv2filt_GULF_res_tab_v2.0.txt', sep = '')

atl_res_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                          'natv2filt_ATL_res_tab_v2.0.txt', sep = '')

mf_res_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                    'natv2filt_MISFIT_res_tab_v2.0.txt', sep = '')
############

## Midwest

grp_mw2_1 <- which(mw_res_tab$mw_tot_PC1 < -5 & mw_res_tab$mw_tot_PC2 < -7.5)
grp_mw2_2 <- which(mw_res_tab$mw_tot_PC1 > -5 & mw_res_tab$mw_tot_PC1 < -2.5 &
                     mw_res_tab$mw_tot_PC2 > 0)
grp_mw2_3 <- which(mw_res_tab$mw_tot_PC1 < -9.5 & mw_res_tab$mw_tot_PC2 > 4)
grp_mw2_4 <- which(mw_res_tab$mw_tot_PC1 < -5 & mw_res_tab$mw_tot_PC2 > -7 &
                     mw_res_tab$mw_tot_PC2 < 4) # inspired by k=3 admix results
grp_mw2_5 <- which(mw_res_tab$mw_tot_PC1 > -10 & mw_res_tab$mw_tot_PC1 < -5 &
                     mw_res_tab$mw_tot_PC2 > 5)
grp_mw3_1 <- which(mw_res_tab$mw_tot_PC1 > -2.5 & mw_res_tab$mw_tot_PC1 < 0)
grp_mw3_2 <- which(mw_res_tab$mw_tot_PC1 > 0 & mw_res_tab$mw_tot_PC1 < 5
                & mw_res_tab$mw_tot_PC3 > 0)
grp_mw3_3 <- which(mw_res_tab$mw_tot_PC1 > 5 & mw_res_tab$mw_tot_PC1 < 9 
                & mw_res_tab$mw_tot_PC3 > 0)
grp_mw3_4 <- which(mw_res_tab$mw_tot_PC1 > 2.5 & mw_res_tab$mw_tot_PC1 < 8
                & mw_res_tab$mw_tot_PC3 < 0)
grp_mw4_1 <- which(mw_res_tab$mw_tot_PC1 > 10 & mw_res_tab$mw_tot_PC4 > 0)
grp_mw4_2 <- which(mw_res_tab$mw_tot_PC1 < 10 & mw_res_tab$mw_tot_PC1 > 7.5 &
                      mw_res_tab$mw_tot_PC4 > 0)
grp_mw5_1 <- which(mw_res_tab$mw_tot_PC5 < -10)
grp_mw5_2 <- which(mw_res_tab$mw_tot_PC1 > 5 & mw_res_tab$mw_tot_PC5 < -5 &
                     mw_res_tab$mw_tot_PC5 > -10)
grp_mw5_3 <- which(mw_res_tab$mw_tot_PC1 > 8.5 & mw_res_tab$mw_tot_PC1 < 11 &
                     mw_res_tab$mw_tot_PC5 > 3.5)

tot_mw_grps <- c(grp_mw2_1, grp_mw2_2, grp_mw2_3, grp_mw2_4, grp_mw2_5,
                 grp_mw3_1, grp_mw3_2, grp_mw3_3, grp_mw3_4, 
                 grp_mw4_1, grp_mw4_2, grp_mw5_1, grp_mw5_2, grp_mw5_3)
sum(duplicated(tot_mw_grps))
#[1] 0

mw_res_tab[, grps := 'NA']
mw_res_tab[grp_mw2_1, grps := 'MW_01'] # West/Cosmo 8X
mw_res_tab[grp_mw2_2, grps := 'MW_02']
mw_res_tab[grp_mw2_3, grps := 'MW_03'] # East 8X
mw_res_tab[grp_mw2_4, grps := 'MW_04']
mw_res_tab[grp_mw2_5, grps := 'MW_05']
mw_res_tab[grp_mw3_1, grps := 'MW_06']
mw_res_tab[grp_mw3_2, grps := 'MW_07']
mw_res_tab[grp_mw3_3, grps := 'MW_08']
mw_res_tab[grp_mw3_4, grps := 'MW_09']
mw_res_tab[grp_mw4_1, grps := 'MW_10']
mw_res_tab[grp_mw4_2, grps := 'MW_11']
mw_res_tab[grp_mw5_1, grps := 'MW_12']
mw_res_tab[grp_mw5_2, grps := 'MW_13']
mw_res_tab[grp_mw5_3, grps := 'MW_14']

summary(mw_res_tab$mw_tot_PC1[which(mw_res_tab$grps == 'NA')])
# all MW samples assigned to a gene pool

mw_res_tab$grps <- as.factor(mw_res_tab$grps)

levels(mw_res_tab$grps) <- sort(levels(mw_res_tab$grps))

gg_mw_tmp <- ggplot(mw_res_tab, aes(x = mw_tot_PC1, y = mw_tot_PC2, 
                                    color = grps, 
                                        shape = ploidy)) +
  geom_point()

gg_mw_tmp

fwrite(mw_res_tab, file = mw_res_out, sep = '\t')

# MW ADMIXTURE K=3 Interpretations
# mw_tot_k3_pop_1 = MW-8X-West/Cosmopolitan + related 4X
# mw_tot_k3_pop_2 = most MW-4X
# mw_tot_k3_pop_3 = MW-8X-East; not present in any 4X samples

# MW Group Interpretations
# MW_01 = MW-8X-West/Cosmopolitan
## 0.86 to 1.0 pop_1
## AR (15) CO (2) IL (1) MI (1) MN (4) NE (5) NM (2) OH (1) OK (1) TX (3)

# MW_02 = MW-8X-Mixture; mixture of all 3 gene pools
## 0.27 pop_1, 0.28 pop_2; 0.45 pop_3
## IN (1)

# MW_03 = main MW-8X-East group
## 0.91 to 1.0 pop_3; 0.0 to 0.09 pop_1
## PA (30) NY (14) WV (2) MA (2) VA (2) KY (1) OH (1)

# MW_04 = mainly MW-8X-East, but some 8X-West
## 0.32 to 0.80 pop_3; 0.14 to 0.68 pop_1
## PA (3) MO (1) NE (1) NY (2) ME (1) IA (1)

# MW_05 = MW-8X-East with some 4X
## 0.8 to 0.91 pop_3; 0.09 to 0.19 pop_2
## OH (6)

# MW_06 = MW-4X closely related to MW-8X West; 
## 0.8 to 1.0 pop_1; 0.0 to 0.21 pop_2
## MO (6) MI (2)

# MW_07 = MW-4X close to MW-8X West
## 0.52 to 0.69 pop_1; 0.31 to 0.48 pop_2
## MO (5) AK (4)

# MW_08 = MW-4X distant to MW-8X West
## 0.17 to 0.35 pop_1; 0.64 to 0.83 pop_2
## MO (18)

# MW_09 = MW-4X distant to MW-8X West
## 0.19 to 0.31 pop_1; 0.69 to 0.81 pop_2
## MO (2) NE (2) MN (2) WI (3) NY (1)

# MW_10 = MW-4X-GreatLakes
## pure pop_2
## MI (23) IL (16) WI (6) IN (3) 

# MW_11
## 0.94 pop_2
## IN (1)

# MW_12 = MW-4X-East1
## pure pop_2 
## OH (10)

# MW_13 = MW-4X-East2
## pure pop_2
## OH (6)

# MW_14 = MW-4X-UpperMW
## 0.85 to 1.0 pop_2; 0.0 to 0.15 pop_1
## MN (10) WI (2) IA (1) NY (1)

# NOTE: MW_10, MW_12, and MW_13 all have the same patterns from the
## ADMIXTURE results, but show some distinct clustering in PCA, 
## though only at the higher PCs (4 and 5)
# I could be conviced to group all of MW_4_X and MW_5_X together if we want to

###

# Gulf subgropus
grp_g2_1 <- which(gulf_res_tab$gulf_tot_PC2 > 15)
grp_g2_2 <- which(gulf_res_tab$gulf_tot_PC1 > 16 & 
                    gulf_res_tab$gulf_tot_PC2 < 15)
grp_g2_3 <- which(gulf_res_tab$gulf_tot_PC1 < -17)
grp_g3_1 <- which(gulf_res_tab$gulf_tot_PC3 > 10)
grp_g3_2 <- which(gulf_res_tab$gulf_tot_PC3 > 0 & 
                    gulf_res_tab$gulf_tot_PC3 < 10 & 
                    gulf_res_tab$gulf_tot_PC1 > -5 & 
                    gulf_res_tab$gulf_tot_PC1 < 0)
grp_g3_3 <- which(gulf_res_tab$gulf_tot_PC3 < -10 & 
                    gulf_res_tab$gulf_tot_PC3 > -16 &
                    gulf_res_tab$gulf_tot_PC1 > 5 & 
                    gulf_res_tab$gulf_tot_PC1 < 10)
grp_g3_4 <- which(gulf_res_tab$gulf_tot_PC1 > 8 & 
                    gulf_res_tab$gulf_tot_PC1 < 12 &
                    gulf_res_tab$gulf_tot_PC3 < -15)
grp_g3_5 <- which(gulf_res_tab$gulf_tot_PC1 > -7 & 
                    gulf_res_tab$gulf_tot_PC1 < 0 &
                    gulf_res_tab$gulf_tot_PC3 < -20)
grp_g3_6 <- which(gulf_res_tab$gulf_tot_PC1 > -17.5 & 
                    gulf_res_tab$gulf_tot_PC1 < -6)
grp_g3_extra <- which(gulf_res_tab$gulf_tot_PC1 > 0 & 
                    gulf_res_tab$gulf_tot_PC1 < 16 &
                    gulf_res_tab$gulf_tot_PC3 > -10 & 
                    gulf_res_tab$gulf_tot_PC3 < 5)
# split up grp_g3_extra using ADMIXTURE results
grp_AM_b <- intersect(grp_g3_extra, 
                      which(gulf_res_tab$gulf_tot_k3_pop_2 > 0.4 & 
                               gulf_res_tab$gulf_tot_k3_pop_3 > 0.4))
grp_AM_1 <- intersect(grp_g3_extra, 
                     which(gulf_res_tab$gulf_tot_k3_pop_2 > 0.4 & 
                             gulf_res_tab$gulf_tot_k3_pop_3 < 0.4 & 
                             gulf_res_tab$gulf_tot_k3_pop_1 > 0.01 &
                             gulf_res_tab$gulf_tot_PC3 > 0))
grp_AM_2 <- intersect(grp_g3_extra, 
                      which(gulf_res_tab$gulf_tot_k3_pop_2 > 0.4 & 
                              gulf_res_tab$gulf_tot_k3_pop_3 < 0.4 & 
                              gulf_res_tab$gulf_tot_k3_pop_1 > 0.01 &
                              gulf_res_tab$gulf_tot_PC3 < 0))
grp_AM_3 <- intersect(grp_g3_extra, 
                      which(gulf_res_tab$gulf_tot_k3_pop_2 > 0.4 & 
                              gulf_res_tab$gulf_tot_k3_pop_3 < 0.4 & 
                              gulf_res_tab$gulf_tot_k3_pop_1 < 0.01))
grp_AM_4 <- intersect(grp_g3_extra, 
                      which(gulf_res_tab$gulf_tot_k3_pop_2 < 0.4 & 
                              gulf_res_tab$gulf_tot_k3_pop_3 > 0.4))

gulf_allgrps <- c(grp_g2_1, grp_g2_2, grp_g2_3, grp_g3_1, grp_g3_2, grp_g3_3,
                  grp_g3_4, grp_g3_5, grp_g3_6, grp_AM_b, grp_AM_1, grp_AM_2,
                  grp_AM_3, grp_AM_4)
sum(duplicated(gulf_allgrps))
#[1] 0

gulf_res_tab[, grps := 'NA']
gulf_res_tab[grp_g2_1, grps := 'GULF_01']
gulf_res_tab[grp_g2_2, grps := 'GULF_02']
gulf_res_tab[grp_g2_3, grps := 'GULF_03']
gulf_res_tab[grp_g3_1, grps := 'GULF_04']
gulf_res_tab[grp_g3_2, grps := 'GULF_05']
gulf_res_tab[grp_g3_3, grps := 'GULF_06']
gulf_res_tab[grp_g3_4, grps := 'GULF_07']
gulf_res_tab[grp_g3_5, grps := 'GULF_08']
gulf_res_tab[grp_g3_6, grps := 'GULF_09']
gulf_res_tab[grp_AM_b, grps := 'GULF_10']
gulf_res_tab[grp_AM_1, grps := 'GULF_11']
gulf_res_tab[grp_AM_2, grps := 'GULF_12']
gulf_res_tab[grp_AM_3, grps := 'GULF_13']
gulf_res_tab[grp_AM_4, grps := 'GULF_14']

which(gulf_res_tab$gulf_tot_PC1 > -200 & gulf_res_tab$grps == 'NA')
# all Gulf samples accounted for

gulf_res_tab$grps <- as.factor(gulf_res_tab$grps)

levels(gulf_res_tab$grps) <- sort(levels(gulf_res_tab$grps))

gg_gulf_tmp <- ggplot(gulf_res_tab, aes(x = gulf_tot_PC1, y = gulf_tot_PC2, 
                                        color = grps, 
                                      shape = ploidy)) +
  geom_point()

gg_gulf_tmp

fwrite(gulf_res_tab, file = gulf_res_out, sep = '\t')

###### GULF Total K=3 Explanation
# Pop 1 = Gulf Coast 1
# Pop 2 = Gulf Coast 2
# Pop 3 = TX/MX

###### GULF subgroup explanations
# GULF_01 - pure TX/MX (pop_3)
## MX (11)

# GULF_02 - pure TX/MX (pop_3)
## TX(8)

# GULF_03 - Pure GC 1 (pop_1)
## MS(20) FL (8) LA (7)

# GULF_04 - Pure GC 2 (pop_2)
## AR (12) LA (3)

# GULF_05 - Mainly GC 2 (pop_2) with a little GC 1
## 0.82 to 0.96 pop_2; 0.04 to 0.18 pop_1
## LA(6) TX (1)

# GULF_06 - Mainly GC 2 (pop_2) with a little TX/MX
## 0.78 to 0.92 pop_2; 0.08 to 0.22 pop_3
## TX (11)

# GULF_07 - Mainly GC 2 (pop_2) with a little TX/MX
## 0.74 to 0.81 pop_2; 0.19 to 0.26 pop_3
## TX (4)

# GULF_08 - Mainly GC 2 (pop_2) with some GC 1
## 0.64 to 0.71 pop_2; 0.29 to 0.36 pop_1
## FL (3)

# GULF_09 - Majority Gulf Coast 1 (pop_1) with some GC 2 and less TX/MX
## 0.49 to 0.96 pop_1; 0.04 to 0.51 pop_2; 0.0 to 0.09 pop_3
## LA (9) FL (4)

# GULF_10 - Majority GC 2 (pop_2) with TX/MX and a bit GC 1
## 0.72 to 0.86 pop_2; 0.07 to 0.23 pop_3; 0.02 to 0.06 pop_1
## TX (4)

# GULF_11 - Majority GC 2 (pop 2) with GC 1 and/or TX/MX
## 0.71 to 0.72 pop_2; 0.13 to 0.20 pop_3, 0.07 to 0.17 pop_1
## FL (3)

# GULF_12 - Majority GC 2 (pop_2) with some TX/MX
## 0.63 to 0.89 pop_2, 0.11 to 0.37 pop_3
## TX (18)

# GULF_13 - TX/MX (pop_3) with some GC 2 (pop_2)
## 0.63 to 0.80 pop_3; 0.20 to 0.37 pop_2
## TX (18) SC (2)

# GULF_14 - roughly half TX/MX (pop_3) and GC 2 (pop_2) 
## 0.46 to 0.59 pop_3; 0.41 to 0.54 pop_2
## TX (5)

#########

### Atlantic sub-groups
grp_a2_2 <- which(atl_res_tab$atl_tot_PC2 > 15)
grp_a2_3 <- which(atl_res_tab$atl_tot_PC2 > 9 & 
                    atl_res_tab$atl_tot_PC2 < 11)
grp_a2_4 <- which(atl_res_tab$atl_tot_PC1 > 5 & 
                    atl_res_tab$atl_tot_PC2 > 2.5 &
                    atl_res_tab$atl_tot_PC2 < 9)
grp_a3_1 <- which(atl_res_tab$atl_tot_PC3 < -7)
grp_a3_2 <- which(atl_res_tab$atl_tot_PC1 > -7.5 & 
                    atl_res_tab$atl_tot_PC1 < -2.5 &
                    atl_res_tab$atl_tot_PC3 > 0 & 
                    atl_res_tab$atl_tot_PC3 < 6)
grp_a3_3 <- which((atl_res_tab$atl_tot_PC3 > 0 & 
                     atl_res_tab$atl_tot_PC2 < -5) |
                    atl_res_tab$atl_tot_PC3 > 10)
grp_a3_4 <- which(atl_res_tab$atl_tot_PC2 > 11 & atl_res_tab$atl_tot_PC2 < 15)
grp_a4_1 <- which(atl_res_tab$atl_tot_PC4 < -5 & 
                    atl_res_tab$atl_tot_PC1 > -15 &
                    atl_res_tab$atl_tot_PC1 < -5)
grp_a4_2 <- which(atl_res_tab$atl_tot_PC1 > -11 & 
                    atl_res_tab$atl_tot_PC1 < -8 & 
                    atl_res_tab$atl_tot_PC4 > 0 & 
                    atl_res_tab$atl_tot_PC4 < 2.5)
grp_a4_3 <- which((atl_res_tab$atl_tot_PC1 < -10 & 
                    atl_res_tab$atl_tot_PC1 > -13 & 
                    atl_res_tab$atl_tot_PC4 > 2.5) | 
                    (atl_res_tab$atl_tot_PC1 < -11 & 
                       atl_res_tab$atl_tot_PC1 > -13 & 
                       atl_res_tab$atl_tot_PC4 > 1.5))
grp_a2_1 <- setdiff(which(atl_res_tab$atl_tot_PC1 < -7.5), 
                    c(grp_a4_1, grp_a4_2, grp_a4_3))
grp_a2_5 <- which(atl_res_tab$atl_tot_PC1 > 4 &
                    atl_res_tab$atl_tot_PC1 < 8 &
                    atl_res_tab$atl_tot_PC2 < 2.5 &
                    atl_res_tab$atl_tot_PC3 > -5.5 &
                    atl_res_tab$atl_tot_k4_pop_3 > 0.01 &
                    atl_res_tab$atl_tot_k4_pop_4 < 0.3)
grp_a2_6 <- which(atl_res_tab$atl_tot_PC1 > 4 &
                    atl_res_tab$atl_tot_PC1 < 8 &
                    atl_res_tab$atl_tot_PC2 < 2.5 &
                    atl_res_tab$atl_tot_PC3 > -5.5 &
                    atl_res_tab$atl_tot_k4_pop_3 < 0.01)
grp_a2_7 <- which(atl_res_tab$atl_tot_PC1 > 4 &
                    atl_res_tab$atl_tot_PC1 < 8 &
                    atl_res_tab$atl_tot_PC2 < 2.5 &
                    atl_res_tab$atl_tot_PC3 > -5.5 &
                    atl_res_tab$atl_tot_k4_pop_3 > 0.01 &
                    atl_res_tab$atl_tot_k4_pop_4 > 0.3)

grp_a_all <- c(grp_a2_1, grp_a2_2, grp_a2_3, grp_a2_4, 
               grp_a2_5, grp_a2_6, grp_a2_7, 
               grp_a3_1, grp_a3_2, grp_a3_3, grp_a3_4, 
               grp_a4_1, grp_a4_2, grp_a4_3)

sum(duplicated(grp_a_all))
# 0

# grp_a_rest <- setdiff(seq(nrow(atl_res_tab)), grp_a_all)

atl_res_tab[, grps := 'NA']
atl_res_tab[grp_a2_1, grps := 'ATL_01']
atl_res_tab[grp_a2_2, grps := 'ATL_02']
atl_res_tab[grp_a2_3, grps := 'ATL_03']
atl_res_tab[grp_a2_4, grps := 'ATL_04']
atl_res_tab[grp_a2_5, grps := 'ATL_05']
atl_res_tab[grp_a2_6, grps := 'ATL_06']
atl_res_tab[grp_a2_7, grps := 'ATL_07']
atl_res_tab[grp_a3_1, grps := 'ATL_08']
atl_res_tab[grp_a3_2, grps := 'ATL_09']
atl_res_tab[grp_a3_3, grps := 'ATL_10']
atl_res_tab[grp_a3_4, grps := 'ATL_11']
atl_res_tab[grp_a4_1, grps := 'ATL_12']
atl_res_tab[grp_a4_2, grps := 'ATL_13']
atl_res_tab[grp_a4_3, grps := 'ATL_14']

which(atl_res_tab$atl_tot_PC1 > -200 & atl_res_tab$grps == 'NA')
# 0

atl_res_tab$grps <- as.factor(atl_res_tab$grps)
levels(atl_res_tab$grps) <- sort(levels(atl_res_tab$grps))

gg_atl_tmp <- ggplot(atl_res_tab, aes(x = atl_tot_PC1, y = atl_tot_PC2, 
                                      color = grps, 
                                      shape = ploidy)) +
  geom_point()

gg_atl_tmp

fwrite(atl_res_tab, file = atl_res_out, sep = '\t')

### K=4 Explanation
# Pop_1 = NY & MidAtlantic
## DE (13) ML (2) NJ (24) NY (32) NC (13) VA (4)
# Pop_2 = New England and NY
## CT (8) ME (3) MA (9) NH (3) NY (42) RI (20)
# Pop_3 = FL and SE Atlantic Coast
## FL (21) GA (3) NC (16) SC (19)
# Pop_4 = MidAtlantic
# MD (24) VA (15)

### Subgroup Explanation
# ATl_01 - New England and NY
## 0.75 to 1.0 pop_2
## CT (8) ME (3) MA (3) NH (3) NY (20) RI (20)

# ATL_02 - FL-1
## pure pop_3
## FL (9)

# ATL_03 - SE Atlantic Coast
## Pure pop_3
## FL (1) GA (3) SC (2)

# ATL_04 - Carolina Coast
## 0.63 to 0.97 pop_3; 0.03 to 0.24 pop_4
## NC (11) SC (17)

# ATL_05 - NC
## 0.28 to 0.61 pop_1, 0.26 to 0.62 pop_3; 0.08 to 0.21 pop_4
## NC (25)

# ATL_06 - VA
## 0.53 to 0.68 pop_4; 0.32 to 0.46 pop_1, 
## VA (8)

# ATL_07 - Southern MidAtlantic
## 0.37 to 0.76 pop_4; 0.16 to 0.53 pop_1; 0.05 to 0.29 pop_3
## MD (1) NC (3) VA (2)

# ATL_08 - MidAtlantic
## 0.82 to 1.0 pop_4
## MD (23) VA (7)

# ATL_09 - NY-1
## 0.58 to 0.70 pop_1, 0.30 to 0.42 pop_2
## NY (5)

# ATL_10 - Northern MidAtlantic
## 0.73 to 1.0 pop_1; 0.0 to 0.27 pop_4
## DE (13) MD (2) NJ (24) NY (27) VA (3)

# ATL_11 - FL-2
## pure pop_3
## FL (11)

# ATL_12 - MA
## 0.76 to 0.86 pop_2
## MA (6)

# ATL_13 - NY-2
## 0.64 to 0.71 pop_2 
## NY (8)

# ATL_14 - NY-3
## 0.80 to 0.99 pop_2
## NY (14)

########

## Misfit splits

# Generate Misfit table
misfit_inds <- c()
for(i in rownames(misfit_res$scores)){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  misfit_inds <- c(misfit_inds, tmp_ind)
}

misfit_tab <- data.table(samp_names = rownames(misfit_res$scores), 
                          ploidy = samp_meta$NQUIRE_PLOIDY[misfit_inds],
                          PC1 = misfit_res$scores[,1], 
                          PC2 = misfit_res$scores[,2],
                          PC3 = misfit_res$scores[,3],
                          PC4 = misfit_res$scores[,4],
                          PC5 = misfit_res$scores[,5],
                          PC6 = misfit_res$scores[,6])

gg_misfit_tmp <- ggplot(misfit_tab, aes(x = PC1, y = PC2)) +
  geom_point(color = misfit_tab$ploidy)

gg_misfit_tmp

## Split misfit groups

grp_mf1_1 <- which(misfit_tab$PC1 > 30)
grp_mf1_2 <- which(misfit_tab$PC1 > 20 & misfit_tab$PC1 < 30)
grp_mf1_3 <- which(misfit_tab$PC1 > -2.5 & misfit_tab$PC1 < 1)

grp_mf2_1 <- which(misfit_tab$PC2 < -10)
grp_mf2_2 <- which(misfit_tab$PC2 > 10 & misfit_tab$PC2 < 15)
grp_mf2_3 <- which(misfit_tab$PC1 > 10 & misfit_tab$PC1 < 15 &
                     misfit_tab$PC2 > 0)
#grp_mf3_1 <- which(misfit_tab$PC3 > 20)
grp_mf3_2 <- which(misfit_tab$PC1 < -1.5 & misfit_tab$PC3 > -4 & 
                     misfit_tab$PC3 < 0)
grp_mf3_3 <- which(misfit_tab$PC1 < -5 & misfit_tab$PC3 < -4)
#grp_mf6_1 <- which(misfit_tab$PC1 < -7 & misfit_tab$PC1 < -2 &
#                     misfit_tab$PC6 < -2)
grp_tmp_all <- c(grp_mf1_1, grp_mf1_2, grp_mf1_3, grp_mf2_1, grp_mf2_2, 
                 grp_mf2_3, grp_mf3_2, grp_mf3_3)
grp_mf_rest <- setdiff(seq(nrow(misfit_tab)), grp_tmp_all)

length(grp_mf_rest)
# 0

misfit_tab[, grps := 'NA']
misfit_tab[grp_mf1_1, grps := 'MF_01']
misfit_tab[grp_mf1_2, grps := 'MF_02']
misfit_tab[grp_mf1_3, grps := 'MF_03']
misfit_tab[grp_mf2_1, grps := 'MF_04']
misfit_tab[grp_mf2_2, grps := 'MF_05']
misfit_tab[grp_mf2_3, grps := 'MF_06']
#misfit_tab[grp_mf3_1, grps := 'MF_07']
misfit_tab[grp_mf3_2, grps := 'MF_07']
misfit_tab[grp_mf3_3, grps := 'MF_08']
#misfit_tab[grp_mf_rest, grps := 'MF_REST']

misfit_tab$grps <- as.factor(misfit_tab$grps)
levels(misfit_tab$grps) <- sort(levels(misfit_tab$grps))

gg_misfit_tmp <- ggplot(misfit_tab, aes(x = PC1, y = PC2, 
                                        color = grps,
                                        shape = ploidy)) +
  geom_point()

gg_misfit_tmp

fwrite(misfit_tab, file = mf_res_out, sep = '\t')

# Explanation of Misfit groups

# MF_01 - FL-4X
## 4X
## FL (2)

# MF_02 - NC-4X
## 4X
## NC (2) 

# MF_03 - GA-8X-1
## 8X
## GA (2) NC (1)

# MF_04 - NY/NJ
## 4X (6) 8X (1)
## NJ (3) NY (4)

# MF_05 - GA-8X-2
## 8X (3)
## GA (3)

# MF_06 - FL-8X
## 8X (1)
## FL (1)

# MF_07 - GA-8X-3
## 8X (7)
## GA (5) IL (1) TN (1)

# MF_08 - MS-8X
## 8X (11) 6X (1) 
## AR (1) FL (1) MS (10)
