
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

##
# Notes on splitting up samples
# First grouping is defnitely to split into 3 - maybe use ADMIX MW_tot K=3
### results - I think this will break up the 8X into 2
# For Splitting up 4X:
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
mw_res_tab[grp_mw2_1, grps := 'MW_2_1'] # West/Cosmo 8X
mw_res_tab[grp_mw2_2, grps := 'MW_2_2']
mw_res_tab[grp_mw2_3, grps := 'MW_2_3'] # East 8X
mw_res_tab[grp_mw2_4, grps := 'MW_2_4']
mw_res_tab[grp_mw2_5, grps := 'MW_2_5']
mw_res_tab[grp_mw3_1, grps := 'MW_3_1']
mw_res_tab[grp_mw3_2, grps := 'MW_3_2']
mw_res_tab[grp_mw3_3, grps := 'MW_3_3']
mw_res_tab[grp_mw3_4, grps := 'MW_3_4']
mw_res_tab[grp_mw4_1, grps := 'MW_4_1']
mw_res_tab[grp_mw4_2, grps := 'MW_4_2']
mw_res_tab[grp_mw5_1, grps := 'MW_5_1']
mw_res_tab[grp_mw5_2, grps := 'MW_5_2']
mw_res_tab[grp_mw5_3, grps := 'MW_5_3']

summary(mw_res_tab$mw_tot_PC1[which(mw_res_tab$grps == 'NA')])
# all MW samples assigned to a gene pool

mw_res_tab$grps <- as.factor(mw_res_tab$grps)

levels(mw_res_tab$grps) <- sort(levels(mw_res_tab$grps))

gg_mw_tmp <- ggplot(mw_res_tab, aes(x = mw_tot_PC1, y = mw_tot_PC2, 
                                    color = grps, 
                                        shape = ploidy)) +
  geom_point()

gg_mw_tmp

summary(mw_res_tab$mw_main_k2_pop_2[mw_res_tab$grps == 'MW_2_1'])

# MW-main K=2
for(pop_name in sort(unique(mw_res_tab$grps))){
  print(pop_name)
  print(summary(mw_res_tab$mw_main_k2_pop_2[mw_res_tab$grps == pop_name]))
}

summary(mw_res_tab$mw_main_k2_pop_2[mw_res_tab$grps == 'MW_2_1'])
table(mw_res_tab$grps[which(mw_res_tab$mw_main_k2_pop_1 > 0.5)])
table(mw_res_tab$grps[which(mw_res_tab$mw_main_k2_pop_2 > 0.5)])

# mw_main_k2_pop_1 > 0.5 = 8X samples - MW_2_X
# mw_main_k2_pop_2 > 0.5 = 4x samples (and 1 8X) - remaining samples

# MW-main K=3
summary(mw_res_tab$mw_main_k3_pop_2[mw_res_tab$grps == 'MW_4_2'])

for(pop_name in sort(unique(mw_res_tab$grps))){
  print(pop_name)
  print(summary(mw_res_tab$mw_main_k3_pop_2[mw_res_tab$grps == pop_name]))
}

table(mw_res_tab$grps[which(mw_res_tab$mw_main_k3_pop_1 > 0.5)])
table(mw_res_tab$grps[which(mw_res_tab$mw_main_k3_pop_2 > 0.5)])
table(mw_res_tab$grps[which(mw_res_tab$mw_main_k3_pop_3 > 0.5)])

# mw_main_k3_pop_1 > 0.5 is MW_2_3, MW_2_4, MW_2_5 (and MW_2_2 = .509)
# main_k3_pop_1 = 8X east
## MW_2_3 (PA + others) - almost pure pop_1
## MW_2_4 (MO, NE, NY) - 0.37 to 0.8
## MW_2_5 (OH) - 0.8 to 0.91

# mw_main_k3_pop_2 > 0.5 is 4X groups except MW_3_1 and MW_3_2

# mw_main_k3_pop_3 > 0.5 is MW_2_1, a few MW_2_4, and MW_3_1 and MW_3_2

# MW_3_1 is "pure" mw_main_k3_pop_3
# MW_3_2 is 0.6 to 0.8 mw_main_k3_pop_3, the rest is ...pop_2 (4X)
# MW_3_3 is 0.27 to 0.45 ..pop_3; rest pop_2
# MW_3_4 is 0.3 to 0.47 pop_3; rest pop_2
# MW_5_3 is 0.0 to 0.26 pop_3;

table(samp_meta$STATE[samp_meta$VCF_NAME %in% mw_res_tab$samp_name[which(
  mw_res_tab$grps == 'MW_2_5')]])
# MW_3_1 is Missouri (6) and Michigan (2) samples
# MW_3_2 is Arkansas (3) and Missouri (5)
# MW_3_3 is Missouri (18)
# MW_3_4 is all over - MO (2) MN (2) NE (2) NY (1) WI (3)
# MW_5_3 is MN (10) IA (1) NY (1) and WI (2)
# MW_2_1 (8X) is the western/cosmopolitan 8X group
 
for(pop_name in unique(mw_res_tab$grps)){
  print(pop_name)
  print(table(samp_meta$STATE[samp_meta$PLANT_ID %in% 
                                mw_res_tab$samp_name[which(
                                  mw_res_tab$grps == pop_name)]]))
}

# MW-tot k=2
table(mw_res_tab$grps[which(mw_res_tab$mw_tot_k2_pop_1 > 0.5)])
table(mw_res_tab$grps[which(mw_res_tab$mw_tot_k2_pop_2 > 0.5)])

for(pop_name in sort(unique(mw_res_tab$grps))){
  print(pop_name)
  print(summary(mw_res_tab$mw_tot_k2_pop_2[mw_res_tab$grps == pop_name]))
}

# mw_tot_k2_pop_1 = 4X (MW_3_X, MW_4_X, MW_5_X)
# .._pop_2 = 8X (MW_2_X)

# MW_2_1 (8X West) is .73 to 1 pop_2
# MW_2_5 (8X OH) is 0.8 to 0.91 pop_2
# MW_3_1 is 0.34 to 0.44 pop_2
# MW_3_2 is 0.18 to 0.42 pop_2
# MW_3_3 and MW_3_4 0 to less than 0.1 pop_2

# MW-tot k=3
table(mw_res_tab$grps[which(mw_res_tab$mw_tot_k3_pop_1 > 0.5)])
# MW_2_1, a few MW_2_4 (3), MW_3_1, MW_3_2
table(mw_res_tab$grps[which(mw_res_tab$mw_tot_k3_pop_2 > 0.5)])
# Many of the 4X grps, but not MW_3_1 or MW_3_2
table(mw_res_tab$grps[which(mw_res_tab$mw_tot_k3_pop_3 > 0.5)])
# MW_2_3 (8X East) only

for(pop_name in sort(unique(mw_res_tab$grps))){
  print(pop_name)
  print(summary(mw_res_tab$mw_tot_k3_pop_3[mw_res_tab$grps == pop_name]))
}

for(pop_name in sort(unique(mw_res_tab$grps))){
  print(pop_name)
  print(table(samp_meta$STATE[
    samp_meta$VCF_NAME %in% mw_res_tab$samp_name[which(
      mw_res_tab$grps == pop_name)]]))
}

test_inds <- intersect(which(mw_res_tab$grp == 'MW_2_3'), 
          which(mw_res_tab$mw_tot_k3_pop_1 > 0.1))

mw_res_tab$mw_tot_PC2[test_inds]

# MW ADMIXTURE K=3 Interpretations
# mw_tot_k3_pop_1 = MW-8X-West/Cosmopolitan + related 4X
# mw_tot_k3_pop_2 = most MW-4X
# mw_tot_k3_pop_3 = MW-8X-East; not present in any 4X samples

# MW Group Interpretations
# MW_2_1 = MW-8X-West/Cosmopolitan
## 0.86 to 1.0 pop_1
## AR (15) CO (2) IL (1) MI (1) MN (4) NE (5) NM (2) OH (1) OK (1) TX (3)

# MW_2_2 = MW-8X-Mixture; mixture of all 3 gene pools
## 0.27 pop_1, 0.28 pop_2; 0.45 pop_3
## IN (1)

# MW_2_3 = main MW-8X-East group
## 0.91 to 1.0 pop_3; 0.0 to 0.09 pop_1
## PA (30) NY (14) WV (2) MA (2) VA (2) KY (1) OH (1)

# MW_2_4 = mainly MW-8X-East, but some 8X-West
## 0.32 to 0.80 pop_3; 0.14 to 0.68 pop_1
## PA (3) MO (1) NE (1) NY (2) ME (1) IA (1)

# MW_2_5 = MW-8X-East with some 4X
## 0.8 to 0.91 pop_3; 0.09 to 0.19 pop_2
## OH (6)

# MW_3_1 = MW-4X closely related to MW-8X West; 
## 0.8 to 1.0 pop_1; 0.0 to 0.21 pop_2
## MO (6) MI (2)

# MW_3_2 = MW-4X close to MW-8X West
## 0.52 to 0.69 pop_1; 0.31 to 0.48 pop_2
## MO (5) AK (4)

# MW_3_3 = MW-4X distant to MW-8X West
## 0.17 to 0.35 pop_1; 0.64 to 0.83 pop_2
## MO (18)

# MW_3_4 = MW-4X distant to MW-8X West
## 0.19 to 0.31 pop_1; 0.69 to 0.81 pop_2
## MO (2) NE (2) MN (2) WI (3) NY (1)

# MW_4_1 = MW-4X-GreatLakes
## pure pop_2
## MI (23) IL (16) WI (6) IN (3) 

# MW_4_2
## 0.94 pop_2
## IN (1)

# MW_5_1 = MW-4X-East1
## pure pop_2 
## OH (10)

# MW_5_2 = MW-4X-East2
## pure pop_2
## OH (6)

# MW_5_3 = MW-4X-UpperMW
## 0.85 to 1.0 pop_2; 0.0 to 0.15 pop_1
## MN (10) WI (2) IA (1) NY (1)

# NOTE: MW_4_1, MW_5_1, and MW_5_2 all have the same patterns from the
## ADMIXTURE results, but show some distinct clustering in PCA, 
## though only at the higher PCs (4 and 5)
# I could be conviced to group all of MW_4_X and MW_5_X together if we want to

###

# NOTES ON SPLITTING UP GULF
# ADMIX Gulf_main gives K=2; Gulf_tot gives K=3, maybe K=4

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
gulf_res_tab[grp_g2_1, grps := 'GULF_2_1']
gulf_res_tab[grp_g2_2, grps := 'GULF_2_2']
gulf_res_tab[grp_g2_3, grps := 'GULF_2_3']
gulf_res_tab[grp_g3_1, grps := 'GULF_3_1']
gulf_res_tab[grp_g3_2, grps := 'GULF_3_2']
gulf_res_tab[grp_g3_3, grps := 'GULF_3_3']
gulf_res_tab[grp_g3_4, grps := 'GULF_3_4']
gulf_res_tab[grp_g3_5, grps := 'GULF_3_5']
gulf_res_tab[grp_g3_6, grps := 'GULF_3_6']
gulf_res_tab[grp_AM_b, grps := 'GULF_AM_b']
gulf_res_tab[grp_AM_1, grps := 'GULF_AM_1']
gulf_res_tab[grp_AM_2, grps := 'GULF_AM_2']
gulf_res_tab[grp_AM_3, grps := 'GULF_AM_3']
gulf_res_tab[grp_AM_4, grps := 'GULF_AM_4']

which(gulf_res_tab$gulf_tot_PC1 > -200 & gulf_res_tab$grps == 'NA')
# all Gulf samples accounted for

gulf_res_tab$grps <- as.factor(gulf_res_tab$grps)

levels(gulf_res_tab$grps) <- sort(levels(gulf_res_tab$grps))

gg_gulf_tmp <- ggplot(gulf_res_tab, aes(x = gulf_tot_PC1, y = gulf_tot_PC2, 
                                        color = grps, 
                                      shape = ploidy)) +
  geom_point()

gg_gulf_tmp

# Look at State distribution
for(pop_name in sort(unique(gulf_res_tab$grps))){
  print(pop_name)
  print(table(samp_meta$STATE[
    samp_meta$VCF_NAME %in% gulf_res_tab$samp_name[which(
      gulf_res_tab$grps == pop_name)]]))
}

############
# Gulf main K=2 results
table(gulf_res_tab$grps[which(gulf_res_tab$gulf_main_k2_pop_1 > 0.5)])
# GULF_2_1, GULF_2_2 - top-right corner of PC1vsPC2 plot
table(gulf_res_tab$grps[which(gulf_res_tab$gulf_main_k2_pop_2 > 0.5)])
# GULF_3_3, GULF_3_4, GULF_3_7 - main bottom-right cluster

for(pop_name in sort(unique(gulf_res_tab$grps))){
  print(pop_name)
  print(summary(gulf_res_tab$gulf_main_k2_pop_1[
    gulf_res_tab$grps == pop_name]))
}

########
# Gulf tot K=3 results
table(gulf_res_tab$grps[which(gulf_res_tab$gulf_tot_k3_pop_1 > 0.5)])
# GULF_2_3 (MS, FL, LA), GULF_3_6 (FL, LA)
table(gulf_res_tab$grps[which(gulf_res_tab$gulf_tot_k3_pop_2 > 0.5)])
# GULF_3_1 (AR, LA), GULF_3_2 (LA, TX), GULF_3_3 (TX) GULF_3_4 (TX) 
# GULF_3_5 (FL), GULF_3_6 (FL, LA), GULF_AM_1 (TX), GULF_AM_2 (FL) 
# GULF_AM_3 (TX), GULF_AM_b (TX)
table(gulf_res_tab$grps[which(gulf_res_tab$gulf_tot_k3_pop_3 > 0.5)])
# GULF_2_1 (MX) GULF_2_2 (TX) GULF_AM_3 (TX, SC), GULF_AM_b (TX)

## GOOD REASON TO SPLIT UP GULF_3_7 is that is split between pop2 and pop3?

for(pop_name in sort(unique(gulf_res_tab$grps))){
  print(pop_name)
  print(summary(gulf_res_tab$gulf_tot_k3_pop_3[
    gulf_res_tab$grps == pop_name]))
}

###### GULF Total K=3 Explanation
# Pop 1 = Gulf Coast 1
# Pop 2 = Gulf Coast 2
# Pop 3 = TX/MX

###### GULF subgroup explanations
# GULF_2_1 - pure TX/MX (pop_3)
## MX (11)

# GULF_2_2 - pure TX/MX (pop_3)
## TX(8)

# GULF_2_3 - Pure GC 1 (pop_1)
## MS(20) FL (8) LA (7)

# GULF_3_1 - Pure GC 2 (pop_2)
## AR (12) LA (3)

# GULF_3_2 - Mainly GC 2 (pop_2) with a little GC 1
## 0.82 to 0.96 pop_2; 0.04 to 0.18 pop_1
## LA(6) TX (1)

# GULF_3_3 - Mainly GC 2 (pop_2) with a little TX/MX
## 0.78 to 0.92 pop_2; 0.08 to 0.22 pop_3
## TX (11)

# GULF_3_4 - Mainly GC 2 (pop_2) with a little TX/MX
## 0.74 to 0.81 pop_2; 0.19 to 0.26 pop_3
## TX (4)

# GULF_3_5 - Mainly GC 2 (pop_2) with some GC 1
## 0.64 to 0.71 pop_2; 0.29 to 0.36 pop_1
## FL (3)

# GULF_3_6 - Majority Gulf Coast 1 (pop_1) with some GC 2 and less TX/MX
## 0.49 to 0.96 pop_1; 0.04 to 0.51 pop_2; 0.0 to 0.09 pop_3
## LA (9) FL (4)

# GULF_AM_1 - Majority GC 2 (pop_2) with TX/MX and a bit GC 1
## 0.72 to 0.86 pop_2; 0.07 to 0.23 pop_3; 0.02 to 0.06 pop_1
## TX (4)

# GULF_AM_2 - Majority GC 2 (pop 2) with GC 1 and/or TX/MX
## 0.71 to 0.72 pop_2; 0.13 to 0.20 pop_3, 0.07 to 0.17 pop_1
## FL (3)

# GULF_AM_3 - Majority GC 2 (pop_2) with some TX/MX
## 0.63 to 0.89 pop_2, 0.11 to 0.37 pop_3
## TX (18)

# GULF_AM_3 - TX/MX (pop_3) with some GC 2 (pop_2)
## 0.63 to 0.80 pop_3; 0.20 to 0.37 pop_2
## TX (18) SC (2)

# GULF_AM_b - roughly half TX/MX (pop_3) and GC 2 (pop_2) 
## 0.46 to 0.59 pop_3; 0.41 to 0.54 pop_2
## TX (5)


#########
atl_tot_res$eig[1:10]/sum(atl_tot_res$eig)

atl_meta_inds <- c()
for(i in rownames(atl_tot_res$scores)){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  atl_meta_inds <- c(atl_meta_inds, tmp_ind)
}

res_atl_tab <- data.table(samp_names = rownames(atl_tot_res$scores), 
                           ploidy = samp_meta$NQUIRE_PLOIDY[atl_meta_inds],
                           PC1 = atl_tot_res$scores[,1], 
                           PC2 = atl_tot_res$scores[,2],
                           PC3 = atl_tot_res$scores[,3],
                           PC4 = atl_tot_res$scores[,4],
                           PC5 = atl_tot_res$scores[,5],
                           PC6 = atl_tot_res$scores[,6])

gg_atl_tmp <- ggplot(res_atl_tab, aes(x = PC1, y = PC2)) +
  geom_point(color = res_atl_tab$ploidy)

gg_atl_tmp

# NOTES FOR SPLITTING ATL
# ADMIXTURE ATL_main is K=3; ATL_tot is K=3 or K=4

grp_g2_1 <- which(gulf_res_tab$gulf_tot_PC2 > 15)


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
atl_res_tab[grp_a2_1, grps := 'ATL_2_1']
atl_res_tab[grp_a2_2, grps := 'ATL_2_2']
atl_res_tab[grp_a2_3, grps := 'ATL_2_3']
atl_res_tab[grp_a2_4, grps := 'ATL_2_4']
atl_res_tab[grp_a2_5, grps := 'ATL_2_5']
atl_res_tab[grp_a2_6, grps := 'ATL_2_6']
atl_res_tab[grp_a2_7, grps := 'ATL_2_7']
atl_res_tab[grp_a3_1, grps := 'ATL_3_1']
atl_res_tab[grp_a3_2, grps := 'ATL_3_2']
atl_res_tab[grp_a3_3, grps := 'ATL_3_3']
atl_res_tab[grp_a3_4, grps := 'ATL_3_4']
atl_res_tab[grp_a4_1, grps := 'ATL_4_1']
atl_res_tab[grp_a4_2, grps := 'ATL_4_2']
atl_res_tab[grp_a4_3, grps := 'ATL_4_3']

which(atl_res_tab$atl_tot_PC1 > -200 & atl_res_tab$grps == 'NA')
# 0

atl_res_tab$grps <- as.factor(atl_res_tab$grps)
levels(atl_res_tab$grps) <- sort(levels(atl_res_tab$grps))

gg_atl_tmp <- ggplot(atl_res_tab, aes(x = atl_tot_PC1, y = atl_tot_PC2, 
                                      color = grps, 
                                      shape = ploidy)) +
  geom_point()

gg_atl_tmp

# Look at State distribution
for(pop_name in sort(unique(atl_res_tab$grps))){
  print(pop_name)
  print(table(samp_meta$STATE[
    samp_meta$VCF_NAME %in% atl_res_tab$samp_name[which(
      atl_res_tab$grps == pop_name)]]))
}

############
# Atl main K=3 results
table(atl_res_tab$grps[which(atl_res_tab$atl_main_k3_pop_1 > 0.5)])
# ATL_2_1, ATL_4_1, ATL_4_2, ATL_4_3 (NE)
table(atl_res_tab$grps[which(atl_res_tab$atl_main_k3_pop_2 > 0.5)])
# ATL_2_3, ATL_2_4, ATL_2_5, ATL_3_1 (Mid-Atlantic and Carolinas)
table(atl_res_tab$grps[which(atl_res_tab$atl_main_k3_pop_3 > 0.5)])
# ATL_3_2, ATL_3_3 (NY + Mid-Atlantic)

for(pop_name in sort(unique(atl_res_tab$grps))){
  print(pop_name)
  print(summary(atl_res_tab$atl_main_k3_pop_3[
    atl_res_tab$grps == pop_name]))
}

### ATL Tot K=3 results
table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k3_pop_1 > 0.5)])
# ATL_2_1, ATL_4_1, ATL_4_2, ATL_4_3 (NE + NY)

table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k3_pop_2 > 0.5)])
# ATL_2_2, ATL_2_3, ATL_2_4, ATL_2_5, ATL_3_1, ATL_3_4 (FL, NC/SC, MidAtlantic)

table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k3_pop_3 > 0.5)])
# ATL_2_5 ATL_3_1, ATL_3_2, ATL_3_3 (MidAtlantic + NY)

for(pop_name in sort(unique(atl_res_tab$grps))){
  print(pop_name)
  print(summary(atl_res_tab$atl_tot_k3_pop_2[
    atl_res_tab$grps == pop_name]))
}

### ATL Tot K=4 results
table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k4_pop_1 > 0.5)])
# ATL_2_5, ATL_3_2, ATL_3_3 (NY, MidAtlantic, NC)
table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k4_pop_2 > 0.5)])
# ATL_2_1, ATL_4_1, ATL_4_2, ATL_4_3 (NE + NY)
table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k4_pop_3 > 0.5)])
# ATL_2_2, ATL_2_3, ATL_2_4, ATL_2_5, ATL_3_4 (FL and SE coast)
table(atl_res_tab$grps[which(atl_res_tab$atl_tot_k4_pop_4 > 0.5)])
# ATL_2_6, ATL_3_1 (MidAtlantic)

for(pop_name in sort(unique(atl_res_tab$grps))){
  print(pop_name)
  print(summary(atl_res_tab$atl_tot_k4_pop_4[
    atl_res_tab$grps == pop_name]))
}

table(samp_meta$STATE[
  samp_meta$VCF_NAME %in% 
    atl_res_tab$samp_name[which(atl_res_tab$atl_tot_k4_pop_4 > 0.5)]])

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
# ATl_2_1 - New England and NY
## 0.75 to 1.0 pop_2
## CT (8) ME (3) MA (3) NH (3) NY (20) RI (20)

# ATL_2_2 - FL-1
## pure pop_3
## FL (9)

# ATL_2_3 - SE Atlantic Coast
## Pure pop_3
## FL (1) GA (3) SC (2)

# ATL_2_4 - Carolina Coast
## 0.63 to 0.97 pop_3; 0.03 to 0.24 pop_4
## NC (11) SC (17)

# ATL_2_5 - NC
## 0.28 to 0.61 pop_1, 0.26 to 0.62 pop_3; 0.08 to 0.21 pop_4
## NC (25)

# ATL_2_6 - VA
## 0.53 to 0.68 pop_4; 0.32 to 0.46 pop_1, 
## VA (8)

# ATL_2_7 - Southern MidAtlantic
## 0.37 to 0.76 pop_4; 0.16 to 0.53 pop_1; 0.05 to 0.29 pop_3
## MD (1) NC (3) VA (2)

# ATL_3_1 - MidAtlantic
## 0.82 to 1.0 pop_4
## MD (23) VA (7)

# ATL_3_2 - NY-1
## 0.58 to 0.70 pop_1, 0.30 to 0.42 pop_2
## NY (5)

# ATL_3_3 - Northern MidAtlantic
## 0.73 to 1.0 pop_1; 0.0 to 0.27 pop_4
## DE (13) MD (2) NJ (24) NY (27) VA (3)

# ATL_3_4 - FL-2
## pure pop_3
## FL (11)

# ATL_4_1 - MA
## 0.76 to 0.86 pop_2
## MA (6)

# ATL_4_2 - NY-2
## 0.64 to 0.71 pop_2 
## NY (8)

# ATL_4_3 - NY-3
## 0.80 to 0.99 pop_2
## NY (14)

###

misfit_res$eig[1:10]/sum(misfit_res$eig)

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

## NOTES FOR SPLITTING MISFITS

misfit_tab[which(misfit_tab$PC1 > 20)]
grp_mf1_1 <- which(misfit_tab$PC1 > 30)
grp_mf1_2 <- which(misfit_tab$PC1 > 20 & misfit_tab$PC1 < 30)
grp_mf1_3 <- which(misfit_tab$PC1 > -2.5 & misfit_tab$PC1 < 1)

grp_mf2_1 <- which(misfit_tab$PC2 < -10)
grp_mf2_2 <- which(misfit_tab$PC2 > 10 & misfit_tab$PC2 < 15)
grp_mf2_3 <- which(misfit_tab$PC1 > 10 & misfit_tab$PC1 < 15 &
                     misfit_tab$PC2 > 0)
# grp_mf3_1 <- which(misfit_tab$PC3 > 20)
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
misfit_tab[grp_mf1_1, grps := 'MF_1_1']
misfit_tab[grp_mf1_2, grps := 'MF_1_2']
misfit_tab[grp_mf1_3, grps := 'MF_1_3']
misfit_tab[grp_mf2_1, grps := 'MF_2_1']
misfit_tab[grp_mf2_2, grps := 'MF_2_2']
misfit_tab[grp_mf2_3, grps := 'MF_2_3']
#misfit_tab[grp_mf3_1, grps := 'MF_3_1']
misfit_tab[grp_mf3_2, grps := 'MF_3_2']
misfit_tab[grp_mf3_3, grps := 'MF_3_3']
misfit_tab[grp_mf_rest, grps := 'MF_REST']

misfit_tab$grps <- as.factor(misfit_tab$grps)
levels(misfit_tab$grps) <- sort(levels(misfit_tab$grps))

gg_misfit_tmp <- ggplot(misfit_tab, aes(x = PC1, y = PC2, 
                                        color = grps,
                                        shape = ploidy)) +
  geom_point()

gg_misfit_tmp

for(pop_name in sort(unique(misfit_tab$grps))){
  print(pop_name)
  print(table(samp_meta$STATE[
    samp_meta$VCF_NAME %in% misfit_tab$samp_name[which(
      misfit_tab$grps == pop_name)]]))
}

for(pop_name in sort(unique(misfit_tab$grps))){
  print(pop_name)
  print(table(misfit_tab$ploidy[which(misfit_tab$grps == pop_name)]))
}

# MF_1_1 / MF_01 - FL-4X
## 4X
## FL (2)

# MF_1_2 / MF_02 - NC-4X
## 4X
## NC (2) 

# MF_1_3 / MF_03 - GA-8X-1
## 8X
## GA (2) NC (1)

# MF_2_1 / MF_04 - NY/NJ
## 4X (6) 8X (1)
## NJ (3) NY (4)

# MF_2_2 / MF_05 - GA-8X-2
## 8X (3)
## GA (3)

# MF_2_3 / MF_06 - FL-8X
## 8X (1)
## FL (1)

# MF_3_2 / MF_07 - GA-8X-3
## 8X (7)
## GA (5) IL (1) TN (1)

# MF_3_3 - MS-8X
## 8X (11) 6X (1) 
## AR (1) FL (1) MS (10)