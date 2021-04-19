### I had to move to a v2.0 of this document to keep it from getting too
###   unweildy

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###
#mw_tot_pca_res_file <- '/Users/grabowsk/data/sg_8X_analysis/GW.100kSNPs.tetrasomic.CDS.natv2filt.MW_tot.genlight.PCAresults.rds'
#mw_tot_res <- readRDS(mw_tot_pca_res_file)

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

#gulf_tot_pca_res_file <- '/Users/grabowsk/data/sg_8X_analysis/GW.100kSNPs.tetrasomic.CDS.natv2filt.GULF_tot.genlight.PCAresults.rds'
#gulf_tot_res <- readRDS(gulf_tot_pca_res_file)

#atl_tot_pca_res_file <- '/Users/grabowsk/data/sg_8X_analysis/GW.100kSNPs.tetrasomic.CDS.natv2filt.ATL_tot.genlight.PCAresults.rds'
#atl_tot_res <- readRDS(atl_tot_pca_res_file)

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
### Using PCA results before compiled all results ###
mw_tot_res$eig[1:10]/sum(mw_tot_res$eig)

mw_meta_inds <- c()
for(i in rownames(mw_tot_res$scores)){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  mw_meta_inds <- c(mw_meta_inds, tmp_ind)
}

res_mw_tab <- data.table(samp_names = rownames(mw_tot_res$scores), 
                      ploidy = samp_meta$NQUIRE_PLOIDY[mw_meta_inds],
                      PC1 = mw_tot_res$scores[,1], 
                      PC2 = mw_tot_res$scores[,2],
                      PC3 = mw_tot_res$scores[,3],
                      PC4 = mw_tot_res$scores[,4],
                      PC5 = mw_tot_res$scores[,5],
                      PC6 = mw_tot_res$scores[,6])

gg_tmp <- ggplot(res_mw_tab, aes(x=PC1, y=PC5)) + 
  geom_point(color = res_mw_tab$ploidy)

gg_tmp

# Notes on splitting up samples
# First grouping is defnitely to split into 3 - maybe use ADMIX MW_tot K=3
### results - I think this will break up the 8X into 2
# For Splitting up 4X:
grp_mw2_1 <- which(res_mw_tab$PC1 < -5 & res_mw_tab$PC2 < -7.5)
grp_mw2_2 <- which(res_mw_tab$PC1 > -5 & res_mw_tab$PC1 < -2.5 &
                     res_mw_tab$PC2 > 0)
grp_mw2_3 <- which(res_mw_tab$PC1 < -5 & res_mw_tab$PC2 > -7)
grp_mw3_1 <- which(res_mw_tab$PC1 > -2.5 & res_mw_tab$PC1 < 0)
grp_mw3_2 <- which(res_mw_tab$PC1 > 0 & res_mw_tab$PC1 < 5
                & res_mw_tab$PC3 > 0)
grp_mw3_3 <- which(res_mw_tab$PC1 > 5 & res_mw_tab$PC1 < 9 
                & res_mw_tab$PC3 > 0)
grp_mw3_4 <- which(res_mw_tab$PC1 > 2.5 & res_mw_tab$PC1 < 8
                & res_mw_tab$PC3 < 0)
grp_mw4_1 <- which(res_mw_tab$PC1 > 10 & res_mw_tab$PC4 > 0)
grp_mw5_1 <- which(res_mw_tab$PC5 < -10)

tot_mw_grps <- c(grp_mw2_1, grp_mw2_2, grp_mw2_3, grp_mw3_1, grp_mw3_2, 
                 grp_mw3_3, grp_mw3_4, grp_mw4_1, grp_mw5_1)
sum(duplicated(tot_mw_grps))
#[1] 0

res_mw_tab[, grps := 'MW_rest']
res_mw_tab[grp_mw2_1, grps := 'MW_2_1']
res_mw_tab[grp_mw2_2, grps := 'MW_2_2']
res_mw_tab[grp_mw2_3, grps := 'MW_2_3']
res_mw_tab[grp_mw3_1, grps := 'MW_3_1']
res_mw_tab[grp_mw3_2, grps := 'MW_3_2']
res_mw_tab[grp_mw3_3, grps := 'MW_3_3']
res_mw_tab[grp_mw3_4, grps := 'MW_3_4']
res_mw_tab[grp_mw4_1, grps := 'MW_4_1']
res_mw_tab[grp_mw5_1, grps := 'MW_5_1']

res_mw_tab$grps <- as.factor(res_mw_tab$grps)

gg_mw_tmp <- ggplot(res_mw_tab, aes(x = PC1, y = PC2, color = grps, 
                                        shape = ploidy)) +
  geom_point()

gg_mw_tmp

###
###
###
gulf_tot_res$eig[1:10]/sum(gulf_tot_res$eig)

gulf_meta_inds <- c()
for(i in rownames(gulf_tot_res$scores)){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  gulf_meta_inds <- c(gulf_meta_inds, tmp_ind)
}

res_gulf_tab <- data.table(samp_names = rownames(gulf_tot_res$scores), 
                         ploidy = samp_meta$NQUIRE_PLOIDY[gulf_meta_inds],
                         PC1 = gulf_tot_res$scores[,1], 
                         PC2 = gulf_tot_res$scores[,2],
                         PC3 = gulf_tot_res$scores[,3],
                         PC4 = gulf_tot_res$scores[,4],
                         PC5 = gulf_tot_res$scores[,5],
                         PC6 = gulf_tot_res$scores[,6])

gg_gulf_tmp <- ggplot(res_gulf_tab, aes(x = PC1, y = PC2)) +
  geom_point(color = res_gulf_tab$ploidy)

gg_gulf_tmp

# NOTES ON SPLITTING UP GULF
# ADMIX Gulf_main gives K=2; Gulf_tot gives K=3, maybe K=4
grp_g2_1 <- which(res_gulf_tab$PC2 > 15)
grp_g2_2 <- which(res_gulf_tab$PC1 > 16 & res_gulf_tab$PC2 < 15)
grp_g2_3 <- which(res_gulf_tab$PC1 < -17)
grp_g3_1 <- which(res_gulf_tab$PC3 > 10)
grp_g3_2 <- which(res_gulf_tab$PC3 > 0 & res_gulf_tab$PC3 < 10 & 
                    res_gulf_tab$PC1 > -5 & res_gulf_tab$PC1 < 0)
grp_g3_3 <- which(res_gulf_tab$PC3 < -10 & res_gulf_tab$PC3 > -16 &
                    res_gulf_tab$PC1 > 5 & res_gulf_tab$PC1 < 10)
grp_g3_4 <- which(res_gulf_tab$PC1 > 8 & res_gulf_tab$PC1 < 12 &
                    res_gulf_tab$PC3 < -15)
grp_g3_5 <- which(res_gulf_tab$PC1 > -7 & res_gulf_tab$PC1 < 0 &
                    res_gulf_tab$PC3 < -20)
grp_g3_6 <- which(res_gulf_tab$PC1 > -17.5 & res_gulf_tab$PC1 < -6)
grp_g3_7 <- which(res_gulf_tab$PC1 > 0 & res_gulf_tab$PC16 &
                    res_gulf_tab$PC3 > -10 & res_gulf_tab$PC3 < 5)

gulf_allgrps <- c(grp_g2_1, grp_g2_2, grp_g2_3, grp_g3_1, grp_g3_2, grp_g3_3,
                  grp_g3_4, grp_g3_5, grp_g3_6, grp_g3_7)
sum(duplicated(gulf_allgrps))
#[1] 0

res_gulf_tab[, grps := 'GULF_rest']
res_gulf_tab[grp_g2_1, grps := 'GULF_2_1']
res_gulf_tab[grp_g2_2, grps := 'GULF_2_2']
res_gulf_tab[grp_g2_3, grps := 'GULF_2_3']
res_gulf_tab[grp_g3_1, grps := 'GULF_3_1']
res_gulf_tab[grp_g3_2, grps := 'GULF_3_2']
res_gulf_tab[grp_g3_3, grps := 'GULF_3_3']
res_gulf_tab[grp_g3_4, grps := 'GULF_3_4']
res_gulf_tab[grp_g3_5, grps := 'GULF_3_5']
res_gulf_tab[grp_g3_6, grps := 'GULF_3_6']
res_gulf_tab[grp_g3_7, grps := 'GULF_3_7']

res_gulf_tab$grps <- as.factor(res_gulf_tab$grps)

gg_gulf_tmp <- ggplot(res_gulf_tab, aes(x = PC1, y = PC2, color = grps, 
                                      shape = ploidy)) +
  geom_point()

gg_gulf_tmp

######

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
grp_a2_2 <- which(res_atl_tab$PC2 > 15)
grp_a2_3 <- which(res_atl_tab$PC2 > 9 & res_atl_tab$PC2 < 11)
grp_a2_4 <- which(res_atl_tab$PC1 > 5 & res_atl_tab$PC2 > 2.5 &
                    res_atl_tab$PC2 < 9)
grp_a3_1 <- which(res_atl_tab$PC3 < -7)
grp_a3_2 <- which(res_atl_tab$PC1 > -7.5 & res_atl_tab$PC1 < -2.5 &
                    res_atl_tab$PC3 > 0 & res_atl_tab$PC3 < 6)
grp_a3_3 <- which((res_atl_tab$PC3 > 0 & res_atl_tab$PC2 < -5) |
                    res_atl_tab$PC3 > 10)
grp_a3_4 <- which(res_atl_tab$PC2 > 11 & res_atl_tab$PC2 < 15)
grp_a4_1 <- which(res_atl_tab$PC4 < -5 & res_atl_tab$PC1 > -15 &
                    res_atl_tab$PC1 < -5)
grp_a4_2 <- which(res_atl_tab$PC1 > -11 & res_atl_tab$PC1 < -8
                  & res_atl_tab$PC4 > 0 & res_atl_tab$PC4 < 2.5)
grp_a4_3 <- which(res_atl_tab$PC1 < -10 & res_atl_tab$PC1 > -13 & 
                    res_atl_tab$PC4 > 2.5 | 
                    (res_atl_tab$PC1 < -11 & res_atl_tab$PC1 > -13 & 
                       res_atl_tab$PC4 > 1.5))
grp_a2_1 <- setdiff(which(res_atl_tab$PC1 < -7.5), 
                    c(grp_a4_1, grp_a4_2, grp_a4_3))

grp_a_all <- c(grp_a2_1, grp_a2_2, grp_a2_3, grp_a2_4,grp_a3_1, grp_a3_2, 
               grp_a3_3, grp_a3_4, grp_a4_1, grp_a4_2, grp_a4_3)
grp_a_rest <- setdiff(seq(nrow(res_atl_tab)), grp_a_all)

res_atl_tab[, grps := 'ATL_rest']
res_atl_tab[grp_a2_1, grps := 'ATL_2_1']
res_atl_tab[grp_a2_2, grps := 'ATL_2_2']
res_atl_tab[grp_a2_3, grps := 'ATL_2_3']
res_atl_tab[grp_a2_4, grps := 'ATL_2_4']
res_atl_tab[grp_a3_1, grps := 'ATL_3_1']
res_atl_tab[grp_a3_2, grps := 'ATL_3_2']
res_atl_tab[grp_a3_3, grps := 'ATL_3_3']
res_atl_tab[grp_a3_4, grps := 'ATL_3_4']
res_atl_tab[grp_a4_1, grps := 'ATL_4_1']
res_atl_tab[grp_a4_2, grps := 'ATL_4_2']
res_atl_tab[grp_a4_3, grps := 'ATL_4_3']

res_atl_tab$grps <- as.factor(res_atl_tab$grps)

gg_atl_tmp <- ggplot(res_atl_tab, aes(x = PC1, y = PC4, color = grps, 
                                      shape = ploidy)) +
  geom_point()

gg_atl_tmp

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
grp_mf1_2 <- which(misfit_tab$PC1 > 30 & misfit_tab$PC1 < 30)
grp_mf1_3 <- which(misfit_tab$PC1 > -2.5 & misfit_tab$PC1 < 1)

grp_mf2_1 <- which(misfit_tab$PC2 < -10)
grp_mf2_2 <- which(misfit_tab$PC2 > 10 & misfit_tab$PC2 < 15)
grp_mf2_3 <- which(misfit_tab$PC1 > 10 & misfit_tab$PC1 < 15 &
                     misfit_tab$PC2 > 0)
grp_mf3_1 <- which(misfit_tab$PC3 > 30)
grp_mf6_1 <- which(misfit_tab$PC1 < -7 & misfit_tab$PC1 < -2 &
                     misfit_tab$PC6 < -2)
grp_tmp_all <- c(grp_mf1_1, grp_mf1_2, grp_mf1_3, grp_mf2_1, grp_mf2_2, 
                 grp_mf2_3, grp_mf3_1, grp_mf6_1)
grp_mf_rest <- setdiff(seq(nrow(misfit_tab)), grp_tmp_all)

