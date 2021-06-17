# Characterizing and plots of subgroups and full-data info

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

### INPUT DATA ###

samp_meta_file <- paste('/Users/grabowsk/Analysis/sg_8X/metadata_management/',
  'meta_files/', 'sg_8X_metadata_v3.1.csv', sep = '')
samp_meta <- fread(samp_meta_file)

res_tab_file <- paste('/Users/grabowsk/data/sg_8X_analysis/', 
  'natv2filt_res_tab_v4.0.txt', sep = '')
res_tab <- fread(res_tab_file)

### SET OUTPUT ###

atl_plot_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                      'natv2filt_PC1v2_ATL_colors.pdf', sep = '')
gulf_plot_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                      'natv2filt_PC1v2_GULF_colors.pdf', sep = '')
mw_plot_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                       'natv2filt_PC1v2_MW_colors.pdf', sep = '')
mf_plot_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                     'natv2filt_PC1v2_MISFIT_colors.pdf', sep = '')

#####################

table(res_tab$sub_grp[res_tab$full_admix_k3_ATL > 0.5])
# ATL subgroups and MF_01, MF_02, MF_04, MF_06

table(res_tab$sub_grp[res_tab$full_admix_k3_GULF > 0.5])
# Only GULF subgroups

table(res_tab$sub_grp[res_tab$full_admix_k3_MW > 0.5])
# MW subgroups and MF_04, MF_07< MF_08

# Note: MF_04 found in both ATL > 0.5 and MW > 0.5

################
# ATL subgroups looking at Full K=3
all_groups <- sort(unique(res_tab$sub_grp))

atl_groups <- all_groups[grep('ATL', all_groups)]

for(pop_name in atl_groups){
  print(paste('#####', pop_name, '######'))
  print(paste(pop_name, 'K=3 ATL'))
  print(summary(res_tab$full_admix_k3_ATL[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 GULF'))
  print(summary(res_tab$full_admix_k3_GULF[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 MW'))
  print(summary(res_tab$full_admix_k3_MW[res_tab$sub_grp == pop_name]))
}

### ATL K=4 Explanation
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
## pure Full_K3_ATL
## 0.75 to 1.0 ATL pop_2
## CT (8) ME (3) MA (3) NH (3) NY (20) RI (20)

# ATL_02 - FL-1
## 0.9 to 0.94 Full_K3_ATL; 0.06 to 0.1 Full_K3_GULF
## pure ATL pop_3
## FL (9)

# ATL_03 - SE Atlantic Coast
## 0.74 to 0.97 Full_K3_ATL; 0.03 to 0.26 Full_K3_GULF
## Pure ATL pop_3
## FL (1) GA (3) SC (2)

# ATL_04 - Carolina Coast
## Almost pure Full_K3_ATL; trace amounts of FULL_K3_GULF 
## 0.63 to 0.97 ATL pop_3; 0.03 to 0.24 ATL pop_4
## NC (11) SC (17)

# ATL_05 - NC
## pure Full_K3_ATL
## 0.28 to 0.61 ATL pop_1, 0.26 to 0.62 ATL pop_3; 0.08 to 0.21 ATL pop_4
## NC (25)

# ATL_06 - VA
## pure Full_K3_ATL
## 0.53 to 0.68 ATL pop_4; 0.32 to 0.46 ATL pop_1, 
## VA (8)

# ATL_07 - Southern MidAtlantic
## 0.68 to 0.98 Full_K3_ATL; 0.02 to 0.32 Full_K3_GULF
## 0.37 to 0.76 ATL pop_4; 0.16 to 0.53 ATL pop_1; 0.05 to 0.29 ATL pop_3
## MD (1) NC (3) VA (2)

# ATL_08 - MidAtlantic
## 0.8 to 1.0 Full_K3_ATL; 0.0 to 0.20 Full_K3_GULF
## 0.82 to 1.0 ATL pop_4
## MD (23) VA (7)

# ATL_09 - NY-1
## pure Full_K3_ATL
## 0.58 to 0.70 ATL pop_1, 0.30 to 0.42 ATL pop_2
## NY (5)

# ATL_10 - Northern MidAtlantic
## almost pure Full_K3_ATL with traces amounts of both ..GULF and ..MW
## 0.73 to 1.0 ATL pop_1; 0.0 to 0.27 ATL pop_4
## DE (13) MD (2) NJ (24) NY (27) VA (3)

# ATL_11 - FL-2
## 0.7 to 0.91 Full_K3_ATL; 0.09 to 0.3 Full_K3_GULF, trace amounds of ...MW
## pure ATL pop_3
## FL (11)

# ATL_12 - MA
## pure Full_K3_ATL
## 0.76 to 0.86 ATL pop_2
## MA (6)

# ATL_13 - NY-2
## pure Full_K3_ATL
## 0.64 to 0.71 ATL pop_2 
## NY (8)

# ATL_14 - NY-3
## pure Full_K3_ATL
## 0.80 to 0.99 ATL pop_2
## NY (14)

########

res_tab[, ATL_GROUPS := 'OTHER']
atl_inds <- grep('ATL', res_tab$sub_grp)
res_tab[atl_inds, ATL_GROUPS := res_tab$sub_grp[atl_inds]]
res_tab$ATL_GROUPS <- as.factor(res_tab$ATL_GROUPS)
levels(res_tab$ATL_GROUPS) <- sort(unique(res_tab$ATL_GROUPS))

atl_palette <- rainbow(15)
atl_palette[15] <- 'black'
gg_atl <- ggplot(res_tab, aes(x = full_pc01_stand, y = full_pc02_stand, 
                              color = ATL_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('All Samples, Atlantic sub group colors') +
  scale_color_manual(values = atl_palette)

gg_atl

gg_atl_sub <- ggplot(res_tab, aes(x = ATL_pc01_stand, y = ATL_pc02_stand,
                                  color = ATL_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('Atlantic Samples') +
  scale_color_manual(values = atl_palette)

gg_atl_sub

pdf(atl_plot_out, height = 6, width = 13)
gg_atl + gg_atl_sub
dev.off()

#######
###
# GULF subgroups looking at Full K=3

gulf_groups <- all_groups[grep('GULF', all_groups)]

for(pop_name in gulf_groups){
  print(paste('#####', pop_name, '######'))
  print(paste(pop_name, 'K=3 GULF'))
  print(summary(res_tab$full_admix_k3_GULF[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 ATL'))
  print(summary(res_tab$full_admix_k3_ATL[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 MW'))
  print(summary(res_tab$full_admix_k3_MW[res_tab$sub_grp == pop_name]))
}

###### GULF Total K=3 Explanation
# Pop 1 = Gulf Coast 1
# Pop 2 = Gulf Coast 2 (perhaps some shared MW ancestry?)
# Pop 3 = TX/MX

###### GULF subgroup explanations
# GULF_01 - pure TX/MX
## pure Full_K3_GULF
## pure TX/MX (GULF pop_3)
## MX (11)

# GULF_02 - pure TX/MX
## pure Full_K3_GULF
## pure TX/MX (GULF pop_3)
## TX(8)

# GULF_03 - Pure GC 1
## 0.93 to 1.0 Full_K3_GULF; 0.0 to 0.07 ...ATL; trace ...MW
## Pure GC 1 (GULF pop_1)
## MS(20) FL (8) LA (7)

# GULF_04 - Pure GC 2
## 0.80 to 0.92 Full_K3_GULF; 0.08 to 0.20 Full_K3_MW
## Pure GC 2 (GULF pop_2)
## AR (12) LA (3)

# GULF_05 - Mainly GC 2 with a little GC 1
## 0.76 to 0.81 Full_K3_GULF; 0.19 to 0.24 Full_K3_MW
## 0.82 to 0.96 GULF pop_2; 0.04 to 0.18 GULF pop_1
## LA(6) TX (1)

# GULF_06 - Mainly GC 2 with a little TX/MX
## Almost pure Full_K3_GULF; trace amounts of Full_K3_MW
## 0.78 to 0.92 GULF pop_2; 0.08 to 0.22 GULF pop_3
## TX (11)

# GULF_07 - Mainly GC 2 (pop_2) with a little TX/MX
## Pure Full_K3_GULF
## 0.74 to 0.81 GULF pop_2; 0.19 to 0.26 GULF pop_3
## TX (4)

# GULF_08 - Mainly GC 2 with some GC 1
## 0.93 to 0.95 Full_K3_GULF; 0.05 to 0.07 Full_K3_ATL
## 0.64 to 0.71 GULF pop_2; 0.29 to 0.36 GULF pop_1
## FL (3)

# GULF_09 - Majority GC 1 with some GC 2 and less TX/MX
## 0.75 to 0.96 Full_K3_GULF; 0.0 to 0.25 ..ATL; 0.0 to 0.06 ..MW
## 0.49 to 0.96 GULF pop_1; 0.04 to 0.51 GULF pop_2; GULF 0.0 to 0.09 pop_3
## LA (9) FL (4)

# GULF_10 - Majority GC 2 with TX/MX and a bit GC 1
## pure Full_K3_GULF
## 0.72 to 0.86 GULF pop_2; 0.07 to 0.23 GULF pop_3; 0.02 to 0.06 GULF pop_1
## TX (4)

# GULF_11 - Majority GC 2 with GC 1 and/or TX/MX
## almost pure Full_K3_GULF with trace amounts of ..MW
## 0.71 to 0.72 GULF pop_2; 0.13 to 0.20 GULF pop_3, 0.07 to 0.17 GULF pop_1
## FL (3)

# GULF_12 - Majority GC 2 with some TX/MX
## 0.72 to 0.85 Full_K3_GULF; 0.15 to 0.28 Full_K3_ATL
## 0.63 to 0.89 GULF pop_2, 0.11 to 0.37 GULF pop_3
## TX (18)

# GULF_13 - TX/MX with some GC 2
## 0.90 to 1.0 Full_K3_GULF; 0.0 to 0.10 Full_K3_MW
## 0.63 to 0.80 GULF pop_3; 0.20 to 0.37 GULF pop_2
## TX (18) SC (2)

# GULF_14 - roughly half TX/MX and GC 2 
## Almost pure Full_K3_GULF with trace amounts of Full_K3_MW
## 0.46 to 0.59 GULF pop_3; 0.41 to 0.54 GULF pop_2
## TX (5)

###
res_tab[, GULF_GROUPS := 'OTHER']
gulf_inds <- grep('GULF', res_tab$sub_grp)
res_tab[gulf_inds, GULF_GROUPS := res_tab$sub_grp[gulf_inds]]
res_tab$GULF_GROUPS <- as.factor(res_tab$GULF_GROUPS)
levels(res_tab$GULF_GROUPS) <- sort(unique(res_tab$GULF_GROUPS))

gulf_palette <- rainbow(15)
gulf_palette[15] <- 'black'
gg_gulf <- ggplot(res_tab, aes(x = full_pc01_stand, y = full_pc02_stand, 
                              color = GULF_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('All Samples, Gulf subgroup colors') + 
  scale_color_manual(values = gulf_palette)

gg_gulf

gg_gulf_sub <- ggplot(res_tab, aes(x = GULF_pc01_stand, y = GULF_pc02_stand,
                                  color = GULF_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('Gulf Samples') +
  scale_color_manual(values = gulf_palette)

gg_gulf_sub

pdf(gulf_plot_out, height = 6, width = 13)
gg_gulf + gg_gulf_sub
dev.off()

###################
# MW results
mw_groups <- all_groups[grep('MW', all_groups)]

for(pop_name in mw_groups){
  print(paste('#####', pop_name, '######'))
  print(paste(pop_name, 'K=3 MW'))
  print(summary(res_tab$full_admix_k3_MW[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 GULF'))
  print(summary(res_tab$full_admix_k3_GULF[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 ATL'))
  print(summary(res_tab$full_admix_k3_ATL[res_tab$sub_grp == pop_name]))
}

# MW ADMIXTURE K=3 Interpretations
# mw_tot_k3_pop_1 = MW-8X-West/Cosmopolitan + related 4X
# mw_tot_k3_pop_2 = most MW-4X
# mw_tot_k3_pop_3 = MW-8X-East; not present in any 4X samples

# MW Group Interpretations

# MW_01 = MW-8X-West/Cosmopolitan
## 0.72 to 1.0 Full_K3_MW; 0.0 to 0.28 Full_K3_GULF; trace ...ATL
## 0.86 to 1.0 MW pop_1
## AR (15) CO (2) IL (1) MI (1) MN (4) NE (5) NM (2) OH (1) OK (1) TX (3)

# MW_02 = MW-8X-Mixture; mixture of all 3 gene pools
## pure Full_K3_MW 
## 0.27 MW pop_1, 0.28 MW pop_2; 0.45 MW pop_3
## IN (1)

# MW_03 = main MW-8X-East group
## 0.93 to 1.0 Full_K3_MW; 0.0 to 0.07 ..ATL; trace ..GULF
## 0.91 to 1.0 MW pop_3; 0.0 to 0.09 MW pop_1
## PA (30) NY (14) WV (2) MA (2) VA (2) KY (1) OH (1)

# MW_04 = mainly MW-8X-East, but some 8X-West
## 0.97 to 1.0 Full_K3_MW; 0.0 to 0.03 Full_K3_GULF
## 0.32 to 0.80 MW pop_3; 0.14 to 0.68 MW pop_1
## PA (3) MO (1) NE (1) NY (2) ME (1) IA (1)

# MW_05 = MW-8X-East with some 4X
## pure Full_K3_MW
## 0.8 to 0.91 MW pop_3; 0.09 to 0.19 MW pop_2
## OH (6)

# MW_06 = MW-4X closely related to MW-8X West; 
## almost pure Full_K3_MW; trace amounts of Full_K3_GULF
## 0.8 to 1.0 pop_1; 0.0 to 0.21 pop_2
## MO (6) MI (2)

# MW_07 = MW-4X close to MW-8X West
## pure Full_K3_MW except occasional traces of Full_K3_GULF
## 0.52 to 0.69 MW pop_1; 0.31 to 0.48 MWpop_2
## MO (5) AR (4)

# MW_08 = MW-4X distant to MW-8X West
## pure Full_K3_MW
## 0.17 to 0.35 MW pop_1; 0.64 to 0.83 MW pop_2
## MO (18)

# MW_09 = MW-4X distant to MW-8X West
## pure Full_K3_MW except occasional traces of Full_K3_ATL
## 0.19 to 0.31 MW pop_1; 0.69 to 0.81 MW pop_2
## MO (2) NE (2) MN (2) WI (3) NY (1)

# MW_10 = MW-4X-GreatLakes
## pure Full_K3_MW
## pure MW pop_2
## MI (23) IL (16) WI (6) IN (3) 

# MW_11
# 0.96 Full_K3_MW; 0.04 Full_K3_ATL
## 0.94 pop_2
## IN (1)

# MW_12 = MW-4X-East1
## pure Full_K3_MW
## pure MW pop_2 
## OH (10)

# MW_13 = MW-4X-East2
## pure Full_K3_MW
## pure MW pop_2
## OH (6)

# MW_14 = MW-4X-UpperMW
## pure Full_K3_MW
## 0.85 to 1.0 MW pop_2; 0.0 to 0.15 MW pop_1
## MN (10) WI (2) IA (1) NY (1)

# NOTE: MW_10, MW_12, and MW_13 all have the same patterns from the
## ADMIXTURE results, but show some distinct clustering in PCA, 
## though only at the higher PCs (4 and 5)
# I could be conviced to group all of MW_4_X and MW_5_X together if we want to

res_tab[, MW_GROUPS := 'OTHER']
mw_inds <- grep('MW', res_tab$sub_grp)
res_tab[mw_inds, MW_GROUPS := res_tab$sub_grp[mw_inds]]
res_tab$MW_GROUPS <- as.factor(res_tab$MW_GROUPS)
levels(res_tab$MW_GROUPS) <- sort(unique(as.character(res_tab$MW_GROUPS)))

mw_palette <- rainbow(15)
mw_palette[15] <- 'black'
names(mw_palette) < sort(unique(as.character(res_tab$MW_GROUPS)))
gg_mw <- ggplot(res_tab, aes(x = full_pc01_stand, y = full_pc02_stand, 
                               color = MW_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('All Samples, MIDWEST subgroup colors') +
  scale_color_manual(values = mw_palette)

gg_mw

gg_mw_sub <- ggplot(res_tab, aes(x = MW_pc01_stand, y = MW_pc02_stand,
                                   color = MW_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('MIDWEST Samples') +
  scale_color_manual(values = mw_palette)

gg_mw_sub

pdf(mw_plot_out, height = 6, width = 13)
gg_mw + gg_mw_sub
dev.off()

###################
# Misfit results
mf_groups <- all_groups[grep('MF', all_groups)]

for(pop_name in mf_groups){
  print(paste('#####', pop_name, '######'))
  print(paste(pop_name, 'K=3 MW'))
  print(summary(res_tab$full_admix_k3_MW[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 GULF'))
  print(summary(res_tab$full_admix_k3_GULF[res_tab$sub_grp == pop_name]))
  print(paste(pop_name, 'K=3 ATL'))
  print(summary(res_tab$full_admix_k3_ATL[res_tab$sub_grp == pop_name]))
}

# Explanation of Misfit groups

# MF_01 - FL-4X
## p.54 to 0.54 Full_K3_ATL; 0.46 to 0.46 Full_K3_GULF
## 4X
## FL (2)

# MF_02 - NC-4X
## 0.51 to 0.61 Full_K3_ATL; 0.39 to 0.49 Full_K3_GULF
## 4X
## NC (2) 

# MF_03 - GA-8X-1
## 0.34 to 0.36 Full_K3_MW; 0.30 to 0.31 ..GULF; 0.33 to 0.35 ..ATL
## 8X
## GA (2) NC (1)

# MF_04 - NY/NJ
## 0.41 to 0.77 Full_K3_ATL; 0.23 to 0.52 ..MW; 0.0 to 0.31 ..GULF
## 4X (6) 8X (1)
## NJ (3) NY (4)

# MF_05 - GA-8X-2
## 0.47 to 0.48 Full_K3_MW; 0.29 to 0.30 ..ATL; 0.22 to 0.23 ..GULF
## 8X (3)
## GA (3)

# MF_06 - FL-8X
## 0.53 Full_K3_ATL; 0.40 ..GULF; 0.07 ..MW
## 8X (1)
## FL (1)

# MF_07 - GA-8X-3
## 0.44 to 0.93 Full_K3_MW; 0.04 to 0.32 ..GULF; 0.02 to 0.25 ..ATL
## 8X (7)
## GA (5) IL (1) TN (1)

# MF_08 - MS-8X
## 0.53 to 0.64 Full_K3_MW; 0.33 to 0.42 ..GULF; 0.0 to 0.07 ..ATL
## 8X (11) 6X (1) 
## AR (1) FL (1) MS (10)

res_tab[, MF_GROUPS := 'OTHER']
mf_inds <- grep('MF', res_tab$sub_grp)
res_tab[mf_inds, MF_GROUPS := res_tab$sub_grp[mf_inds]]
res_tab$MF_GROUPS <- as.factor(res_tab$MF_GROUPS)
levels(res_tab$MF_GROUPS) <- sort(unique(res_tab$MF_GROUPS))

mf_palette <- rainbow(9)
mf_palette[9] <- 'black'
gg_mf <- ggplot(res_tab, aes(x = full_pc01_stand, y = full_pc02_stand, 
                             color = MF_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('All Samples, "Misfit" subgroup colors') +
  scale_color_manual(values = mf_palette)

gg_mf

gg_mf_sub <- ggplot(res_tab, aes(x = MF_pc01_stand, y = MF_pc02_stand,
                                 color = MF_GROUPS, shape = ploidy)) +
  geom_point() +
  ggtitle('Misfit Samples') +
  scale_color_manual(values = mf_palette)

gg_mf_sub

pdf(mf_plot_out, height = 6, width = 13)
gg_mf + gg_mf_sub
dev.off()
