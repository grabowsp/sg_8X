# bash
# source activate R_analysis

library(data.table)

meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_file)

MW_train_libs_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/Midwest_libnames.txt'
MW_train_libs <- fread(MW_train_libs_file, header = F)

ATL_train_libs_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/Atlantic_libnames.txt'
ATL_train_libs <- fread(ATL_train_libs_file, header = F)

GULF_train_libs_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/Gulf_libnames.txt'
GULF_train_libs <- fread(GULF_train_libs_file, header = F)


length(intersect(MW_train_libs$V1, res_tab$samp_name))
# 37

length(intersect(ATL_train_libs$V1, res_tab$samp_name))
# 38

length(intersect(GULF_train_libs$V1, res_tab$samp_name))
# 24

res_tab[which(res_tab$samp_name %in% MW_train_libs$V1), .N, by = subgrp_v2]
# 1:     MW_14  8
# 2:     MW_09  5
# 3:     MW_10 16
# 4:     MW_08  4
# 5:     MW_13  3
# 6:     MW_12  1

res_tab[grep('MW_', res_tab$subgrp_v2), .N, by = subgrp_v2]
    subgrp_v2  N
 1:     MW_10 48
 2:     MW_14 14
# 3:     MW_09 10
# 4:     MW_11  1
# 5:     MW_01 21
# 6:     MW_08 18
# 7:     MW_03 52
# 8:     MW_07  9
# 9:  MW_01_hi 14
10:     MW_12 10
11:     MW_13  6
#12:     MW_05  6
#13:     MW_04 10
#14:     MW_02  1
#15:     MW_06  8

# sum(c(48,14,10,6)) = 78 samples from "non-introgressed clade" for next
#   the training sets...

sum(res_tab$full_admix_k3_MW > 0.99)
# 190

res_tab[full_admix_k3_MW > 0.99, .N, by = subgrp_v2]

summary(res_tab[full_admix_k3_MW > 0.99, full_pc01_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8416  0.9001  0.9375  0.9378  0.9796  1.0000 

summary(res_tab[grep('MW_', res_tab$subgrp_v2), full_pc01_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7006  0.8926  0.9126  0.9169  0.9759  1.0000 

summary(res_tab[samp_name %in% MW_train_libs$V1, full_pc01_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.8975  0.9670  0.9766  0.9726  0.9824  0.9929

res_tab[full_pc01_stand > 0.9759, .N, by = subgrp_v2]
#   subgrp_v2  N
#1:     MW_10 39
#2:     MW_14  6
#3:     MW_12  8
#4:     MW_13  4

mw_new_train_samps <- res_tab[full_pc01_stand > 0.9759, samp_name]

samp_meta[VCF_NAME %in% mw_new_train_samps, .N, by = STATE]
samp_meta[VCF_NAME %in% mw_new_train_samps, .N, by = UNI_ACC]


summary(res_tab[full_pc01_stand > 0.9759, full_admix_k3_MW])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      1       1       1       1       1       1

setdiff(MW_train_libs$V1, res_tab$samp_name)
# [1] "J036.A" "J636.B" "DAC6"
# all 3 are part of Dacotah cluser, which is similar to MW_09

# so 12 of 40 training samples from groups that show some level of Lowland
#  introgression or Dacotah which seems like it might, too

sort(res_tab[subgrp_v2 == 'MW_09', samp_name])

#####

res_tab[which(res_tab$samp_name %in% ATL_train_libs$V1), .N, by = subgrp_v2]
   subgrp_v2  N
1:    ATL_10 15
2:    ATL_05  4
3:    ATL_04  5
4:    ATL_08 13
5:    ATL_06  1

res_tab[full_admix_k3_ATL > 0.99, .N, by = subgrp_v2]

summary(res_tab[full_admix_k3_ATL > 0.99, full_pc01_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02577 0.03259 0.03658 0.04451 0.13020
summary(res_tab[grep('ATL_', res_tab$subgrp_v2), full_pc01_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.02714 0.03562 0.04303 0.05106 0.14824
summary(res_tab[samp_name %in% ATL_train_libs$V1, full_pc01_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.02266 0.02977 0.03657 0.04112 0.04189 0.10617 

res_tab[full_pc01_stand < 0.02714, .N, by = subgrp_v2]
# 71 4X from 8 different subgroups

setdiff(ATL_train_libs$V1, res_tab$samp_name)
# [1] "J598.A" "J602.A" - these were removed because probably a sample switch,
#  but shouldn't have affected the results

#######

res_tab[which(res_tab$samp_name %in% GULF_train_libs$V1), .N, by = subgrp_v2]
   subgrp_v2 N
1:   GULF_10 3
2:   GULF_07 4
3:   GULF_01 2
4:   GULF_02 3
5:   GULF_13 5
6:   GULF_14 4
7:   GULF_06 3

res_tab[full_admix_k3_GULF > 0.99, .N, by = subgrp_v2]

summary(res_tab[full_admix_k3_GULF > 0.99, full_pc02_stand])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.03104 0.07348 0.09115 0.11080 0.32095 
summary(res_tab[grep('GULF_', res_tab$subgrp_v2), full_pc02_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.07751 0.22088 0.19552 0.31226 0.41773
summary(res_tab[samp_name %in% GULF_train_libs$V1, full_pc02_stand])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01702 0.03211 0.07327 0.07616 0.09560 0.18996

res_tab[full_pc02_stand < 0.07751, .N, by = subgrp_v2]
res_tab[full_pc02_stand < 0.1, .N, by = ploidy]

sort(res_tab[subgrp_v2 == 'GULF_10', samp_name])

