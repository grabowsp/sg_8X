# Script for comparing structure and ADMIXTURE results for 'geobig' samples

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###

data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/'

admix_k2_res_file <- paste(data_dir, 'GW_50k_geobig.2.results.txt', sep = '')
admix_k2 <- fread(admix_k2_res_file, header = F)

admix_k3_res_file <- paste(data_dir, 'GW_50k_geobig.3.results.txt', sep = '')
admix_k3 <- fread(admix_k3_res_file, header = F)

admix_k4_res_file <- paste(data_dir, 'GW_50k_geobig.4.results.txt', sep = '')
admix_k4 <- fread(admix_k4_res_file, header = F)

struc_k2_res_file <- paste(data_dir, 
  'geobig_alltet_2_combo_k2.clumpp.processed', sep = '')
struc_k2 <- fread(struc_k2_res_file, header = F)

struc_k3_res_file <- paste(data_dir, 
  'geobig_alltet_3_combo_k3.clumpp.processed', sep = '')
struc_k3 <- fread(struc_k3_res_file, header = F)

struc_k4_res_file <- paste(data_dir, 
  'geobig_alltet_4_combo_k4.clumpp.processed', sep = '')
struc_k4 <- fread(struc_k4_res_file, header = F)


#######

struc_k2_reorder <- c()
for(i in seq(nrow(admix_k2))){
  tmp_ind <- which(struc_k2$V1 == admix_k2$V1[i])
  struc_k2_reorder <- c(struc_k2_reorder, tmp_ind)
}

struc_k2_2 <- struc_k2[struc_k2_reorder, ]

cor(struc_k2_2$V2, admix_k2$V3)
# [1] 0.9992777
# essentially perferctly correlated

struc_k3_reorder <- c()
for(i in seq(nrow(admix_k3))){
  tmp_ind <- which(struc_k3$V1 == admix_k3$V1[i])
  struc_k3_reorder <- c(struc_k3_reorder, tmp_ind)
}

struc_k3_2 <- struc_k3[struc_k3_reorder, ]

cor(struc_k3_2$V2, admix_k3$V3)
#[1] 0.9989913
cor(struc_k3_2$V3, admix_k3$V4)
# [1] 0.9993618
cor(struc_k3_2$V4, admix_k3$V2)
# [1] 0.9994477

struc_k4_reorder <- c()
for(i in seq(nrow(admix_k4))){
  tmp_ind <- which(struc_k4$V1 == admix_k4$V1[i])
  struc_k4_reorder <- c(struc_k4_reorder, tmp_ind)
}

struc_k4_2 <- struc_k4[struc_k4_reorder, ]

cor(struc_k4_2$V2, admix_k4$V3)
# [1] 0.8913076
cor(struc_k4_2$V3, admix_k4$V2)
# [1] 0.9983528
cor(struc_k4_2$V4, admix_k4$V4)
# [1] 0.9125295
cor(struc_k4_2$V5, admix_k4$V5)
# [1] 0.2569367

# since structure and admixture don't pull out the same 4th group, seems like
#  a good reason to stick with k=3 when using the full geographic patterns
quit(save = 'no')

