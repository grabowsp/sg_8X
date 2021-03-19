# Compare Sujan's list to libraries for Adam

### LOAD PACKAGES ###
library(tidyverse)
library('data.table')

### INPUT DATA ###

libs_for_Adam <- '/Users/grabowsk/data/sg_8X_analysis/kmer_4x8x_libs_v1.csv'
test_libs <- fread(libs_for_Adam)

grp_3_file <- '/Users/grabowsk/data/sg_8X_analysis/kmer_grp3_comp_libs_v1.csv'
grp3_libs <- fread(grp_3_file)

sujan_file_1 <- '/Users/grabowsk/Downloads/Pvirgatum_Paul_Adam_preps.txt'
sfloc_1 <- fread(sujan_file_1, header = F)
sfloc_2 <- unique(sfloc_1$V1)
sfloc_3 <- sfloc_2[-grep('bicolor', sfloc_2)]

sujan_file_2 <- '/Users/grabowsk/Downloads/Pvirgatum_Paul_Adam_preps_notAvail.txt'
miss_libs_1 <- fread(sujan_file_2, header = F)
miss_libs_2 <- unique(miss_libs_1$V1)

####
test_libs[, LIB_LOC := as.character(NA)]
miss_libs <- which(test_libs$LIB %in% miss_libs_2)

test_libs[miss_libs, LIB_LOC := 'PREP_UNAVAILABLE']

for(i in seq(nrow(test_libs))){
  if(is.na(test_libs$LIB_LOC[i])){
    tmp_inds <- grep(test_libs$LIB[i], sfloc_3)
    test_libs[i, LIB_LOC := sfloc_3[tmp_inds]]
  }
}

grp3_libs[, LIB_LOC := as.character(NA)]
g3_miss_libs <- which(grp3_libs$LIB %in% miss_libs_2)

grp3_libs[g3_miss_libs, LIB_LOC := 'PREP_UNAVAILABLE']

for(j in seq(nrow(grp3_libs))){
  if(is.na(grp3_libs$LIB_LOC[j])){
    tmp_inds <- grep(grp3_libs$LIB[j], sfloc_3)
    grp3_libs[j, LIB_LOC := sfloc_3[tmp_inds]]
  }
}
