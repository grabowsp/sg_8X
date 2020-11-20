# Code for quick geosamps analysis using adegenet

#module load python/3.7-anaconda-2019.07
#source activate adegenet_2_env

args = commandArgs(trailingOnly = TRUE)

### LOAD PACKAGES ###

library(adegenet)
library(parallel)

### INPUT DATA ###
data_file <- args[1]
#data_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.rds'

tot_gl <- readRDS(data_file)

### SET OUTPUT ###
out_file <- gsub('.rds', '.PCAresults.rds', data_file)

### SET VARIABLES ###

############
n_eig <- nInd(tot_gl) - 1

tot_pca <- glPca(tot_gl, nf = n_eig, loadings = F, alleleAsUnit = F, useC = F)

saveRDS(tot_pca, out_file)

quit(save = 'no')

