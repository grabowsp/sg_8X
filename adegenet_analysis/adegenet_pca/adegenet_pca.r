# Code for quick geosamps analysis using adegenet

#module load python/3.7-anaconda-2019.07
#source activate adegenet_2_env

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

# get locations of ancillary github files
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.dir <- dirname(script.name)

gen.dir <- sub('/adegenet_analysis/adegenet_pca',
  '/general_r_tools', script.dir)
gen_function_file <- file.path(gen.dir, 'general_functions.r')
source(gen_function_file)

### LOAD PACKAGES ###

library(adegenet)
library(parallel)

### INPUT DATA ###
data_file <- args[1]
#data_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/genlight_objs/Combo.595K.polyploid.CDS.geosamps.genlight.rds'

tot_gl <- readRDS(data_file)

### SET OUTPUT ###
out_dir <- args[2]

out_dir <- add_slash(out_dir)

base_in <- basename(data_file)

out_file_short <- gsub('.rds', '.PCAresults.rds', base_in)

out_file <- paste(out_dir, out_file_short, sep = '')

### SET VARIABLES ###

############
n_eig <- nInd(tot_gl) - 1

tot_pca <- glPca(tot_gl, nf = n_eig, loadings = F, alleleAsUnit = F, useC = F)

saveRDS(tot_pca, out_file)

quit(save = 'no')

