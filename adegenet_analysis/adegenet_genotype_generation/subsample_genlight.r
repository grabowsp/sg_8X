# Script for subsampling the chromosome genlight objects to make a single
#   object including a random set of SNPs from all chromosomes

#module load python/3.7-anaconda-2019.07
#source activate r_adegenet_env

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###

library(adegenet)
library(parallel)

# get locations of ancillary github files
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.dir <- dirname(script.name)
gen.dir <- sub('/adegenet_analysis/adegenet_genotype_generation',
  '/general_r_tools', script.dir)
gen_function_file <- file.path(gen.dir, 'general_functions.r')
source(gen_function_file)

#gen_function_file <- '/global/homes/g/grabowsp/tools/sg_ha_effort/polyploid_genos/general_functions.r'
#source(gen_function_file)

### INPUT DATA ###

data_dir <- args[1]
#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps'

data_dir <- add_slash(data_dir)

file_sub <- args[2]
#file_sub <- '*geosamps.genlight.rds'

file_ls <- system(paste('ls ', data_dir, file_sub, sep = ''), intern = T)

#sub_in <- 'Chr01K.polyploid.CDS.geosamps.genlight.rds'

tmp_gl <- readRDS(file_ls[1])

### SET OUTPUT ###
out_short <- args[3]
#out_short <- 'Combo.595K.polyploid.CDS.geosamps.genlight.rds'
out_full <- paste(data_dir, out_short, sep = '')

### SET VARIABLES ###
per_subsamp <- as.numeric(args[4])
#per_subsamp <- 0.07

##############

keep_inds <- sort(sample(seq(nLoc(tmp_gl)), size = nLoc(tmp_gl) * per_subsamp))

tot_gl <- tmp_gl[, keep_inds]

for(fn in c(2:length(file_ls))){
  tmp_gl <- readRDS(file_ls[fn])
  keep_inds <- sort(sample(seq(nLoc(tmp_gl)), 
    size = nLoc(tmp_gl) * per_subsamp))
  sub_gl <- tmp_gl[, keep_inds]
  tot_gl <- cbind(tot_gl, sub_gl)
  print(fn)
}

saveRDS(tot_gl, out_full)

quit(save = 'no')

