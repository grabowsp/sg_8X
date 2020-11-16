# Script to get the total number of SNPs in the genlight object for each
#   chromosome
#  Can be used to calculate at what rate to subsample SNPs from each file

#module load python/3.7-anaconda-2019.07
#source activate adegenet_2_env

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###

library(adegenet)
#library(parallel)

# get locations of ancillary github files
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.dir <- dirname(script.name)
gen.dir <- sub('/adegenet_analysis/adegenet_genotype_generation',
  '/general_r_tools', script.dir)
gen_function_file <- file.path(gen.dir, 'general_functions.r')
source(gen_function_file)

#gen_function_file <- '/global/homes/g/grabowsp/tools/sg_8X/general_r_tools/general_functions.r'
#source(gen_function_file)

### INPUT DATA ###

data_dir <- args[1]
#data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs'

data_dir <- add_slash(data_dir)

file_sub <- args[2]
#file_sub <- '*geobig.genlight.rds'

file_ls <- system(paste('ls ', data_dir, file_sub, sep = ''), intern = T)

#sub_in <- 'Chr01K.tetrasomic.CDS.geobig.genlight.rds'
#tmp_gl <- readRDS(file_ls[1])

### SET OUTPUT ###
out_short <- args[3]
#out_short <- 'geobig.SNPcount.txt'
out_suf <- '.SNPcount.txt'
out_full <- paste(data_dir, out_short, out_suf, sep = '')

### SET VARIABLES ###


##############

nSNP_vec <- c()

for(i in seq(length(file_ls))){
  print(i)
  tmp_gl <- readRDS(file_ls[i])
  tmp_nSNP <- nLoc(tmp_gl)
  nSNP_vec <- c(nSNP_vec, tmp_nSNP)
}

nSNP_df <- data.frame(file = file_ls, nSNPs = nSNP_vec, stringsAsFactors = F)

write.table(nSNP_df, file = out_full, quote = F, sep = '\t', row.names = F,
  col.names = T)

quit(save = 'no')


