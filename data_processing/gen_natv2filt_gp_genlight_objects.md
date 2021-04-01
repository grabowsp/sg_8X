# Generate the genlight objects for each of the natv2filt temp genepool groups

# R script
```
# module load python/3.7-anaconda-2019.10
# source activate adegenet_2_env

### LOAD PACKAGES ###
library(adegenet)
library(data.table)
library(parallel)

### LOAD INPUTS ###
genlight_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/GW.100kSNPs.tetrasomic.CDS.natv2filt.genlight.rds'
tot_genlight <- readRDS(genlight_in)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

mw_main_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_MW_main_names.txt'
mw_main_names <- fread(mw_main_file, header = F)

mw_tot_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_MW_total_names.txt'
mw_tot_names <- fread(mw_tot_file, header = F)

gulf_main_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_GULF_main_names.txt'
gulf_main_names <- fread(gulf_main_file, header = F)

gulf_tot_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_GULF_total_names.txt'
gulf_tot_names <- fread(gulf_tot_file, header = F)

atl_main_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_ATL_main_names.txt'
atl_main_names <- fread(atl_main_file, header = F)

atl_tot_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_ATL_total_names.txt'
atl_tot_names <- fread(atl_tot_file, header = F)

misfit_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_MISFIT_names.txt'
misfit_names <- fread(misfit_file, header = F)

### SET OUTPUTS ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/'

out_pre <- 'GW.100kSNPs.tetrasomic.CDS.natv2filt.'
out_suf <- '.genlight.rds'

mw_main_out <- paste(out_dir, out_pre, 'MW_main', out_suf, sep = '')
mw_tot_out <- paste(out_dir, out_pre, 'MW_tot', out_suf, sep = '')

gulf_main_out <- paste(out_dir, out_pre, 'GULF_main', out_suf, sep = '')
gulf_tot_out <- paste(out_dir, out_pre, 'GULF_tot', out_suf, sep = '')

atl_main_out <- paste(out_dir, out_pre, 'ATL_main', out_suf, sep = '')
atl_tot_out <- paste(out_dir, out_pre, 'ATL_tot', out_suf, sep = '')

misfit_out <- paste(out_dir, out_pre, 'MISFIT', out_suf, sep = '')

### Analysis Functions ###

make_sub_genlight <- function(tot_gl_object, sub_names, nMA_cut = 5){
  tmp_gl <- tot_gl_object[which(indNames(tot_gl_object) %in% sub_names), ]
  tmp_novar <- which(glSum(tmp_gl) < nMA_cut | 
    glSum(tmp_gl) > (sum(ploidy(tmp_gl))-5))
  res_gl <- tmp_gl[, -tmp_novar]
  return(res_gl)
}

##################

mw_main_gl <- make_sub_genlight(tot_gl_object = tot_genlight, 
  sub_names = mw_main_names$V1)
nLoc(mw_main_gl)
# [1] 48924
saveRDS(mw_main_gl, file = mw_main_out)

mw_tot_gl <- make_sub_genlight(tot_gl_object = tot_genlight,
  sub_names = mw_tot_names$V1)
nLoc(mw_tot_gl)
# 53019
saveRDS(mw_tot_gl, file = mw_tot_out)

gulf_main_gl <- make_sub_genlight(tot_gl_object = tot_genlight,
  sub_names = gulf_main_names$V1)
nLoc(gulf_main_gl)
# [1] 31385
saveRDS(gulf_main_gl, file = gulf_main_out)

gulf_tot_gl <- make_sub_genlight(tot_gl_object = tot_genlight,
  sub_names = gulf_tot_names$V1)
nLoc(gulf_tot_gl)
# [1] 50477
saveRDS(gulf_tot_gl, file = gulf_tot_out)

atl_main_gl <- make_sub_genlight(tot_gl_object = tot_genlight,
  sub_names = atl_main_names$V1)
nLoc(atl_main_gl)
# [1] 49975
saveRDS(atl_main_gl, file = atl_main_out)

atl_tot_gl <- make_sub_genlight(tot_gl_object = tot_genlight,
  sub_names = atl_tot_names$V1)
nLoc(atl_tot_gl)
# [1] 53965
saveRDS(atl_tot_gl, file = atl_tot_out)

misfit_gl <- make_sub_genlight(tot_gl_object = tot_genlight,
  sub_names = misfit_names$V1, nMA_cut = 3)
nLoc(misfit_gl)
# 47038
saveRDS(misfit_gl, file = misfit_out)

quit(save = 'no')
```
