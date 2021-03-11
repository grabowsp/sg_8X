# Generate lists of PCoA matrices from repeated distance matrices made using
#  different types of genotypes

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
#library(ape, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
#library(phangorn, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
#library(ggplot2)

### INPUT DATA ###
dist_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist/'

dip_files <- system(paste('ls ', dist_dir, '*diploid_DistMat.rds', sep = ''), 
  intern = T)
tet_files <- system(paste('ls ', dist_dir, '*tetraploid_DistMat.rds', 
  sep = ''), intern = T)
poly_files <- system(paste('ls ', dist_dir, '*polyploid_DistMat.rds', 
  sep = ''), intern = T)

dip_list <- list()
for(dip in seq(length(dip_files))){
  tmp_data <- readRDS(dip_files[dip])
  tmp_mat <- as.matrix(tmp_data[['euclidean_dist']])
  rm_inds <- grep('TET', rownames(tmp_mat))
  tmp_dist <- as.dist(tmp_mat[-rm_inds, -rm_inds])
  dip_list[[dip]] <- tmp_dist
}

tet_list <- list()
for(tet in seq(length(tet_files))){
  tmp_data <- readRDS(tet_files[tet])
  tmp_mat <- as.matrix(tmp_data[['euclidean_dist']])
  rm_inds <- grep('TET', rownames(tmp_mat))
  tmp_dist <- as.dist(tmp_mat[-rm_inds, -rm_inds])
  tet_list[[tet]] <- tmp_dist
}

poly_list <- list()
for(pf in seq(length(poly_files))){
  tmp_data <- readRDS(poly_files[pf])
  tmp_mat <- as.matrix(tmp_data[['euclidean_dist']])
  rm_inds <- grep('TET', rownames(tmp_mat))
  tmp_dist <- as.dist(tmp_mat[-rm_inds, -rm_inds])
  poly_list[[pf]] <- tmp_dist
}

### SET OUTPUT ###
res_file <- paste(dist_dir, 'natv2.100k.pcoa_mats.list.rds', sep = '')

################

### Make PCoA result lists

dip_pcoa_list <- lapply(dip_list, function(x) cmdscale(x, k = 20))
tet_pcoa_list <- lapply(tet_list, function(x) cmdscale(x, k = 20))
poly_pcoa_list <- lapply(poly_list, function(x) cmdscale(x, k = 20))

tot_list <- list()
tot_list[['dip_pcoa_list']] <- dip_pcoa_list
tot_list[['tet_pcoa_list']] <- tet_pcoa_list
tot_list[['poly_pcoa_list']] <- poly_pcoa_list

saveRDS(tot_list, file = res_file)

quit(save = 'no')


