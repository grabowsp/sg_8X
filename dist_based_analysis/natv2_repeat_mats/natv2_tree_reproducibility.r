# Analysis of tree reproducibility for different genotypes

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
library(ape, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
library(phangorn, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')

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
  dip_list[[dip]] <- tmp_data[['euclidean_dist']]
}

tet_list <- list()
for(tet in seq(length(tet_files))){
  tmp_data <- readRDS(tet_files[tet])
  tet_list[[tet]] <- tmp_data[['euclidean_dist']]
}

poly_list <- list()
for(pf in seq(length(poly_files))){
  tmp_data <- readRDS(poly_files[pf])
  poly_list[[pf]] <- tmp_data[['euclidean_dist']]
}


### SET OUTPUT ###
res_file <- paste(dist_dir, 'natv2.100k.NJtree.reproducibility.results.rds',
  sep = '')

dip_tree_file <- paste(dist_dir, 'natv2.100k.dip_NJtree.list.v1.rds', sep = '')
tet_tree_file <- paste(dist_dir, 'natv2.100k.tet_NJtree.list.v1.rds', sep = '')
poly_tree_file <- paste(dist_dir, 'natv2.100k.poly_NJtree.list.v1.rds', 
  sep = '')
################

### Disomic genotypes ###
dip_tree_list <- lapply(dip_list, function(x) root(nj(x), 
                                    outgroup = 'ANCESTRAL_TET_GENO',
                                    resolve.root = T))
# for each sub-tree, calculate how many of the trees contain the same
#  edges
dip_reproduc_list <- list()
for(i in seq(length(dip_tree_list))){
  dip_reproduc_list[[i]] <- prop.clades(dip_tree_list[[i]], dip_tree_list, 
    rooted = T)
}

# Subtract 1 from mean and median because tested trees are part of total list
dip_mean_shared_vec <- (unlist(lapply(dip_reproduc_list, mean))-1)/(length(
                                dip_reproduc_list)-1)
mean(dip_mean_shared_vec)
# [1] 0.6951721
dip_median_shared_vec <- (unlist(lapply(dip_reproduc_list, median))-1)/(length(
                                dip_reproduc_list)-1)
mean(dip_median_shared_vec)
# [1] 0.9288889

### Tetrasomic genotypes ###
tet_tree_list <- lapply(tet_list, function(x) root(nj(x), 
                                    outgroup = 'ANCESTRAL_TET_GENO',
                                    resolve.root = T))
tet_reproduc_list <- list()
for(i in seq(length(tet_tree_list))){
  tet_reproduc_list[[i]] <- prop.clades(tet_tree_list[[i]], tet_tree_list, 
    rooted = T)
}

# Subtract 1 from mean and median because tested trees are part of total list
tet_mean_shared_vec <- (unlist(lapply(tet_reproduc_list, mean))-1)/(length(
                                tet_reproduc_list)-1)
mean(tet_mean_shared_vec)
# [1] 0.6905259
tet_median_shared_vec <- (unlist(lapply(tet_reproduc_list, median))-1)/(length(
                                tet_reproduc_list)-1)
mean(tet_median_shared_vec)
# [1] 0.9291919

### Polyploid genotypes ###
poly_tree_list <- lapply(poly_list, function(x) root(nj(x),
                                    outgroup = 'ANCESTRAL_TET_GENO',
                                    resolve.root = T))
poly_reproduc_list <- list()
for(i in seq(length(poly_tree_list))){
  poly_reproduc_list[[i]] <- prop.clades(poly_tree_list[[i]], poly_tree_list,
    rooted = T)
}

# Subtract 1 from mean and median because tested trees are part of total list
poly_mean_shared_vec <- (unlist(lapply(poly_reproduc_list, mean))-1)/(length(
                                poly_reproduc_list)-1)
mean(poly_mean_shared_vec)
# [1] 0.6915566
poly_median_shared_vec <- (unlist(lapply(poly_reproduc_list, median))-1)/(
                            length(poly_reproduc_list)-1)
mean(poly_median_shared_vec)
# [1] 0.9385859

### TAKE HOMES
# Diploid genotypes have highest (slightly) mean shared number of edges
# Polyploid genotypes have highest (slightly) median shared number of edges

### Save results
tot_reproduc_list <- list()
tot_reproduc_list[['dip_mean']] <- dip_mean_shared_vec
tot_reproduc_list[['dip_median']] <- dip_median_shared_vec
tot_reproduc_list[['tet_mean']] <- tet_mean_shared_vec
tot_reproduc_list[['tet_median']] <- tet_median_shared_vec
tot_reproduc_list[['poly_mean']] <- poly_mean_shared_vec
tot_reproduc_list[['poly_median']] <- poly_median_shared_vec

saveRDS(tot_reproduc_list, res_file)

saveRDS(dip_tree_list, dip_tree_file)
saveRDS(tet_tree_list, tet_tree_file)
saveRDS(poly_tree_list, poly_tree_file)

quit(save = 'no')

