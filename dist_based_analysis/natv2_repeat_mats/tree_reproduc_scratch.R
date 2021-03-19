# Use phangorn package to make NeighborNet figure

### INSTALL PACKAGES ###
# install.packages('phangorn', dependencies=T)
# install.packages('devtools', dependencies=T)
# library(devtools)
# install_github('KlausVigo/phangorn')

### LOAD PACKAGES ###
library(phangorn)

### INPUT DATA ###

# dist_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2_ancest_100k_tet_DistMat.rds'
# dist_data <- readRDS(dist_file)

nj_tree_file <- 'PLACE/HOLDER'
tree_list <- readRDS(nj_tree_file)

### SET VARIABLE ###
n_tests <- 10

########

# test <- prop.part(tree_list, check.labels = T)  

reproduc_list <- list()
for(i in seq(n_tests)){
  reproduc_list[[i]] <- prop.clades(tree_list[[i]], tree_list, rooted = T)
}

##########

nexus_file <- '/Users/grabowsk/Desktop/GW.allsamps.ancestral.sub.min4.nexus'
test <- read.nexus.data(nexus_file)

test_1 <- phyDat(test)

tree <- nj(dist.hamming(test_1))

bs <- bootstrap.phyDat(test_1, FUN = function(x)nj(dist.hamming(x)), bs = 100)

####

tree_list <- lapply(dist_list, function(x) root(nj(x), 
                                                outgroup = 'ANCESTRAL_TET_GENO',
                                                resolve.root = T))

test <- prop.part(tree_list, check.labels = T)  
  
test_2 <- prop.clades(tree_list[[1]], tree_list, rooted = T)

nn_1 <- neighborNet(dist_data[['euclidean_dist']])


tree <- root(nj(dist_data[['euclidean_dist']]), 
             outgroup = 'ANCESTRAL_TET_GENO',
             resolve.root = T)

bs <- bootstrap.phyDat(dist_data[['euclidean_dist']], FUN = function(x) nj(x),
                       bs=100)
