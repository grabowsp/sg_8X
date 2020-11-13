# module load python/3.7-anaconda-2019.07
# source activate r_adegenet_env

library(adegenet)
library(parallel)

data_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/Chr01N.tetrasomic.CDS.allsamps.genlight.rds'

data <- readRDS(data_file)

tet_inds <- which(ploidy(data) == 2)
oct_inds <- which(ploidy(data) == 4)

# figure out number of samples with het genotype
test <- lapply(lapply(seploc(data, n.block = 10, parallel = T), function(x) 
  as.matrix(x[tet_inds,])), function(y) 
  apply(y, 2, function(z) sum(z == 1, na.rm = T)))

# count NAs at each SNP
test_isna <- lapply(lapply(seploc(data, n.block = 10, parallel = T), function(x)
  as.matrix(x[tet_inds,])), function(y) 
  apply(y, 2, function(z) sum(is.na(z))))

# calculate n samps for each SNP
test_nsamps <- lapply(test_isna, function(x) length(tet_inds) - x)

# calculate the cutoff for minimum n het samps; if n het samps is above
#   this number, then SNP het is too high and should be removed
test_cutoff <- lapply(test_nsamps, function(x) 
  (qbinom(0.001, size = x, prob = 0.5, lower.tail = F)))

test_excess_het <- list()
for(i in seq(length(test))){
  test_excess_het[[i]] <- test[[i]] - test_cutoff[[i]]
}

sum(unlist(lapply(test_excess_het, function(x) sum(x > 0))))
# 0


test_2 <- lapply(lapply(seploc(data, n.block = 10, parallel = T), function(x)
  as.matrix(x[oct_inds,])), function(y) 
  apply(y, 2, function(z) sum(z == 1, na.rm = T)))

test_isna_2 <- lapply(lapply(seploc(data, n.block = 10, parallel = T), function(x)
  as.matrix(x[oct_inds,])), function(y)
  apply(y, 2, function(z) sum(is.na(z))))

test_nsamps_2 <- lapply(test_isna_2, function(x) length(oct_inds) - x)

test_cutoff_2 <- lapply(test_nsamps_2, function(x) 
  (qbinom(0.001, size = x, prob = 0.875, lower.tail = F)))

test_excess_het_2 <- list()
for(i in seq(length(test_2))){
  test_excess_het_2[[i]] <- test_2[[i]] - test_cutoff_2[[i]]
}

sum(unlist(lapply(test_excess_het_2, function(x) sum(x > 0))))
# 0


