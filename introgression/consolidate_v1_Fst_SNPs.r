# Script for selecting and consolidating the high-Fst SNPs to be used
#  for introgression analysis

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/'

data_files <- system(paste('ls ', res_dir, '*weir.fst', sep = ''), intern = T)

### SET OUTPUT ###
out_file <- paste(res_dir, 'hi_fst_SNPS_v1.txt', sep = '')

### SET VARIABLES ###
fst_cut <- 0.75

#########################
hi_fst_list <- list()

for(fdf in data_files){
  tmp_data <- fread(fdf)
  tmp_filename <- sub(res_dir, '', fdf)
  tmp_comp <- unlist(strsplit(tmp_filename, split = '_'))[1]
  tmp_subdata <- tmp_data[WEIR_AND_COCKERHAM_FST >= fst_cut, ]
  tmp_subdata[, COMP := tmp_comp]  
  hi_fst_list[[tmp_filename]] <- tmp_subdata
}

combo_fst_list <- list()

for(i in seq(ncol(hi_fst_list[[1]]))){
  combo_fst_list[[i]] <- unlist(lapply(hi_fst_list, 
    function(x) x[, i, with = F]))
}

hi_fst_big <- data.table(CHROM = combo_fst_list[[1]], POS = combo_fst_list[[2]],
  SNP_NAME = paste(combo_fst_list[[1]], combo_fst_list[[2]], sep = '_'), 
  FST_VAL = combo_fst_list[[3]], COMPARISON = combo_fst_list[[4]])

table(hi_fst_big$COMPARISON[hi_fst_big$SNP_NAME %in% 
  names(which(table(hi_fst_big$SNP_NAME) == 1))])
#         ALL3 ATLANTICvGULF   ATLANTICvMW       GULFvMW 
#            5          3217         19791         15711 
# Including Fst calculated using all 3 groups only added 5 unique SNPs

hi_fst_unique <- unique(hi_fst_big, by = 'SNP_NAME')

keycol <- c('CHROM', 'POS')
setorderv(hi_fst_unique, keycol)

out_tab <- hi_fst_unique[, c('CHROM', 'POS')]

fwrite(out_tab, out_file, sep = '\t', col.names = F)

quit(save = 'no')

