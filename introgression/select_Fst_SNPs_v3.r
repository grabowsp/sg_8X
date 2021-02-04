# Script for selecting and consolidating the high-Fst SNPs to be used
#  for introgression analysis

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/'

AG_data_files <- system(paste('ls ', res_dir, 'ATLANTICvGULF*weir.fst', 
  sep = ''), intern = T)
AM_data_files <- system(paste('ls ', res_dir, 'ATLANTICvMW*weir.fst', 
  sep = ''), intern = T)
GM_data_files <- system(paste('ls ', res_dir, 'GULFvMW*weir.fst',  
  sep = ''), intern = T)

### SET OUTPUT ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/'

AG_out_file <- paste(out_dir, 'AtlanticVsGulf_hiFst_SNPs_v3.txt', sep = '')
MG_out_file <- paste(out_dir, 'GulfVsMW_hiFst_SNPs_v3.txt', sep = '')
MA_out_file <- paste(out_dir, 'AtlanticVsMW_hiFst_SNPs_v3.txt', sep = '')

#out_file <- paste(res_dir, 'hi_fst_SNPS_v2.txt', sep = '')

### SET VARIABLES ###
fst_cut <- 0.5

#########################
hi_fst_list <- list()

test <- fread(AG_data_files[1])

get_fst_snps <- function(file_vec, fst_cut){
  hi_fst_list <- list()
  for(fdf in file_vec){
    tmp_data <- fread(fdf)
    tmp_chr <- rev(unlist(strsplit(fdf, split = '_')))[2]
    tmp_subdata <- tmp_data[WEIR_AND_COCKERHAM_FST >= fst_cut, ]
    hi_fst_list[[tmp_chr]] <- tmp_subdata
  }
  tot_fst_tab <- rbindlist(hi_fst_list)
  return(tot_fst_tab)
}

AG_hi_Fst <- get_fst_snps(file_vec = AG_data_files, fst_cut = fst_cut)
# 67057 with Fst > 0.5

MG_hi_Fst_0 <- get_fst_snps(file_vec = GM_data_files, fst_cut = fst_cut)
# 150649 with Fst > 0.5
MG_sub_inds <- sort(sample(seq(nrow(MG_hi_Fst_0)), size = 1e5))
MG_hi_Fst <- MG_hi_Fst_0[MG_sub_inds, ]

MA_hi_Fst_0 <- get_fst_snps(file_vec = AM_data_files, fst_cut = fst_cut)
# 156090 with Fst > 0.5
MA_sub_inds <- sort(sample(seq(nrow(MA_hi_Fst_0)), size = 1e5))
MA_hi_Fst <- MA_hi_Fst_0[MA_sub_inds]

keycol <- c('CHROM', 'POS')

AG_out_tab <- setorderv(AG_hi_Fst, keycol)[, ..keycol]

MG_out_tab <- setorderv(MG_hi_Fst, keycol)[, ..keycol]

MA_out_tab <- setorderv(MA_hi_Fst, keycol)[, ..keycol]

fwrite(AG_out_tab, AG_out_file, sep = '\t', col.names = F)
fwrite(MG_out_tab, MG_out_file, sep = '\t', col.names = F)
fwrite(MA_out_tab, MA_out_file, sep = '\t', col.names = F)

quit(save = 'no')

