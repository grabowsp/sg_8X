# Script for selecting and consolidating the high-Fst SNPs to be used
#  for introgression analysis

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/'

AG_data_files <- system(paste('ls ', res_dir, 'ATLANTICvGULF*weir.fst', 
  sep = ''), intern = T)
AM_data_files <- system(paste('ls ', res_dir, 'ATLANTICvMW*weir.fst', 
  sep = ''), intern = T)
GM_data_files <- system(paste('ls ', res_dir, 'GULFvMW*weir.fst',  
  sep = ''), intern = T)

### SET OUTPUT ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/'

AG_out_file_1 <- paste(out_dir, 
  'AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.txt', sep = '')
AG_out_file_2 <- paste(out_dir, 
  'AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.txt', sep = '')

MG_out_file_1 <- paste(out_dir, 'GulfVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt', 
  sep = '')
MG_out_file_2 <- paste(out_dir, 'GulfVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt',
  sep = '')

MA_out_file_1 <- paste(out_dir, 
  'AtlanticVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt', sep = '')
MA_out_file_2 <- paste(out_dir,
  'AtlanticVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt', sep = '')

#out_file <- paste(res_dir, 'hi_fst_SNPS_v2.txt', sep = '')

### SET VARIABLES ###
fst_cut_1 <- 0.5
fst_cut_2 <- 1.0
#########################
#hi_fst_list <- list()

#test <- fread(AG_data_files[1])

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

AG_hi_Fst_1 <- get_fst_snps(file_vec = AG_data_files, fst_cut = fst_cut_1)
# 65318 with Fst >= 0.5
AG_hi_Fst_2 <- get_fst_snps(file_vec = AG_data_files, fst_cut = fst_cut_2)
# 569 with Fst = 1

MG_hi_Fst_1 <- get_fst_snps(file_vec = GM_data_files, fst_cut = fst_cut_1)
# 149467 with Fst >= 0.5
MG_hi_Fst_2 <- get_fst_snps(file_vec = GM_data_files, fst_cut = fst_cut_2)
# 3848 with Fst = 1

MA_hi_Fst_1 <- get_fst_snps(file_vec = AM_data_files, fst_cut = fst_cut_1)
# 167400
MA_hi_Fst_2 <- get_fst_snps(file_vec = AM_data_files, fst_cut = fst_cut_2)
# 8435 with Fst = 1

# immediately, I'm more confident in the MW_vs_ATL than MW_vs_GULF SNPs 

keycol <- c('CHROM', 'POS')

AG_out_tab_1 <- setorderv(AG_hi_Fst_1, keycol)[, ..keycol]
AG_out_tab_2 <- setorderv(AG_hi_Fst_2, keycol)[, ..keycol]

MG_out_tab_1 <- setorderv(MG_hi_Fst_1, keycol)[, ..keycol]
MG_out_tab_2 <- setorderv(MG_hi_Fst_2, keycol)[, ..keycol]

MA_out_tab_1 <- setorderv(MA_hi_Fst_1, keycol)[, ..keycol]
MA_out_tab_2 <- setorderv(MA_hi_Fst_2, keycol)[, ..keycol]

fwrite(AG_out_tab_1, AG_out_file_1, sep = '\t', col.names = F)
fwrite(AG_out_tab_2, AG_out_file_2, sep = '\t', col.names = F)

fwrite(MG_out_tab_1, MG_out_file_1, sep = '\t', col.names = F)
fwrite(MG_out_tab_2, MG_out_file_2, sep = '\t', col.names = F)

fwrite(MA_out_tab_1, MA_out_file_1, sep = '\t', col.names = F)
fwrite(MA_out_tab_2, MA_out_file_2, sep = '\t', col.names = F)

quit(save = 'no')

