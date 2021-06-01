# Generate sample name sets for further MW introgression analysis

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
samp_info_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
samp_info <- fread(samp_info_file)

### SET OUTPUT ###
out_dir <- '/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/mw_introgression_sets/'

mw8X_all_out <- paste(out_dir, 'MW8X_all_names.txt', sep = '')
mw8X_low_out <- paste(out_dir, 'MW8X_reg_names.txt', sep = '')
mw8X_hi_out <- paste(out_dir, 'MW8X_hi_names.txt', sep = '')
mw_01_out <- paste(out_dir, 'MW_01_names.txt', sep = '')
mw_03_out <- paste(out_dir, 'MW_03_names.txt', sep = '')

mw4X_intro_out <- paste(out_dir, 'MW4X_intro_names.txt', sep = '')
########

mw8X_all_names <- samp_info[grep('MW_8X', samp_info$MW_grp_1), list(samp_name)]
fwrite(mw8X_all_names, file = mw8X_all_out, col.names = F)

mw8X_low_names <- samp_info[which(samp_info$MW_grp_1 == 'MW_8X'), 
  list(samp_name)]
fwrite(mw8X_low_names, file = mw8X_low_out, col.names = F)

mw8X_hi_names <- samp_info[which(samp_info$MW_grp_1 == 'MW_8X_hi'), 
  list(samp_name)]
fwrite(mw8X_hi_names, file = mw8X_hi_out, col.names = F)

mw_01_names <- samp_info[which(samp_info$subgrp_v2 == 'MW_01'), 
  list(samp_name)]
fwrite(mw_01_names, file = mw_01_out, col.names = F)

mw_03_names <- samp_info[which(samp_info$subgrp_v2 == 'MW_03'), 
  list(samp_name)]
fwrite(mw_03_names, file = mw_03_out, col.names = F)

mw4X_intro_names <- samp_info[which(samp_info$MW_grp_1 == 'MW_4X_intro'),
  list(samp_name)]
fwrite(mw4X_intro_names, file = mw4X_intro_out, col.names = F)

quit(save = 'no')

