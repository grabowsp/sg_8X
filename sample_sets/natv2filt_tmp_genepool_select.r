# make "temp" genepool selections from natv2filt PCA
#  cutoffs chosen in selecting_natv2filt_subgroups.R

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate /global/homes/g/grabowsp/.conda/envs/R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_tab_file <- '/global/homes/g/grabowsp/data/switchgrass/results_tables_8X/natv2filt_res_tab_v1.0.txt'
res_tab <- fread(res_tab_file)

### SET OUTPUTS ###
out_dir <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/'

mw_main_out <- paste(out_dir, 'natv2filt_tmp_MW_main_names.txt', sep = '')
mw_tot_out <- paste(out_dir, 'natv2filt_tmp_MW_total_names.txt', sep = '')

gulf_main_out <- paste(out_dir, 'natv2filt_tmp_GULF_main_names.txt', sep = '')
gulf_tot_out <- paste(out_dir, 'natv2filt_tmp_GULF_total_names.txt', sep = '')

atl_main_out <- paste(out_dir, 'natv2filt_tmp_ATL_main_names.txt', sep = '')
atl_tot_out <- paste(out_dir, 'natv2filt_tmp_ATL_total_names.txt', sep = '')

misfit_out <- paste(out_dir, 'natv2filt_tmp_MISFIT_names.txt', sep = '')

### SET VARIABLES ###
mw_cut <- 0.7
gulf_cut <- 0.7
atl_cut <- 0.65

###########

mw_main_inds <- which((res_tab$full_pc01_raw > 22 & 
                         res_tab$full_pc02_raw > 2.5) | 
                      (res_tab$full_pc01_raw > 26 & 
                         res_tab$full_pc02_raw > 1))
mw_main_names <- res_tab$samp_name[mw_main_inds]

mw_tot_inds <- setdiff(which(res_tab$full_admix_k3_MW > mw_cut), 
                 which(res_tab$samp_name == 'J581.A'))
mw_tot_names <- res_tab$samp_name[mw_tot_inds]

###

gulf_main_inds <- which(res_tab$full_pc01_raw < -3 & 
                         res_tab$full_pc02_raw < -20)
gulf_main_names <- res_tab$samp_name[gulf_main_inds]

gulf_tot_inds <- which(res_tab$full_admix_k3_GULF > gulf_cut)
gulf_tot_names <- res_tab$samp_name[gulf_tot_inds]

###

atl_main_inds <- which(res_tab$full_pc01_raw < -15 & 
                         res_tab$full_pc02_raw > 4)
atl_main_names <- res_tab$samp_name[atl_main_inds]

atl_tot_inds <- setdiff(which(res_tab$full_admix_k3_ATL > atl_cut),
                         c(which(res_tab$samp_name == 'J534.A'), 
                           which(res_tab$samp_name == 'J530.B')))
atl_tot_names <- res_tab$samp_name[atl_tot_inds]

###

misfit_inds <- union(which(res_tab$full_admix_k3_MW < mw_cut & 
                            res_tab$full_admix_k3_GULF < gulf_cut &
                            res_tab$full_admix_k3_ATL < atl_cut),
                     which(res_tab$samp_name %in% 
                            c('J581.A', 'J534.A', 'J530.B')))
misfit_names <- res_tab$samp_name[misfit_inds]

# check that all samples are included
setdiff(res_tab$samp_name, c(mw_tot_names, gulf_tot_names, atl_tot_names, 
  misfit_names))
# character(0)

### Write name files

write.table(mw_main_names, file = mw_main_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)
write.table(mw_tot_names, file = mw_tot_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

write.table(gulf_main_names, file = gulf_main_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)
write.table(gulf_tot_names, file = gulf_tot_out, quote = F, sep = '\t',  
  row.names = F, col.names = F)

write.table(atl_main_names, file = atl_main_out, quote = F, sep = '\t',
  row.names = F, col.names = F)
write.table(atl_tot_names, file = atl_tot_out, quote = F, sep = '\t',
  row.names = F, col.names = F)

write.table(misfit_names, file = misfit_out, quote = F, sep = '\t',
  row.names = F, col.names = F)

quit(save = 'no')
