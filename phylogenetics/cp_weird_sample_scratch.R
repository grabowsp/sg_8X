# scratch space for analysis of weird samples while network is down

### LOAD LIBRARIES ###
library(data.table)

### INPUT DATA ###
res_tab_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_res_tab_v4.0.txt'
res_tab <- fread(res_tab_file)

meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

weird_samp_file <- '/Users/grabowsk/Downloads/List_of_upland_cp_LIBs.csv'
weird_libs <- fread(weird_samp_file)

######

weird_samps <- samp_meta[which(samp_meta$LIB %in% weird_libs$LIB), VCF_NAME]
res_weird_inds <- which(res_tab$samp_name %in% weird_samps)

weird_gulf_inds <- intersect(res_weird_inds, grep('GULF_', res_tab$sub_grp))
weird_gulf_ord <- order(res_tab[weird_gulf_inds, samp_name])
res_tab[weird_gulf_inds[weird_gulf_ord], list(samp_name, ploidy)]


              