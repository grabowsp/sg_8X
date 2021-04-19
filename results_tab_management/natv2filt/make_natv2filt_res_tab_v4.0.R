# Combine genepool-level splits with full results into combined table
### LOAD PACKAGES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

res_tab_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_res_tab_v1.0.txt'
res_tab <- fread(res_tab_file)

res_tab[ , sub_grp := as.character(NA)]

mw_res_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                    'natv2filt_MW_res_tab_v2.0.txt', sep = '')
mw_res_tab <- fread(mw_res_file)
# For MW: PC_1 - PC_5; mw_tot_k3

gulf_res_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                      'natv2filt_GULF_res_tab_v2.0.txt', sep = '')
gulf_res_tab <- fread(gulf_res_file)
# For GULF: PC_1 - PC_3; gulf_tot_k3

atl_res_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                     'natv2filt_ATL_res_tab_v2.0.txt', sep = '')
atl_res_tab <- fread(atl_res_file)
# For ATL: PC_1 - PC_4; atl_tot_k4

mf_res_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                    'natv2filt_MISFIT_res_tab_v2.0.txt', sep = '')
mf_res_tab <- fread(mf_res_file)
# For MISFITs: PC_1 - PC_3

### SET OUTPUT ###

res_out_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_res_tab_v1.0.txt'

### FUNCTIONS FOR ANALYSIS ######
standardize_PC <- function(pc_vals){
  tmp_vals <- pc_vals - min(pc_vals, na.rm = T)
  tmp_vals <- tmp_vals/max(tmp_vals, na.rm = T)
  return(tmp_vals)
}

########################

# Midwest

sum(res_tab$samp_name == mw_res_tab$samp_name)
# [1] 706

# For MW:
# PC_1 - PC_5
# mw_tot_k3
res_tab[ , MW_pc01_raw := mw_res_tab$mw_tot_PC1]
res_tab[ , MW_pc01_stand := standardize_PC(mw_res_tab$mw_tot_PC1)]

res_tab[ , MW_pc02_raw := mw_res_tab$mw_tot_PC2]
res_tab[ , MW_pc02_stand := standardize_PC(mw_res_tab$mw_tot_PC2)]

res_tab[ , MW_pc03_raw := mw_res_tab$mw_tot_PC3]
res_tab[ , MW_pc03_stand := standardize_PC(mw_res_tab$mw_tot_PC3)]

res_tab[ , MW_pc04_raw := mw_res_tab$mw_tot_PC4]
res_tab[ , MW_pc04_stand := standardize_PC(mw_res_tab$mw_tot_PC4)]

res_tab[ , MW_pc05_raw := mw_res_tab$mw_tot_PC5]
res_tab[ , MW_pc05_stand := standardize_PC(mw_res_tab$mw_tot_PC5)]

res_tab[ , MW_admix_k3_MW8X_West := mw_res_tab$mw_tot_k3_pop_1]
res_tab[ , MW_admix_k3_MW4X := mw_res_tab$mw_tot_k3_pop_2]
res_tab[ , MW_admix_k3_MW8X_East := mw_res_tab$mw_tot_k3_pop_3]

for(i in seq(nrow(mw_res_tab))){
  if(is.na(mw_res_tab$grps[i]) == F)
  res_tab[i, sub_grp := mw_res_tab$grps[i]]
}

#############

# Gulf
sum(res_tab$samp_name == gulf_res_tab$samp_name)
# [1] 706

# For GULF:
# PC_1 - PC_3
# gulf_tot_k3

res_tab[ , GULF_pc01_raw := gulf_res_tab$gulf_tot_PC1]
res_tab[ , GULF_pc01_stand := standardize_PC(gulf_res_tab$gulf_tot_PC1)]

res_tab[ , GULF_pc02_raw := gulf_res_tab$gulf_tot_PC2]
res_tab[ , GULF_pc02_stand := standardize_PC(gulf_res_tab$gulf_tot_PC2)]

res_tab[ , GULF_pc03_raw := gulf_res_tab$gulf_tot_PC3]
res_tab[ , GULF_pc03_stand := standardize_PC(gulf_res_tab$gulf_tot_PC3)]

res_tab[ , GULF_admix_k3_GULF_COAST_1 := gulf_res_tab$gulf_tot_k3_pop_1]
res_tab[ , GULF_admix_k3_GULF_COAST_2 := gulf_res_tab$gulf_tot_k3_pop_2]
res_tab[ , GULF_admix_k3_TX_MX := gulf_res_tab$gulf_tot_k3_pop_3]

for(i in seq(nrow(gulf_res_tab))){
  if(is.na(gulf_res_tab$grps[i]) == F)
    res_tab[i, sub_grp := gulf_res_tab$grps[i]]
}

##############

# Atlantic

sum(res_tab$samp_name == atl_res_tab$samp_name)
# 706

# For ATL:
# PC_1 - PC_4
# atl_tot_k4

res_tab[ , ATL_pc01_raw := atl_res_tab$atl_tot_PC1]
res_tab[ , ATL_pc01_stand := standardize_PC(atl_res_tab$atl_tot_PC1)]

res_tab[ , ATL_pc02_raw := atl_res_tab$atl_tot_PC2]
res_tab[ , ATL_pc02_stand := standardize_PC(atl_res_tab$atl_tot_PC2)]

res_tab[ , ATL_pc03_raw := atl_res_tab$atl_tot_PC3]
res_tab[ , ATL_pc03_stand := standardize_PC(atl_res_tab$atl_tot_PC3)]

res_tab[ , ATL_pc04_raw := atl_res_tab$atl_tot_PC4]
res_tab[ , ATL_pc04_stand := standardize_PC(atl_res_tab$atl_tot_PC4)]

res_tab[ , ATL_admix_k4_MID_ATL:= atl_res_tab$atl_tot_k4_pop_1]
# Mid atlantic
res_tab[ , ATL_admix_k4_NE:= atl_res_tab$atl_tot_k4_pop_2]
# New England
res_tab[ , ATL_admix_k4_SE:= atl_res_tab$atl_tot_k4_pop_3]
# South East Coast
res_tab[ , ATL_admix_k4_CH_BAY:= atl_res_tab$atl_tot_k4_pop_4]
# Chesapeake Bay

for(i in seq(nrow(atl_res_tab))){
  if(is.na(atl_res_tab$grps[i]) == F)
    res_tab[i, sub_grp := atl_res_tab$grps[i]]
}

###########

# Misfit results

# For MISFITs:
# PC_1 - PC_3

mf_inds <- c()
for(i in mf_res_tab$samp_name){
  tmp_ind <- which(res_tab$samp_name == i)
  mf_inds <- c(mf_inds, tmp_ind)
}

res_tab[ mf_inds, MF_pc01_raw := mf_res_tab$PC1]
res_tab[ mf_inds, MF_pc01_stand := standardize_PC(mf_res_tab$PC1)]

res_tab[ mf_inds, MF_pc02_raw := mf_res_tab$PC2]
res_tab[ mf_inds, MF_pc02_stand := standardize_PC(mf_res_tab$PC2)]
  
res_tab[ mf_inds, MF_pc03_raw := mf_res_tab$PC3]
res_tab[ mf_inds, MF_pc03_stand := standardize_PC(mf_res_tab$PC3)]  

res_tab[mf_inds, sub_grp := mf_res_tab$grps]  

fwrite(res_tab, file = res_out_file, sep = '\t')


