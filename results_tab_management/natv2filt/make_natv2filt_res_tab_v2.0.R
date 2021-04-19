### LOAD MODULES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###
samp_meta_file <- paste('/Users/grabowsk/Analysis/sg_8X/metadata_management/',
  'meta_files/', 'sg_8X_metadata_v3.1.csv', sep = '')
samp_meta <- fread(samp_meta_file)

res_tab_file <- paste('/Users/grabowsk/data/sg_8X_analysis/',
  'natv2filt_res_tab_v1.0.txt', sep = '')
res_tab <- fread(res_tab_file)

### PCA results ###
pca_res_dir <- '/Users/grabowsk/data/sg_8X_analysis/'
mw_tot_pca_res_file <- paste(pca_res_dir, 
  'GW.100kSNPs.tetrasomic.CDS.natv2filt.MW_tot.genlight.PCAresults.rds',
  sep = '')
mw_tot_pca_res <- readRDS(mw_tot_pca_res_file)

gulf_tot_pca_res_file <- paste(pca_res_dir, 
  'GW.100kSNPs.tetrasomic.CDS.natv2filt.GULF_tot.genlight.PCAresults.rds',
  sep = '')
gulf_tot_pca_res <- readRDS(gulf_tot_pca_res_file)

atl_tot_pca_res_file <- paste(pca_res_dir,
  'GW.100kSNPs.tetrasomic.CDS.natv2filt.ATL_tot.genlight.PCAresults.rds',
  sep = '')
atl_tot_pca_res <- readRDS(atl_tot_pca_res_file)

misfit_pca_res_file <- paste(pca_res_dir, 
  'GW.100kSNPs.tetrasomic.CDS.natv2filt.MISFIT.genlight.PCAresults.rds',
  sep = '')
misfit_pca_res <- readRDS(misfit_pca_res_file)

### ADMIXTURE results ###
admix_res_dir <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_gp_admix_res/'

atl_main_admix_k3_file <- paste(admix_res_dir, 
                                'GW_natv2filt_ATL_main.3.results.txt', 
                                sep = '')
atl_main_admix_k3_res <- fread(atl_main_admix_k3_file)
atl_tot_admix_k3_file <- paste(admix_res_dir, 
                                'GW_natv2filt_ATL_tot.3.results.txt', 
                                sep = '')
atl_tot_admix_k3_res <- fread(atl_tot_admix_k3_file)
atl_tot_admix_k4_file <- paste(admix_res_dir, 
                               'GW_natv2filt_ATL_tot.4.results.txt', 
                               sep = '')
atl_tot_admix_k4_res <- fread(atl_tot_admix_k4_file)

gulf_main_admix_k2_file <- paste(admix_res_dir, 
                                 'GW_natv2filt_GULF_main.2.results.txt', 
                                 sep = '')
gulf_main_admix_k2_res <- fread(gulf_main_admix_k2_file)
gulf_tot_admix_k3_file <- paste(admix_res_dir, 
                                 'GW_natv2filt_GULF_tot.3.results.txt', 
                                 sep = '')
gulf_tot_admix_k3_res <- fread(gulf_tot_admix_k3_file)
gulf_tot_admix_k4_file <- paste(admix_res_dir, 
                                'GW_natv2filt_GULF_tot.4.results.txt', 
                                sep = '')
gulf_tot_admix_k4_res <- fread(gulf_tot_admix_k4_file)

mw_main_admix_k2_file <- paste(admix_res_dir, 
                               'GW_natv2filt_mw_main.2.results.txt', 
                               sep = '')
mw_main_admix_k2_res <- fread(mw_main_admix_k2_file)
mw_main_admix_k3_file <- paste(admix_res_dir, 
                               'GW_natv2filt_mw_main.3.results.txt', 
                               sep = '')
mw_main_admix_k3_res <- fread(mw_main_admix_k3_file)
mw_tot_admix_k2_file <- paste(admix_res_dir, 
                               'GW_natv2filt_mw_tot.2.results.txt', 
                               sep = '')
mw_tot_admix_k2_res <- fread(mw_tot_admix_k2_file)
mw_tot_admix_k3_file <- paste(admix_res_dir, 
                              'GW_natv2filt_mw_tot.3.results.txt', 
                              sep = '')
mw_tot_admix_k3_res <- fread(mw_tot_admix_k3_file)

### SET OUTPUTS ###
mw_tot_res_tab_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                            'natv2filt_MW_res_tab_v1.0.txt', sep = '')
gulf_tot_res_tab_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                            'natv2filt_GULF_res_tab_v1.0.txt', sep = '')
atl_tot_res_tab_out <- paste('/Users/grabowsk/data/sg_8X_analysis/',
                              'natv2filt_ATL_res_tab_v1.0.txt', sep = '')

### FUNCTIONS ###
prep_pca_results <- function(pca_res, cn_pre, max_pc){
  # Prepare pca results for incorporating into larger table
  # INPUTS
  # pca_res = PCA results list
  # cn_pre = the prefix for the colnames for the outputted matrix
  # max_pc = the highest PC to be included
  #######
  tmp_mat <- pca_res$scores[, c(1:max_pc)]
  colnames(tmp_mat) <- paste(cn_pre, colnames(tmp_mat), sep = '')
  return(tmp_mat)
}

prep_admix_results <- function(admix_res, cn_pre){
  # Prepare Admix results for incorporating into larger table
  # admix_res = ADMIXTURE results
  # cn_pre = prefix for the colnames for the results columns
  #####
  tmp_tab <- admix_res
  colnames(tmp_tab) <- c('samp_name', paste(cn_pre, 'pop_', 
                                            seq(ncol(admix_res)-1), 
                                            sep = ''))
  return(tmp_tab)
}

start_res_tab <- function(tot_res_tab, pca_res_mat){
  test_tot_tab <- tot_res_tab[, c('samp_name', 'ploidy')]
  test_tot_tab[, colnames(pca_res_mat):= as.numeric(NA)]
  tmp_tot_inds <- c()
  for(i in rownames(pca_res_mat)){
    tmp_ind <- which(test_tot_tab$samp_name == i)
    tmp_tot_inds <- c(tmp_tot_inds, tmp_ind)
  }
  for(j in seq(ncol(pca_res_mat))){
    test_tot_tab[tmp_tot_inds, colnames(pca_res_mat)[j] := pca_res_mat[,j]]
  }
  return(test_tot_tab)
}

add_admix_results <- function(res_tab, admix_res_tab){
  res_tab[, colnames(admix_res_tab)[-1] <- as.numeric(NA)]
  tmp_tot_inds <- c()
  for(i in admix_res_tab$samp_name){
    tmp_ind <- which(res_tab$samp_name == i)
    tmp_tot_inds <- c(tmp_tot_inds, tmp_ind)
  }
  for(j in c(2:ncol(admix_res_tab))){
    res_tab[tmp_tot_inds, 
            colnames(admix_res_tab)[j] := 
              admix_res_tab[ ,colnames(admix_res_tab)[j], with = F]]
  }
  return(res_tab)
}

#########

# Make ATL results
## ATL pca sub mat
## max pc = 4
atl_pca_mat <- prep_pca_results(pca_res = atl_tot_pca_res, cn_pre = 'atl_tot_',
                                max_pc = 4)
## ATL admix results
atl_main_k3_tab <- prep_admix_results(admix_res = atl_main_admix_k3_res,
                                      cn_pre = 'atl_main_k3_')
atl_tot_k3_tab <- prep_admix_results(admix_res = atl_tot_admix_k3_res,
                                     cn_pre = 'atl_tot_k3_')
atl_tot_k4_tab <- prep_admix_results(admix_res = atl_tot_admix_k4_res,
                                     cn_pre = 'atl_tot_k4_')
# make total ATLANTIC results table
atl_res_tab <- start_res_tab(tot_res_tab = res_tab, pca_res_mat = atl_pca_mat)
atl_res_tab <- add_admix_results(res_tab = atl_res_tab,
                                 admix_res_tab = atl_main_k3_tab)
atl_res_tab <- add_admix_results(res_tab = atl_res_tab,
                                 admix_res_tab = atl_tot_k3_tab)
atl_res_tab <- add_admix_results(res_tab = atl_res_tab,
                                 admix_res_tab = atl_tot_k4_tab)

# Make GULF results
## Gulf pca sub mat
## max pc = 3
gulf_pca_mat <- prep_pca_results(pca_res = gulf_tot_pca_res, 
                                 cn_pre = 'gulf_tot_', max_pc = 3)
## Gulf admix results
gulf_main_k2_tab <- prep_admix_results(admix_res = gulf_main_admix_k2_res,
                                       cn_pre = 'gulf_main_k2_')
gulf_tot_k3_tab <- prep_admix_results(admix_res = gulf_tot_admix_k3_res,
                                       cn_pre = 'gulf_tot_k3_')
gulf_tot_k4_tab <- prep_admix_results(admix_res = gulf_tot_admix_k4_res,
                                      cn_pre = 'gulf_tot_k4_')
# make total GULF results table
gulf_res_tab <- start_res_tab(tot_res_tab = res_tab, 
                              pca_res_mat = gulf_pca_mat)
gulf_res_tab <- add_admix_results(res_tab = gulf_res_tab,
                                  admix_res_tab = gulf_main_k2_tab)
gulf_res_tab <- add_admix_results(res_tab = gulf_res_tab,
                                  admix_res_tab = gulf_tot_k3_tab)
gulf_res_tab <- add_admix_results(res_tab = gulf_res_tab,
                                  admix_res_tab = gulf_tot_k4_tab)

# Make MW results
# pca pre-mat, max pc = 5
mw_pca_mat <- prep_pca_results(pca_res = mw_tot_pca_res, 
                               cn_pre = 'mw_tot_', max_pc = 5)
# MW admix results
mw_main_k2_tab <- prep_admix_results(admix_res = mw_main_admix_k2_res,
                                     cn_pre = 'mw_main_k2_')
mw_main_k3_tab <- prep_admix_results(admix_res = mw_main_admix_k3_res,
                                     cn_pre = 'mw_main_k3_')
mw_tot_k2_tab <- prep_admix_results(admix_res = mw_tot_admix_k2_res,
                                     cn_pre = 'mw_tot_k2_')
mw_tot_k3_tab <- prep_admix_results(admix_res = mw_tot_admix_k3_res,
                                     cn_pre = 'mw_tot_k3_')
# make total MW table
mw_res_tab <- start_res_tab(tot_res_tab = res_tab, 
                              pca_res_mat = mw_pca_mat)
mw_res_tab <- add_admix_results(res_tab = mw_res_tab,
                                  admix_res_tab = mw_main_k2_tab)
mw_res_tab <- add_admix_results(res_tab = mw_res_tab,
                                admix_res_tab = mw_main_k3_tab)
mw_res_tab <- add_admix_results(res_tab = mw_res_tab,
                                admix_res_tab = mw_tot_k2_tab)
mw_res_tab <- add_admix_results(res_tab = mw_res_tab,
                                admix_res_tab = mw_tot_k3_tab)

# save results
fwrite(mw_res_tab, file = mw_tot_res_tab_out, sep = '\t')
fwrite(gulf_res_tab, file = gulf_tot_res_tab_out, sep = '\t')
fwrite(atl_res_tab, file = atl_tot_res_tab_out, sep = '\t')
