# Script for analyzing private alleles in the geobig sample set

#module load python/3.7-anaconda-2019.10
#source activate adegenet_2_env

### LOAD PACKAGES ###
library(adegenet)
library(parallel)
library(data.table)

allele_distr_func_file <- '/global/homes/g/grabowsp/tools/sg_8X/allele_distribution/allele_distribution_functions.r'
source(allele_distr_func_file)

### INPUT DATA ###
geno_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.rds'
genos <- readRDS(geno_file)

admix_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt'
admix_res <- fread(admix_res_file, header = F)

admix_k2_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.2.results.txt'
admix_k2_res <- fread(admix_k2_res_file, header = F)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

### SET OUTPUTS ###

out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig.50kSNPs.tetrasomic.CDS.allelepatterns.rds'

### SET VARIABLES ###

group_cutoff <- 0.9

###########

group_inds <- gen_admix_popind_list(admix_res, group_cutoff)

# Group 1
samp_meta[VCF_NAME %in% admix_res$V1[group_inds[['g1']]], .N, by=STATE]
# Atlantic coast
samp_meta[VCF_NAME %in% admix_res$V1[group_inds[['g1']]], .N, by=NQUIRE_PLOIDY]
# almost entirely 4X; 2 8X, 1 6X

# Group 2
samp_meta[VCF_NAME %in% admix_res$V1[group_inds[['g2']]], .N, by=STATE]
# TX/Gulf coast
samp_meta[VCF_NAME %in% admix_res$V1[group_inds[['g2']]], .N, by=NQUIRE_PLOIDY]
# Mainly 4X, 8 8X, 2 6X

# Group 3
samp_meta[VCF_NAME %in% admix_res$V1[group_inds[['g3']]], .N, by=STATE]
# Midwest, that's not very "midwest"
samp_meta[VCF_NAME %in% admix_res$V1[group_inds[['g3']]], .N, by=NQUIRE_PLOIDY]
# 130 4X, 134 8X

two_group_inds <- gen_admix_popind_list(admix_k2_res, group_cutoff)

# Group 1 
samp_meta[VCF_NAME %in% admix_res$V1[two_group_inds[['g1']]], .N, by=STATE]
# lowland
samp_meta[VCF_NAME %in% admix_res$V1[two_group_inds[['g1']]], .N, 
  by=NQUIRE_PLOIDY]
# mainly 4X

# Group 2
samp_meta[VCF_NAME %in% admix_res$V1[two_group_inds[['g2']]], .N, by=STATE]
# upland
samp_meta[VCF_NAME %in% admix_res$V1[two_group_inds[['g2']]], .N, 
  by=NQUIRE_PLOIDY]
# 50/50 4X/8X


three_group_list <- list()
four_group_list <- list()
tet_group_list <- list()
eco_group_list <- list()
eco_ploidy_list <- list()

for(gco in c(0.9, 0.95, 0.99)){
  co_name <- paste('cutoff_', gco, sep = '')
  group_inds <- gen_admix_popind_list(admix_res, gco)
  two_group_inds <- gen_admix_popind_list(admix_k2_res, gco)
  #
  # adjust so has informative names
  three_group_inds <- list()
  three_group_inds[['Atlantic']] <- group_inds[['g1']]
  three_group_inds[['Gulf']] <- group_inds[['g2']]
  three_group_inds[['Midwest']] <- group_inds[['g3']]
  three_group_list[[co_name]] <- three_group_inds
  ###
  # split midwest into 4X and 8X
  meta_8X_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '8X')]
  mw_8X_inds <- intersect(which(admix_res$V1 %in% meta_8X_names), 
    group_inds[['g3']])
  four_group_inds <- list()
  four_group_inds[['Atlantic']] <- group_inds[['g1']]
  four_group_inds[['Gulf']] <- group_inds[['g2']]
  four_group_inds[['MW4X']] <- setdiff(group_inds[['g3']], mw_8X_inds)
  four_group_inds[['MW8X']] <- mw_8X_inds
  four_group_list[[co_name]] <- four_group_inds
  ###
  # look only at 4X samples
  atl_8X_inds <- intersect(which(admix_res$V1 %in% meta_8X_names),
    group_inds[['g1']])
  glf_8X_inds <- intersect(which(admix_res$V1 %in% meta_8X_names),
    group_inds[['g2']])
  tet_group_inds <- list()
  tet_group_inds[['Atlantic4X']] <- setdiff(group_inds[['g1']], atl_8X_inds)
  tet_group_inds[['Gulf4X']] <- setdiff(group_inds[['g2']], glf_8X_inds)
  tet_group_inds[['MW4X']] <- setdiff(group_inds[['g3']], mw_8X_inds)
  tet_group_list[[co_name]] <- tet_group_inds
  # split into "ecotype"
  eco_inds <- list()
  eco_inds[['low']] <- two_group_inds[['g1']]
  eco_inds[['up']] <- two_group_inds[['g2']]
  eco_group_list[[co_name]] <- eco_inds
  # "ecotype" with upland split by ploidy
  up_8X_inds <- intersect(which(admix_k2_res$V1 %in% meta_8X_names),
    eco_inds[['up']])
  eco_ploidy_inds <- list()
  eco_ploidy_inds[['low']] <- eco_inds[['low']]
  eco_ploidy_inds[['up4X']] <- setdiff(eco_inds[['up']], up_8X_inds)
  eco_ploidy_inds[['up8X']] <- up_8X_inds
  eco_ploidy_list[[co_name]] <- eco_ploidy_inds
}

################

tet_group_allele_df <- get_shared_alleles_multicutoff(
  multicut_group_ind_list = tet_group_list, genos = genos, n_test_reps = 10)

three_group_allele_df <- get_shared_alleles_multicutoff(
  multicut_group_ind_list = three_group_list, genos = genos, n_test_reps = 10)

four_group_allele_df <- get_shared_alleles_multicutoff(
  multicut_group_ind_list = four_group_list, genos = genos, n_test_reps = 10)

####

eco_group_allele_df <- get_shared_alleles_multicutoff(
  multicut_group_ind_list = eco_group_list, genos = genos, n_test_reps = 10)

eco_ploidy_allele_df <- get_shared_alleles_multicutoff(
  multicut_group_ind_list = eco_ploidy_list, genos = genos, n_test_reps = 10)


tot_result_list <- list()
tot_result_list[['tet_samps']] <- tet_group_allele_df
tot_result_list[['three_groups']] <- three_group_allele_df
tot_result_list[['four_groups']] <- four_group_allele_df
tot_result_list[['eco_groups']] <- eco_group_allele_df
tot_result_list[['eco_ploidy']] <- eco_ploidy_allele_df

saveRDS(tot_result_list, file = out_file)

quit(save = 'no')

