# Script for analyzing private alleles in the geobig sample set with major
#  genepools further divided

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
# geno_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/GW.500kSNPs.tetrasomic.CDS.geobig.genlight.rds'
genos <- readRDS(geno_file)

# geobig admix results
gb_admix_k3_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt'
gb_admix_k3_res <- fread(gb_admix_k3_res_file, header = F)

gb_admix_k2_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.2.results.txt'
gb_admix_k2_res <- fread(gb_admix_k2_res_file, header = F)

# NorthInland admix results
ni_admix_k3_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_NorthInland_admix/GW_50k_geobigNorthInland.3.results.txt'
ni_admix_k3_res <- fread(ni_admix_k3_res_file)

# SouthCoastal admix results
sc_admix_k2_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.2.results.txt'
sc_admix_k2_res <- fread(sc_admix_k2_res_file)

sc_admix_k4_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.4.results.txt'
sc_admix_k4_res <- fread(sc_admix_k4_res_file)


meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

### SET OUTPUTS ###

out_file_1 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_3NI_vs_2SC.50kSNPs.tetrasomic.CDS.allelepatterns.rds'

### SET VARIABLES ###

group_cutoff <- 0.9

###########

#gb_k2_group_inds <- gen_admix_popind_list(gb_admix_k2_res, group_cutoff)
#gb_k3_group_inds <- gen_admix_popind_list(gb_admix_k3_res, group_cutoff)

#ni_k3_group_inds <- gen_admix_popind_list(ni_admix_k3_res, group_cutoff)

#sc_k2_group_inds <- gen_admix_popind_list(sc_admix_k2_res, group_cutoff)
#sc_k4_group_inds <- gen_admix_popind_list(sc_admix_k4_res, group_cutoff)

gb_k2_group_names <- gen_admix_pop_samp_list(gb_admix_k2_res, group_cutoff)
gb_k3_group_names <- gen_admix_pop_samp_list(gb_admix_k3_res, group_cutoff)

ni_k3_group_names <- gen_admix_pop_samp_list(ni_admix_k3_res, group_cutoff)

sc_k2_group_names <- gen_admix_pop_samp_list(sc_admix_k2_res, group_cutoff)
sc_k4_group_names <- gen_admix_pop_samp_list(sc_admix_k4_res, group_cutoff)

### gb_K2
## Group 1
samp_meta[VCF_NAME %in% gb_k2_group_names[['g1']], .N, by=STATE]
# Lowland/SouthCoastal
samp_meta[VCF_NAME %in% gb_k2_group_names[['g1']], .N, by=NQUIRE_PLOIDY]
# almost entirely 4X; 14 8X, 3 6X
## Group 2
samp_meta[VCF_NAME %in% gb_k2_group_names[['g2']], .N, by=STATE]
# upland/NorthInland
samp_meta[VCF_NAME %in% gb_k2_group_names[['g2']], .N, by=NQUIRE_PLOIDY]
# 50:50 4X:8X

### gb_K3
## Group 1
samp_meta[VCF_NAME %in% gb_k3_group_names[['g1']], .N, by=STATE]
# Atlantic
samp_meta[VCF_NAME %in% gb_k3_group_names[['g1']], .N, by=NQUIRE_PLOIDY]
# Mainly 4X
## Group 2
samp_meta[VCF_NAME %in% gb_k3_group_names[['g2']], .N, by=STATE]
# Gulf Coast/Tx
samp_meta[VCF_NAME %in% gb_k3_group_names[['g2']], .N, by=NQUIRE_PLOIDY]
# Mainly 4X, 8 8X, 2 6X
## Group 3
samp_meta[VCF_NAME %in% gb_k3_group_names[['g3']], .N, by=STATE]
# Midwest
samp_meta[VCF_NAME %in% gb_k3_group_names[['g3']], .N, by=NQUIRE_PLOIDY]
# 50:50 4X:8X

### ni_K3
## Group 1
samp_meta[VCF_NAME %in% ni_k3_group_names[['g1']], .N, by=STATE]
# Cosmopolitan 8X
samp_meta[VCF_NAME %in% ni_k3_group_names[['g1']], .N, by=NQUIRE_PLOIDY]
# Only 8X
## Group 2
samp_meta[VCF_NAME %in% ni_k3_group_names[['g2']], .N, by=STATE]
# Eastern 8X
samp_meta[VCF_NAME %in% ni_k3_group_names[['g2']], .N, by=NQUIRE_PLOIDY]
# Only 8X
## Group 3
samp_meta[VCF_NAME %in% ni_k3_group_names[['g3']], .N, by=STATE]
# Midwest 4X
samp_meta[VCF_NAME %in% ni_k3_group_names[['g3']], .N, by=NQUIRE_PLOIDY]
# Only 4X

### sc_K2
## Group 1
samp_meta[VCF_NAME %in% sc_k2_group_names[['g1']], .N, by=STATE]
# Texas
samp_meta[VCF_NAME %in% sc_k2_group_names[['g1']], .N, by=NQUIRE_PLOIDY]
# Almost only 4X, 1 8X, 2 6X
## Group 2
samp_meta[VCF_NAME %in% sc_k2_group_names[['g2']], .N, by=STATE]
# Atlantic
samp_meta[VCF_NAME %in% sc_k2_group_names[['g2']], .N, by=NQUIRE_PLOIDY]
# Almost only 4X, 2 8X, 1 6X

### sc_K4
## Group 1
samp_meta[VCF_NAME %in% sc_k4_group_names[['g1']], .N, by=STATE]
# Texas
samp_meta[VCF_NAME %in% sc_k4_group_names[['g1']], .N, by=NQUIRE_PLOIDY]
# Almost only 4X, 1 8X, 1 6X
## Group 2
samp_meta[VCF_NAME %in% sc_k4_group_names[['g2']], .N, by=STATE]
# Mississippi
samp_meta[VCF_NAME %in% sc_k4_group_names[['g2']], .N, by=NQUIRE_PLOIDY]
# Exclusively 4X
## Group 3
samp_meta[VCF_NAME %in% sc_k4_group_names[['g3']], .N, by=STATE]
# Northern Atlantic
samp_meta[VCF_NAME %in% sc_k4_group_names[['g3']], .N, by=NQUIRE_PLOIDY]
# Exclusively 4X
## Group 4
samp_meta[VCF_NAME %in% sc_k4_group_names[['g4']], .N, by=STATE]
# Southern Atlantic
samp_meta[VCF_NAME %in% sc_k4_group_names[['g4']], .N, by=NQUIRE_PLOIDY]
# Almost exclusively 4X, 2 8X

gb_k2_groups_in_genos <- lapply(gb_k2_group_names, function(x)
  which(indNames(genos) %in% x))
gb_k3_groups_in_genos <- lapply(gb_k3_group_names, function(x)
  which(indNames(genos) %in% x))



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

