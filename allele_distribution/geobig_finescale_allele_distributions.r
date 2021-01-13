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
#geno_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.rds'
geno_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/GW.500kSNPs.tetrasomic.CDS.geobig.genlight.rds'
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

# with 50k SNPs
#out_file_1 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_manycomps.50kSNPs.tetrasomic.CDS.allelepatterns.rds'
#out_file_2 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_manycomps.50kSNPs.tetrasomic.CDS.allelepatterns_notcombined.rds'

# with 500k SNPs
#out_file_1 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_manycomps.500kSNPs.tetrasomic.CDS.allelepatterns.rds'
#out_file_2 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_manycomps.500kSNPs.tetrasomic.CDS.allelepatterns_notcombined.rds'

out_file_1 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_manycomps.500kSNPs.tetrasomic.CDS.allelepatterns_7groups.rds'
out_file_2 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/geobig_manycomps.500kSNPs.tetrasomic.CDS.allelepatterns_notcombined_7groups.rds'


### SET VARIABLES ###

group_cutoff <- 0.9

num_test_reps <- 25

###########

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

# SET CONTRASTS:

# 1) 4X, 3 groups: MW-4X, Gulf-4X, Atlantic_4X
# 2) 3 groups all ploidy: MW, Gulf, Atlantic
# 3) 4 groups, NorIn split by ploidy: MW-4X, NorInd-8X, Gulf, Atlantic
# 4) 5 groups, NorIn split into 3: MW-4X, Cosmo-8X, ME-8X; Gulf, Atlantic
# 5) 6 groups, SoCo split into 4, NorIn split by ploidy: MW-4X, NorInd-8X, 
#        Gulf-TX, Gulf-MS, SoAtlantic, NoAtlantic
# 6) 7 groups, NorIn split into 3, SoCo split into 4: MW-4X, Cosmo-8X, ME-8X,
#        Gulf-TX, Gulf-MS, SoAtlantic, NoAtlantic 

# 3 groups, only 4X

three_group_4X_list <- list()
three_group_4X_list[['Atlantic_4X']] <- setdiff(gb_k3_group_names[['g1']], 
  samp_meta[NQUIRE_PLOIDY == '8X' | NQUIRE_PLOIDY == '6X', VCF_NAME])
three_group_4X_list[['Gulf_4X']] <- setdiff(gb_k3_group_names[['g2']], 
  samp_meta[NQUIRE_PLOIDY == '8X' | NQUIRE_PLOIDY == '6X', VCF_NAME])
three_group_4X_list[['MW_4X']] <- setdiff(gb_k3_group_names[['g3']],
  samp_meta[NQUIRE_PLOIDY == '8X' | NQUIRE_PLOIDY == '6X', VCF_NAME])

three_group_4X_allele_many <- gen_share_unique_useNames_list(
  group_name_list = three_group_4X_list, genos = genos, 
  n_test_reps = num_test_reps)
three_group_4X_allele_res <- process_share_results(three_group_4X_allele_many)

# 3 groups, all ploidy

three_group_all_list <- list()
three_group_all_list[['Atlantic_all']] <- gb_k3_group_names[['g1']] 
three_group_all_list[['Gulf_all']] <- gb_k3_group_names[['g2']]
three_group_all_list[['MW_all']] <- gb_k3_group_names[['g3']]

three_group_all_allele_many <- gen_share_unique_useNames_list(
  group_name_list = three_group_all_list, genos = genos, 
  n_test_reps = num_test_reps)
three_group_all_allele_res <- process_share_results(three_group_all_allele_many)

# 4 groups, NorIn split by ploidy
four_group_list <- list()
four_group_list[['Atlantic']] <- gb_k3_group_names[['g1']]
four_group_list[['Gulf']] <- gb_k3_group_names[['g2']]
four_group_list[['MW-4X']] <- intersect(gb_k3_group_names[['g3']], 
  samp_meta[NQUIRE_PLOIDY == '4X', VCF_NAME])
four_group_list[['NorthInland-8X']] <- intersect(gb_k3_group_names[['g3']],
  samp_meta[NQUIRE_PLOIDY == '8X' | NQUIRE_PLOIDY == '6X', VCF_NAME])

four_group_allele_many <- gen_share_unique_useNames_list(
  group_name_list = four_group_list, genos = genos, 
  n_test_reps = num_test_reps)
four_group_allele_res <- process_share_results(four_group_allele_many)

# 5 groups, NorIn split into 3
five_group_list <- list()
five_group_list[['Atlantic']] <- gb_k3_group_names[['g1']]
five_group_list[['Gulf']] <- gb_k3_group_names[['g2']]
five_group_list[['MW-4X']] <- ni_k3_group_names[['g3']]
five_group_list[['Cosmo-8X']] <- ni_k3_group_names[['g1']]
five_group_list[['ME-8X']] <- ni_k3_group_names[['g2']]

five_group_allele_many <- gen_share_unique_useNames_list(
  group_name_list = five_group_list, genos = genos, 
  n_test_reps = num_test_reps)
five_group_allele_res <- process_share_results(five_group_allele_many)

# 6 groups, SoCo split into 4, NorIn split into 2
six_group_list <- list()
six_group_list[['MW-4X']] <- intersect(gb_k3_group_names[['g3']],
  samp_meta[NQUIRE_PLOIDY == '4X', VCF_NAME])
six_group_list[['NorthInland-8X']] <- intersect(gb_k3_group_names[['g3']],
  samp_meta[NQUIRE_PLOIDY == '8X' | NQUIRE_PLOIDY == '6X', VCF_NAME])
six_group_list[['Gulf-TX']] <- sc_k4_group_names[['g1']]
six_group_list[['Gulf-MS']] <- sc_k4_group_names[['g2']]
six_group_list[['Atlantic-South']] <- sc_k4_group_names[['g4']]
six_group_list[['Atlantic-North']] <- sc_k4_group_names[['g3']]

six_group_allele_many <- gen_share_unique_useNames_list(
  group_name_list = six_group_list, genos = genos, 
  n_test_reps = num_test_reps)
six_group_allele_res <- process_share_results(six_group_allele_many)

# 7 groups, NorIn split into 3, SoCo split into 4
seven_group_list <- list()
seven_group_list[['MW-4X']] <- ni_k3_group_names[['g3']]
seven_group_list[['Cosmo-8X']] <- ni_k3_group_names[['g1']]
seven_group_list[['ME-8X']] <- ni_k3_group_names[['g2']]
seven_group_list[['Gulf-TX']] <- sc_k4_group_names[['g1']]
seven_group_list[['Gulf-MS']] <- sc_k4_group_names[['g2']]
seven_group_list[['Atlantic-South']] <- sc_k4_group_names[['g4']]
seven_group_list[['Atlantic-North']] <- sc_k4_group_names[['g3']]

seven_group_allele_many <- gen_share_unique_useNames_list(
  group_name_list = seven_group_list, genos = genos, 
  n_test_reps = num_test_reps)
seven_group_allele_res <- process_share_results(seven_group_allele_many)

many_res_list <- list()
combo_res_list <- list()

many_res_list[['1North_2South_4X']] <- three_group_4X_allele_many
combo_res_list[['1North_2South_4X']] <- three_group_4X_allele_res

many_res_list[['1North_v_2South_all']] <- three_group_all_allele_many
combo_res_list[['1North_v_2South_all']] <- three_group_all_allele_res

many_res_list[['North4X_North8X_2South']] <- four_group_allele_many
combo_res_list[['North4X_North8X_2South']] <- four_group_allele_res

many_res_list[['3North_2South']] <- five_group_allele_many
combo_res_list[['3North_2South']] <- five_group_allele_res

many_res_list[['North4X_North8X_4South']] <- six_group_allele_many
combo_res_list[['North4X_North8X_4South']] <- six_group_allele_res

many_res_list[['3North_4South']] <- seven_group_allele_many
combo_res_list[['3North_4South']] <- seven_group_allele_res

saveRDS(combo_res_list, file = out_file_1)
saveRDS(many_res_list, file = out_file_2)

quit(save = 'no')

