# Script for dividing the geo_big sample set into 2 groups
#  "N_Inland" = "North/Inland"; previous "upland" clade, but can't really use 
#                 "upland because that's polyphyletic
#  "S_Coastal" = "SouthCoastal"; previous "lowland" clade, encompases Atlantic 
#                  and Gulf Groups

# module load python/3.7-anaconda-2019.10
# source activate adegenet_2_env

### LOAD PACKAGES #####
library(adegenet)
library(data.table)

allele_distr_func_file <- '/global/homes/g/grabowsp/tools/sg_8X/allele_distribution/allele_distribution_functions.r'
source(allele_distr_func_file)

### INPUT DATA ###

k2_admix_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.2.results.txt'
k2_admix_res <- fread(k2_admix_res_file, header = F)

k3_admix_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt'
k3_admix_res_file <- fread(k3_admix_res_file)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_pca/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.PCAresults.rds'
pca_res <- readRDS(pca_res_file)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)


################3



quit(save = 'no')


