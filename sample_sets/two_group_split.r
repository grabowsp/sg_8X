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

### SET VARIABLES ###

group_cutoff <- 0.9
# I chose 0.9 because higher cutoffs started filtering out samples that
#  were definitely part of the Gulf group; including them in PCA/ADMIXTURE
#  may be informative

### SET OUTPUT ###
out_dir <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/'

n_in_short <- 'geobig_NorthInland_names.txt'
n_in_out <- paste(out_dir, n_in_short, sep = '')

s_coast_short <- 'geobig_SouthCoastal_names.txt'
s_coast_out <- paste(out_dir, s_coast_short, sep = '')

################3


group_inds <- gen_admix_popind_list(k2_admix_res, group_cutoff)

# Group 1
samp_meta[VCF_NAME %in% k2_admix_res$V1[group_inds[['g1']]], .N, by=STATE]
# "S_Coastal" = "Lowland", co_0.9 = 448 samples; co_0.95 = 390; co_0.99 = 294,
#   but no Gulf 
samp_meta[VCF_NAME %in% k2_admix_res$V1[group_inds[['g1']]], .N, 
  by=NQUIRE_PLOIDY]
# Mainly 4X, <20 8X/6X

# Group 2
samp_meta[VCF_NAME %in% k2_admix_res$V1[group_inds[['g2']]], .N, by=STATE]
# "N_Inland" = 'Upland'. co_0.9 = 265 samples; co_0.95 = 259; co_0.99 = 232
samp_meta[VCF_NAME %in% k2_admix_res$V1[group_inds[['g2']]], .N, 
  by=NQUIRE_PLOIDY]
# ~50/50 4X/8X

# look at PCA results
pca_mat <- pca_res$scores

summary(pca_mat[rownames(pca_mat) %in% k2_admix_res$V1[group_inds[['g1']]], 1])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  3.694   5.333  10.532   8.906  11.251  12.211

summary(pca_mat[rownames(pca_mat) %in% k2_admix_res$V1[group_inds[['g2']]], 1])
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -17.064 -16.138 -14.442 -14.679 -13.889  -5.674 

# Save names

s_coast_names <- k2_admix_res$V1[group_inds[['g1']]]

write.table(s_coast_names, s_coast_out, quote = F, sep = '\t', row.names = F,
  col.names = F)

n_in_names <- k2_admix_res$V1[group_inds[['g2']]]

write.table(n_in_names, n_in_out, quote = F, sep = '\t', row.names = F, 
  col.names = F)

quit(save = 'no')


