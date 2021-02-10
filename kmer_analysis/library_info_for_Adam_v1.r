# Steps for generating library info for Adam to try Kmer analysis

# Goal:
# 1) Select groups of related 8X
#  - Use Joes groups but select only those in "native" group
# 2) Select closest 4X
# 3) Select "control" 4X
# 4) Get JGI library ID's for all samples

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)

### IMPORT DATA ###
meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

natv2_samp_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_samp_file, header = F)[, V1]

test_samp_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/'

atlantic_file <- paste(test_samp_dir, 'Atlantic_ref_samps_40.txt', sep = '')
atl_samps <- fread(atlantic_file, header = T)[, ID]

gulf_file <- paste(test_samp_dir, 'Gulf_ref_samps_40.txt', sep = '')
gulf_samps <- fread(gulf_file, header = T)[, ID]

mw_file <- paste(test_samp_dir, 'Midwest_ref_samps_40.txt', sep = '')
mw_samps <- fread(mw_file, header = T)[, ID]

candidate_samp_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/Introgression_target_individuals.txt'
candidate_samps <- fread(candidate_samp_file)

grp3_file <- paste(test_samp_dir, 'Grp4_subgroups.csv', sep = '')
grp3_samps <- fread(grp3_file)

### SET OUTPUTS ###


####################

grp_1_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp1', ID], 
  natv2_samps)

# Atlantic ancestry; most have > 10% Gulf ancestry, but not all...
# 13/13 are in filtered sample set
# J192.L1 is separate from other grp_1 samps and clusters with other samples
#  from same ACC
# J593.A is probably 6X and is separate from other grp_1 samps and clusters
#   with other sample from same ACC
# J180.L2 clusters within its own ACC

grp_2_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp2', ID], 
  natv2_samps)

# grp_2 is 8X-East
# 66/81 are in filtered sample set

grp_4_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp4', ID],
  natv2_samps)
# grp_4 is 8X-West
# 41/70 are in filtered sample set

grp_5_filt <- intersect(candidate_samps[kmeans_grps_five == 'grp5', ID],
  natv2_samps)

# grp_5 is Gulf ancestry
# 18/21 are in filtered sample set

# J462-J309 form a column in the allsamps tree
# J249.A is on it's own - it might be 6X, but nQuire results aren't great;
#  though old results suggest it's 6X
# J188.L1 is on it's own, clustered within J188

grp_3_1_filt <- intersect(grp3_samps[Sub_Group == 'subgrp1', ID],
  natv2_samps)
# 1/1
# J479_B - Arkansas

grp_3_2_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp2', ID],
  natv2_samps)
# 10/10
# Mississippi 8X

grp_3_3_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp3', ID],
  natv2_samps)
# 12/12
# GA (and 1 NC) 8X

grp_3_4_filt <-  intersect(grp3_samps[Sub_Group == 'subgrp4', ID],
  natv2_samps)
# 1/1
# J530.B
# NY, part of cluster with other J530

####

grp_1_1_8X <- setdiff(grp_1_filt, c('J192.L1', 'J593.A', 'J180.L2'))
# FL 8X with Atlantic background
# 55-85% Atlantic, rest mainly Gulf
grp_1_1_4X_1 <- intersect(natv2_samps, 
  c(samp_meta[ACC == 'J187', VCF_NAME],
  samp_meta[ACC == 'J184', VCF_NAME],
  samp_meta[ACC == 'J185', VCF_NAME],
  samp_meta[ACC == 'J186', VCF_NAME])
  )

samp_meta[ACC == 'J497', VCF_NAME],
samp_meta[ACC == 'J496', VCF_NAME],
samp_meta[ACC == 'J188', VCF_NAME],
samp_meta[ACC == 'J504', VCF_NAME]



quit(save = 'no')

