# Make Sample name lists for "all_samp", "geo_perm" and "geo_big" sample sets

# all_samp = all 1058 libraries Sujan used for making VCFs, after removing
#              the duplicated PLANT_IDs
# geo_perm = natural collections, cultivars, and MX samps - for filtering
#               and checking for cultivar escapees
# geo_big = natural collections, but without other filtering that will probably
#            eventually happen

# on NERSC
# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

library(tidyverse)
library(data.table)

### LOAD DATA ###
meta_file <- paste('/global/homes/g/grabowsp/data/switchgrass/metadata_8X/', 
  'sg_8X_metadata_v1.0.csv', sep = '')
samp_meta <- fread(meta_file)

### SET OUTPUTS ###
out_dir <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/'

all_samp_short <- 'all_samp_names.txt'
all_samp_out <- paste(out_dir, all_samp_short, sep = '')

geo_perm_short <- 'geo_perm_names.txt'
geo_perm_out <- paste(out_dir, geo_perm_short, sep = '')

geo_big_short <- 'geo_big_names.txt'
geo_big_out <- paste(out_dir, geo_big_short, sep = '')
################

all_samp_names <- samp_meta$VCF_NAME[-which(samp_meta$PLANT_ID == 'AP13')]

write.table(all_samp_names, all_samp_out, quote = F, sep = '\t', row.names = F,
  col.names = F)

geo_perm_names <- samp_meta[COLLECTION_TYPE == 'Natural Collection' | 
  COLLECTION_TYPE == 'Cultivar' | COLLECTION_TYPE == 'Mexican Samples', 
  VCF_NAME]

write.table(geo_perm_names, geo_perm_out, quote = F, sep = '\t', row.names = F,
  col.names = F)

geo_big_names <- samp_meta[COLLECTION_TYPE == 'Natural Collection', VCF_NAME]

write.table(geo_big_names, geo_big_out, quote = F, sep = '\t', row.names = F,
  col.names = F)

quit(save = 'no')

