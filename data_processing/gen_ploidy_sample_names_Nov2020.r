# Script for generating lists of VCF sample names for each ploidy
#
# will include the 6X samples with 8X samples, for now

# on NERSC
# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)

### LOAD DATA ###

meta_file <- paste('/global/homes/g/grabowsp/data/switchgrass/metadata_8X/',
  'sg_8X_metadata_v1.0.csv', sep = '')
samp_meta <- fread(meta_file)

### SET OUTPUT ###
out_dir <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/'

tet_samp_out <- paste(out_dir, 'tet_samps_Nov2020.txt', sep = '')
oct_samp_out <- paste(out_dir, 'oct_samps_Nov2020.txt', sep = '')

#############

tet_samp_names <- samp_meta[ NQUIRE_PLOIDY == '4X', VCF_NAME]
oct_samp_names <- samp_meta[NQUIRE_PLOIDY == '8X' | NQUIRE_PLOIDY == '6X', 
  VCF_NAME]

write.table(tet_samp_names, file = tet_samp_out, quote = F, sep = '\t',
  row.names = F, col.names = F)

write.table(oct_samp_names, file = oct_samp_out, quote = F, sep = '\t',
  row.names = F, col.names = F)

quit(save = 'no')

