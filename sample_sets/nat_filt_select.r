# Generate 'nat_filt_1' which are "Natural Collection" Samples that were
#  further filtered by PPG based on NJ tree clustering using 'allsamps'

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES #####
library(data.table)

### INPUT DATA ###
#meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.0.csv'
meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

geobig_samp_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_big_names.txt'
geobig_names <- unlist(fread(geobig_samp_file, header = F)[,1])

### SET VARIABLES ###

### SET OUTPUT ###
out_dir <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/'

#nat_short <- 'nat_filt_1_names.txt'
nat_short <- 'nat_filt_v2_names.txt'
nat_out <- paste(out_dir, nat_short, sep = '')

################3
nat_names <- setdiff(geobig_names, samp_meta[PPG_NATIVE == 'F', VCF_NAME])

# Save names

write.table(nat_names, nat_out, quote = F, sep = '\t', row.names = F,
  col.names = F)

quit(save = 'no')


