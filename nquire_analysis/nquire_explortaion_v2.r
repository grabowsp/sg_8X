# Script for generating figures for initial nQuire results

# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)
library(patchwork)

### LOAD DATA ###

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

allsamps_lib_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt'
allsamps_names <- unlist(fread(allsamps_lib_file))

geobig_lib_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_big_names.txt'
geobig_names <- unlist(fread(geobig_lib_file))

### SET OUTPUTS ###

out_dir <- paste('/global/cscratch1/sd/grabowsp/sg_8X_scratch/',
  'nquire_scratch/', sep = '')

allsamp_bar_out <- paste(out_dir, 'allsamps_ploidy_barplot.pdf', sep = '')

geobig_dot_out <- paste(out_dir, 'geobig_ploidy_vs_coverage.pdf', sep = '')

#########

meta_tab_1 <- samp_meta[is.na(NQUIRE_PLOIDY) == F, .N, by = .(NQUIRE_PLOIDY,
  COVERAGE_CATEGORY)]

gg_all_COV_bar <- ggplot(meta_tab_1, aes(x = NQUIRE_PLOIDY, y = N, 
  fill = COVERAGE_CATEGORY)) + geom_bar(stat = 'identity') +
  ggtitle('nQuire-based ploidy for all 1058 samples\nby sample set')

meta_tab_2 <- samp_meta[is.na(NQUIRE_PLOIDY) == F, .N, by = .(NQUIRE_PLOIDY,
  COLLECTION_TYPE)]

gg_all_COL_bar <- ggplot(meta_tab_2, aes(x = NQUIRE_PLOIDY, y = N, 
  fill = COLLECTION_TYPE)) + geom_bar(stat = 'identity') +
  ggtitle('nQuire-based ploidy for all 1058 samples\nby collection type')

pdf(allsamp_bar_out, width = 9, height = 4)
gg_all_COV_bar + gg_all_COL_bar
dev.off()

###

gg_geobig_dot <- ggplot(samp_meta[VCF_NAME %in% geobig_names, ], 
  aes(x = APPROX_COVERAGE, y = NQUIRE_PLOIDY)) +
  geom_jitter(position = position_jitter(width = 0, height = 0.3)) +
  geom_vline(xintercept = 20, linetype = 'dashed', color = 'black') +
  ggtitle(paste('nQuire ploidy calls for 820 Natural Collection ', 
    'samples\nvs Sequencing coverage', sep = ''))

pdf(geobig_dot_out, width = 5, height = 4)
gg_geobig_dot
dev.off()


