# Script for generating v3 of metadata with notes about native samples and
#  contaminants

### LOAD MODULES ###

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v2.0.csv'
samp_meta <- fread(meta_file)

#filt_info_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/SG_8X_Problematic_Samples.tsv'
filt_info_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/SG_8X_Problematic_Samples_v2.tsv'
filt_info <- fread(filt_info_file)

geobig_samp_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_big_names.txt'
geobig_names <- unlist(fread(geobig_samp_file, header = F)[,1])

### SET OUTPUT ###

#out_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.0.csv'
out_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'

########

colnames(filt_info) <- gsub(' ', '_', colnames(filt_info))

filt_info <- filt_info[-which(duplicated(filt_info$Sample_Name))]

filt_info[Sample_Name == 'Cave_in_Rock_WO1', Sample_Name := 'Cave_In_Rock_WO1']
filt_info[Sample_Name == 'J450_A', Sample_Name := 'J450.A']

length(intersect(filt_info$Sample_Name, samp_meta$VCF_NAME))
setdiff(filt_info$Sample_Name, samp_meta$VCF_NAME)

meta_ind <- c()
for(i in seq(nrow(filt_info))){
  tmp_ind <- which(samp_meta$VCF_NAME == filt_info$Sample_Name[i])
  meta_ind <- c(meta_ind, tmp_ind)
}

samp_meta[, PPG_FILT_NOTE := as.character(NA)]
samp_meta[, PPG_NATIVE := as.character(NA)]
samp_meta[, PPG_CONTAMINANT := as.character(NA)]

samp_meta[meta_ind, PPG_FILT_NOTE := filt_info$Description]
samp_meta[meta_ind, PPG_NATIVE := filt_info$PPG_NATIVE]
samp_meta[meta_ind, PPG_CONTAMINANT := filt_info$PPG_CONTAMINANT]

length(intersect(samp_meta[PPG_NATIVE == 'F', VCF_NAME], geobig_names))
# 115

length(setdiff(geobig_names, samp_meta[PPG_NATIVE == 'F', VCF_NAME]) )
# 706 left

fwrite(samp_meta, out_file, na = 'NA')

quit(save = 'no')

