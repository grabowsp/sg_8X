# Steps for metadata management for 8X switchgrass analysis

## Initial metadata file from Nature analysis:
* John's metadata file
  * Updated May 2020; included values used for paper analysis
  * On NERSC:
    * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/genotype.metadata.May2020.rds`
  * On HA:
    * `/home/t4c1/WORK/grabowsk/data/switchgrass/genotype.metadata.May2020.rds`
* Most recent version from Jason
  * I saved both .txt and the original excel versions
  * On NERSC:
    * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/PVDIV_Master_Metadata_File_7-1-20.txt`
    * `/global/homes/g/grabowsp/data/switchgrass/reseq_metadata/PVDIV_Master_Metadata_File_7-1-20.xlsx`
  * on HA:
    * `/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_7-1-20.txt`
    * `/home/t4c1/WORK/grabowsk/data/switchgrass/PVDIV_Master_Metadata_File_7-1-20.xlsx`

## Accessory files used for generating metadata versions
* Notes from all_samp NJ trees about sample filtering
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/SG_8X_Problematic_Samples_v2.tsv`

## Directory with Metadata
* On NERSC
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X`
* On HA
  * `/home/f2p1/work/grabowsk/data/switchgrass`

## Metadata versions
  * `sg_8X_metadata_v1.0.csv`
    * First assembly of info. Includes Library ID's, Plant ID's, names used in the VCF, info from Jason's metadata, sequencing statistics, and ploidy info
  * `sg_8X_metadata_v2.0.csv`
    * Each JXXX accession assigned a unified accession ('UNI_ACC') to lump together accessions that are supposed to be of the same origin.
  * `sg_8X_metadata_v3.1.csv`
    * Added sample filtering calls and notes about if samples are native or contaminants
    * `...v3.0.csv` was a previous version with 5 fewer samples called non-native


