# Steps used for Selecting samples for Different Sample Set

## `all_samp`
* all 1058 libraries for unique PLANT_ID's that Sujan used for making the VCFs
* VCF sample name list:
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt`

## `geo_perm`
* 936 samples
  * 821 'Natural Collection', 73 'Cultivar', and 42 'Mexican samples' in `COLLECTION_TYPE` column
* to be used for seeing if any of the natural collection samples are very close to cultivars and/or are most relavant to highly-used cultivars
* VCF sample name list:
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_perm_names.txt`

## `geo_big`
* 821 'Natural Collection' samples in `COLLECTION_TYPE` column
  * This set will probably get whittled down do a more refined set as we do more filtering
* VCF sample name list:
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_big_names.txt`

## `geobig_southcoastal`
* 448 samples with >0.9 membership in "group 1" in K=2 ADMIXTURE results using `geo_big sample`s
  * These correspond to previous "lowland" genetic designation
* VCF sample name list:
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geobig_SouthCoastal_names.txt` 

## `geobig_northinland`
* 265 samples with >0.9 membership in "group 2" in K=2 ADMIXTURE results using geo_big samples
  * these correspond to previous "upland" genetic designation
* VCF sample name list
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geobig_NorthInland_names.txt`

## `nat_filt_v2`
* 706 samples
* `geobig` samples that were further filtered based on clustering in NJ joining tree using `allsamps`
  * 115 samples removed
  * `nat_filt_1` had 5 fewer samples removed
* VCF sample name list
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt`

## `natv2filt` temporary genepool designations
* used ADMIXTURE and PCA results to id samples for each gene pool
  * "main" = main cluster
  * "total" = expaned cluster
* Midwest main
  * 206 samples
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_MW_main_names.txt`
* Midwest total
  * 228 samples
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_MW_total_names.txt`
* Gulf main
  * 77 samples
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_GULF_main_names.txt`
* Gulf total
  * 157 samples
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_GULF_total_names.txt`
* Atlantic main
  * 257 samples
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_ATL_main_names.txt`
* Atlantic total
  * * 284 samples
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_ATL_total_names.txt`
* "Misfit" samples
  * 37 samples that don't fit in any of the "total" genepool groups
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/natv2filt_tmp_MISFIT_names.txt`

* `natv2filt` subgroup names
  * Genepools and Admixed samples divided into sub groups
    * Groups decided using PCA and ADMIXTURE results
    * 10+ groups in each genepool, 8 Admixed groups
  * git directory: `~sg_8X/sample_sets/samp_set_files/natv2filt_subgrps/`

