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

