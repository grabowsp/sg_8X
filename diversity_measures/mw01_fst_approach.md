# Steps for testing measuring Fst between MW01 and MW01_hi
* To test if higher Gulf ancestry is from second pulse

## Comparisons
* MW01 vs MW01_hi
  * sites private to MW
  * sites shared by MW and Gulf
  * hiFst sites
  * sites shared by MW and Atlantic but not Gulf

## Steps
1. Generate allele frequency files for MW, GULF, and ATL
  1. Generate Sample sets for each gene pool
    * `/home/f2p1/work/grabowsk/data/switchgrass/mw01_analysis/MW_lib_names_all.txt`
    * `/home/f2p1/work/grabowsk/data/switchgrass/mw01_analysis/GULF_lib_names_all.txt`
    * `/home/f2p1/work/grabowsk/data/switchgrass/mw01_analysis/ATL_lib_names_all.txt`
  2. vcfTools to output SNPs with frequency > 0.01
2. Generate SNP sets
  1. Private to MW
  2. Shared between MW and GULF
  3. MWvGulf hiFst sites
  4. sites shared by MW and Atlantic but NOT Gulf
3. Calculate Fst between MW_01 and MW_01_hi for each SNP set
4. Controls
  1. Fst between MW_01 and MW_03
  2. Fst between 2 groups within MW_03

## Generate Genepool SNP sets
### Generate sample name files
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
result_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt', sep = '')
res_tab <- fread(result_file)

### SET OUTPUTS ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/mw01_analysis/'

mw_lib_out <- paste(out_dir, 'MW_lib_names_all.txt', sep = '')
gulf_lib_out <- paste(out_dir, 'GULF_lib_names_all.txt', sep = '')
atl_lib_out <- paste(out_dir, 'ATL_lib_names_all.txt', sep = '')

########
mw_libs <- res_tab[grep('MW_', res_tab$subgrp_v2), list(samp_name)]
# 228
fwrite(mw_libs, file = mw_lib_out, col.names = F)

gulf_libs <- res_tab[grep('GULF_', res_tab$subgrp_v2), list(samp_name)]
# 157
fwrite(gulf_libs, file = gulf_lib_out, col.names = F)

atl_libs <- res_tab[grep('ATL_', res_tab$subgrp_v2), list(samp_name)]
# 284
fwrite(atl_libs, file = atl_lib_out, col.names = F)

```

