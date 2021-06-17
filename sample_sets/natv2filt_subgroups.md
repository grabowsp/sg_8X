# Generating natv2filt subgroup files

* File used to divide the samples into subgroups:
  * `~sg_8X/results_tab_management/natv2filt/make_natv2filt_res_tab_v3.0.R`

* Git directory with name files:
  * `/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/natv2filt_subgrps`

* Generate subgroup name files for vcftools
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_tab_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.0.txt'
res_tab <- fread(res_tab_file)


### SET OUTPUTS ###
git_out_dir <- '/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/natv2filt_subgrps/'

name_file_suf <- '_names.txt'

####

sg_names <- unique(res_tab$sub_grp)

for(sg in sg_names){
  tmp_names <- res_tab[sub_grp == sg , list(samp_name)]
  out_file <- paste(git_out_dir, sg, name_file_suf, sep = '')
  fwrite(tmp_names, out_file, sep = '\t', col.names = F)
}

quit(save = 'no')
```
