# Generate Full chromosome genlight objects with adegenet

#module load python/3.7-anaconda-2019.07
#source activate r_adegenet_env

args = commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = F)

### LOAD PACKAGES ###

library(adegenet)
library(parallel)

# get locations of ancillary github files
file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.dir <- dirname(script.name)

adeg_function_file <- file.path(script.dir, 'adegenet_functions.r')
source(adeg_function_file)

gen.dir <- sub('/adegenet_analysis/adegenet_genotype_generation', 
  '/general_r_tools', script.dir)
gen_function_file <- file.path(gen.dir, 'general_functions.r')
source(gen_function_file)

# on Cori
#gen_function_file <- '/global/homes/g/grabowsp/tools/sg_8X/general_r_tools/general_functions.r'
#source(gen_function_file)

#adeg_function_file <- '/global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/adegenet_functions.r'
#source(adeg_function_file)

### IMPORT DATA ###

data_dir <- args[1]
#data_dir <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/'

data_dir <- add_slash(data_dir)

vcf_search_string <- args[2]
#vcf_search_string <- 'Chr01K.tetrasomic.CDS.allsamps.vcf_*'

vcf_search_command <- paste('ls ', data_dir, vcf_search_string, sep = '')

vcf_files <- system(vcf_search_command, inter = T)

#vcf_files <- system('ls /global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/Chr01K.polyploid.CDS.geosamps.vcf_*', inter = T)

vcf_1 <- read.table(vcf_files[1], header = F, stringsAsFactors = F, 
    sep = '\t')

vcf_header_file <- args[3]
#vcf_header_file <- '/global/cscratch1/sd/grabowsp/sg_ploidy/polyploid_vcfs/CDS_vcfs/geo_samps/CDS.geosamps.vcf.header.txt'

vcf_header <- gsub('#', '', read.table(vcf_header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf_1) <- vcf_header

tet_lib_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/tet_samps_Nov2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs_tmp <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/oct_samps_Nov2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs_1 <- intersect(oct_libs_0, vcf_header)

# try set all libraries as disomic (4X)
tet_libs_1 <- c(tet_libs_tmp, oct_libs_1)

### SET OUTPUTS ###
out_name_short <- gsub('vcf_*', 'genlight.rds', vcf_search_string, fixed = T)
out_name_full <- paste(data_dir, out_name_short, sep = '')

### SET VARIABLES ###
maf_cut_in <- args[4]
#maf_cut <- 0.002
maf_cut_in <- as.numeric(maf_cut_in)
###########################
preobj_1 <- gen_gl_preobj(vcf = vcf_1, oct_libs = oct_libs_1, 
  tet_libs = tet_libs_1, maf_cut = maf_cut_in)

gl_tot <- gen_gl_object(preobj_list = preobj_1)

if(length(vcf_files) > 1){
  for(i in c(2:length(vcf_files))){
    vcf_tmp <- read.table(vcf_files[i], header = T, stringsAsFactors = F,
      sep = '\t')
    colnames(vcf_tmp) <- vcf_header
    preobj_tmp <- gen_gl_preobj(vcf = vcf_tmp, oct_libs = oct_libs_1, 
      tet_libs = tet_libs_1, maf_cut = maf_cut_in)
    gl_tmp <- gen_gl_object(preobj_list = preobj_tmp)
    gl_tot <- cbind(gl_tot, gl_tmp)
    print(i)
  }
}

saveRDS(gl_tot, out_name_full)

quit(save = 'no')

