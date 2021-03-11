# scratch for working out measures of homozygosity

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

args = commandArgs(trailingOnly = T)
rundir_args <- commandArgs(trailingOnly = F)

file.arg.name <- '--file='
script.name <- sub(file.arg.name, '',
  rundir_args[grep(file.arg.name, rundir_args)])
script.basename <- dirname(script.name)

gen.basename <- sub('/diversity_measures', '', script.basename)

### LOAD PACKAGES ###

library(data.table)
library(R.utils)

# div_funct_file <- '/global/homes/g/grabowsp/tools/sg_8X/diversity_measures/div_measures_functions.r'
div_funct_file <- file.path(script.basename, 'div_measures_functions.r')
source(div_funct_file)

# gen_funct_file <- '/global/homes/g/grabowsp/tools/sg_8X/general_r_tools/general_functions.r'
gen_funct_file <- file.path(gen.basename, 
  'general_r_tools/general_functions.r')
source(gen_funct_file)

### LOAD INPUTS ###

vcf_file_in <- args[1]
#vcf_file_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_tet_vcfs/Chr01K.tetrasomic.CDS.natv2.vcf_00'
#vcf_file_2 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/GW.natv2.ancestral.tet.100k.0001.vcf.gz'

vcf <- fread(vcf_in, skip = '#CHROM')
colnames(vcf)[1] <- 'CHROM'

#######
# NEED TO CONTINUE FROM HERE
#######

### SET OUTPUTs ###

out_dir <- args[2]
out_dir <- add_slash(out_dir)

out_pre <- args[3]
# out_pre <- 'Chr01.geobig.CDS'

out_total <- paste(out_dir, out_pre, '.HETvals.rds', sep = '')

### SET VARIABLES ###

block_in <- as.numeric(args[4])
#block_in <- 1e5


################

genos <- readRDS(geno_file_in)

tot_divers_list <- lapply(seploc(genos, block.size = block_in, parallel = T),
  calc_gl_het_stats)

tot_tet_hom <- apply(data.frame(lapply(tot_divers_list, function(x)
  x[['tet_hom']] * x[['n_snps']])), 1, sum) / nLoc(genos)

tot_oct_hom <- apply(data.frame(lapply(tot_divers_list, function(x)
  x[['oct_hom']] * x[['n_snps']])), 1, sum) / nLoc(genos)

tot_correct_oct_het <- apply(data.frame(lapply(tot_divers_list, function(x)
  x[['corrected_oct_het']] * x[['n_snps']])), 1, sum) / nLoc(genos)

compiled_list <- list()
compiled_list[['tot_tet_hom']] <- tot_tet_hom
compiled_list[['tot_oct_hom']] <- tot_oct_hom
compiled_list[['tot_correct_oct_het']] <- tot_correct_oct_het
compiled_list[['N_SNPs']] <- nLoc(genos)

saveRDS(compiled_list, out_total)

quit(save = 'no')


