#
#

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

library(adegenet)

gen_function_file <- '/global/homes/g/grabowsp/tools/sg_8X/general_r_tools/general_functions.r'
source(gen_function_file)

adeg_function_file <- '/global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/adegenet_functions.r'
source(adeg_function_file)



vcf_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_southcoastal_tet_vcfs/Chr09K.tetrasomic.CDS.geobig_southcoastal.vcf_03'

vcf_0 <- read.table(vcf_file, header = F, stringsAsFactors = F)

header_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_southcoastal_tet_vcfs/CDS.tetrasomic.geobig_southcoastal.vcf.sampheader.txt'

vcf_header <- gsub('#', '', read.table(header_file, stringsAsFactors = F,
  sep = '\t', header = F, comment.char = '@'))

colnames(vcf_0) <- vcf_header

tet_lib_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/tet_samps_Nov2020.txt'
tet_libs_0 <- as.vector(read.table(tet_lib_file, header = F,
  stringsAsFactors = F)[,1])
tet_libs_1 <- intersect(tet_libs_0, vcf_header)

oct_lib_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/oct_samps_Nov2020.txt'
oct_libs_0 <- as.vector(read.table(oct_lib_file, header = F,
  stringsAsFactors = F)[,1])
oct_libs_1 <- intersect(oct_libs_0, vcf_header)

tet_libs <- tet_libs_1
oct_libs <- oct_libs_1

vcf <- vcf_0[c(1:1000), ]

tot_samp_vec <- c(oct_libs, tet_libs)
  tot_df_tmp <- vcf[, tot_samp_vec]
  tot_df <- data.frame(apply(tot_df_tmp, 2, function(y)
      unlist(lapply(strsplit(y, split = ':'), function(x) x[2]))),
      stringsAsFactors = F)
    tot_df[tot_df == '0,0'] <- NA




as.numeric(unlist(strsplit(tot_df[1,1], split = ',')))

sample(x = apply(as.numeric(unlist(strsplit(tot_df[1,1], split = ','))), 1, 
  function(x) c(rep(0, times = x[1]), rep(1, times = x[2])))


make_allele_vec <- function(read_vec){
  if(is.na(read_vec)){
    out_vec <- NA
  } else{
    out_vec <- c(rep(0, times = read_vec[1]), rep(1, times = read_vec[2]))
  }
  return(out_vec)
}

sample(
  x = make_allele_vec(as.numeric(unlist(strsplit(tot_df[4,3], split = ',')))),
  size = 1)

test <- apply(tot_df, 2, function(y) sample(
  x = make_allele_vec(as.numeric(unlist(strsplit(y, split = ',')))), size = 1))


assign_pseudohap_geno <- function(count_char){
  if(is.na(count_char)){
    pseudohap_geno <- NA
  } else{
    count_vec <- as.numeric(unlist(strsplit(count_char, split = ',')))
    allele_vec <- c(rep(0, times = count_vec[1]), rep(1, times = count_vec[2]))
    pseudohap_geno <- sample(allele_vec, size = 1)
  }
  return(pseudohap_geno)
}

test[] <- lapply(tot_df, function(x) assign_pseudohap_geno(x))
assign_pseudohap_geno(tot_df[1,1])

test_1 <- apply(tot_df, c(1,2), assign_pseudohap_geno)

maf_cut <- 6/ncol(tot_df)

test <- gen_pseudohap_gl_preobj(vcf = vcf_0[1:2000,], oct_libs = oct_libs_1,
  tet_libs = tet_libs_1, maf_cut = 6/length(c(oct_libs_1, tet_libs_1)))

test_2 <- gen_gl_object(test)
