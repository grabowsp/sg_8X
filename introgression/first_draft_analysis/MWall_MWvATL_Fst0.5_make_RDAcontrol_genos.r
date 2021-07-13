# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

vcf_func_file <- '/home/grabowsky/tools/workflows/sg_8X/introgression/vcf_analysis_functions.r'
source(vcf_func_file)

### INPUT DATA ###

control_vcf <- '/home/f1p1/tmp/switchgrass_8X/MWall_control_vcfs/ATLvMW_Fst_0.5_MWall_RDAcontrol.tot.1.vcf'
vcf_in <- read_vcf(control_vcf)

### SET OUTPUTS ###
out_dir <- '/home/f1p1/tmp/switchgrass_8X/MWall_control_vcfs/'

out_rda_file <- paste(out_dir, 'ATLvMW_Fst_0.5_MWall_RDAcontrol.genos.1.rds',
  sep = '')

############3

# Generate table with genotypes
test_genos <- vcf_many_SNP_genotypes(vcf_in, snp_name_vec = vcf_in$ID)

# adjust genotypes to the number of reference alleles
test_genos[test_genos == './.'] <- NA
test_genos[test_genos == '0/0'] <- 2
test_genos[test_genos == '0/1'] <- 1
test_genos[test_genos == '1/1'] <- 0

# generate matrix in format of rows = Samples, columns = SNPs
geno_mat <- as.matrix(test_genos[, 2:ncol(test_genos)])
# this also transposes the matrix into the correct orientation
geno_mat_1 <- apply(geno_mat, 1, as.numeric)
rownames(geno_mat_1) <- colnames(geno_mat)
colnames(geno_mat_1) <- test_genos$ID

saveRDS(geno_mat_1, file = out_rda_file)

