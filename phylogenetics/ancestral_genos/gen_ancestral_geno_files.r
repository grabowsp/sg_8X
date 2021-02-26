# Script to generate ancestral allele genotype files using data from 
#  John

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### INPUT ARGUEMENTS ###
args = commandArgs(trailingOnly=T)
chr_name <- args[1]

### LOAD PACKAGES ###
library(data.table)

### IMPORT DATA ###
data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/'
ancestral_file <- paste(data_dir, chr_name, '.sg_8x_ancestral_genos.txt', 
  sep = '')
ancest_info <- fread(ancestral_file, header = F)
colnames(ancest_info) <- c('CHR', 'POS', 'ANCESTRAL', 'AP13')
ancest_info <- ancest_info[-which(duplicated(ancest_info$POS))]

ra_file <- paste(data_dir, chr_name, '.ancestral_snp_ref_alt.txt', sep = '')

ra_info <- fread(ra_file, header = F)
colnames(ra_info) <- c('CHR', 'POS', 'ID', 'REF', 'ALT')

# QC #
if(sum(ancest_info$POS != ra_info$POS) > 0){
  print('SNP positions differ between input files')
  quit(save = 'no')
}
if(sum(ancest_info$AP13 != ra_info$REF) >0){
  print('1 or more ancestral file AP13 genotypes do not match REF')
  quit(save = 'no')
}

### SET OUTPUTS ###
combo_file_out <- paste(data_dir, chr_name, 
  '.sg_8X_processed_ancestral_info.rds', sep = '')
dip_vcf_out <- paste(data_dir, chr_name,
  '.ancestral_AP13_VCF_dip_format.txt', sep = '')
tet_vcf_out <- paste(data_dir, chr_name,
  '.ancestral_AP13_VCF_tet_format.txt', sep = '')
kept_snp_out <- paste(data_dir, chr_name,
  '.ancestral_AP13_kept_SNP_positions.txt', sep = '')

### SET VARIABLES ###
dip_hom_ref_vcf <- '0/0:20,0'
dip_hom_alt_vcf <- '1/1:0,20'

tet_hom_ref_vcf <- '4/0:20,0'
tet_hom_alt_vcf <- '0/4:0,20'

##########

combo_info <- merge(ancest_info, ra_info)

saveRDS(combo_info, combo_file_out)

# remove SNPs where Ancestral genotype is neither REF or ALT
anc_alt_inds <- which(combo_info$ANCEST != combo_info$AP13)

mismatch_ancest_inds <- intersect(which(combo_info$ANCEST != combo_info$ALT),
  anc_alt_inds)

combo_info_1 <- combo_info[-mismatch_ancest_inds]

combo_info_1[, ANCESTRAL_TET_GENO := as.character(NA)]
combo_info_1[, ANCESTRAL_DIP_GENO := as.character(NA)]
combo_info_1[, REF_TET_GENO := as.character(NA)]
combo_info_1[, REF_DIP_GENO := as.character(NA)]

combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$REF),
  ANCESTRAL_TET_GENO := tet_hom_ref_vcf]
combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$REF),
  ANCESTRAL_DIP_GENO := dip_hom_ref_vcf]

combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$ALT),
  ANCESTRAL_TET_GENO := tet_hom_alt_vcf]
combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$ALT),
  ANCESTRAL_DIP_GENO := dip_hom_alt_vcf]

combo_info_1[, REF_TET_GENO := tet_hom_ref_vcf]
combo_info_1[, REF_DIP_GENO := dip_hom_ref_vcf]

combo_info_1[, QUAL := '.']
combo_info_1[, FILTER := '.']
combo_info_1[, INFO := '.']
combo_info_1[, FORMAT := 'GT:AD']

tet_out_table <- combo_info_1[, c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
  'FILTER', 'INFO', 'FORMAT', 'ANCESTRAL_TET_GENO', 'REF_TET_GENO')]
colnames(tet_out_table)[1] <- '#CHROM'

fwrite(tet_out_table, file = tet_vcf_out, sep = '\t')

dip_out_table <- combo_info_1[, c('CHR', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
  'FILTER', 'INFO', 'FORMAT', 'ANCESTRAL_DIP_GENO', 'REF_DIP_GENO')]
colnames(dip_out_table)[1] <- '#CHROM'
fwrite(dip_out_table, file = dip_vcf_out, sep = '\t')

snp_pos_table <- combo_info_1[, c('CHR', 'POS')]
fwrite(snp_pos_table, file = kept_snp_out, sep = '\t')

quit(save = 'no')

