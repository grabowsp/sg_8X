# Script used to check if SNPs are same order for disomic and tetrasomic
#  VCFs
#   if not, then will need to sort the tetrasomic VCF

# RESULTS - CANNOT USE DISOMIC INDEX FILES FOR TETRASOMIC FILES
### The order of the SNPs is the same, but there seem to be
###   occasional missing SNPs; ex: 1 SNP missing disomic (?) file

##################

# on NERSC
# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

#############
# check the length of the header

# cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs
# gunzip -kc Chr01K.disomic.CDS.allsamps.vcf.gz | head -648 | \
# tail -1
### samples

# head -648 Chr01K.disomic.CDS.allsamps.vcf_00 | tail -1
### samples

# cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs
# gunzip -kc Chr01K.tetrasomic.CDS.allsamps.vcf.gz | head -634 | \
# tail -1
### samples
### 14 lines shorter header

# head -634 Chr01K.tetrasomic.CDS.allsamps.vcf_00 | tail -1
### samples
#################

library(tidyverse)
library(data.table)

dis_file_00 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/Chr01K.disomic.CDS.allsamps.vcf_00'
dis_00 <- fread(cmd = paste('grep -v "^#" ', dis_file_00, sep = ''))

tet_file_00 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/Chr01K.tetrasomic.CDS.allsamps.vcf_00'
tet_00 <- fread(cmd = paste('grep -v "^#" ', tet_file_00, sep = ''))


#V2 = position
sum(dis_00$V2 != tet_00$V2[c(1:(nrow(tet_00)-14))])
# [1] 0
# all position elements are the same in file 00

dis_file_10 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/Chr01K.disomic.CDS.allsamps.vcf_10'
dis_10 <- fread(cmd = paste('grep -v "^#" ', dis_file_10, sep = ''))

tet_file_10 <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/Chr01K.tetrasomic.CDS.allsamps.vcf_10'
tet_10 <- fread(cmd = paste('grep -v "^#" ', tet_file_10, sep = ''))

sum(dis_10$V2[14:nrow(dis_10)] != tet_10$V2)
# [1] 0

dis_short <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/Chr01K.disomic.CDS.allsamps.vcf_0'
tet_short <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/Chr01K.tetrasomic.CDS.allsamps.vcf_0'

for(i in c(1:9)){
  dis_file <- paste(dis_short, i, sep = '')
  tet_file <- paste(tet_short, i, sep = '')
  tmp_dis <- fread(cmd = paste('grep -v "^#" ', dis_file, sep = ''))
  tmp_tet <- fread(cmd = paste('grep -v "^#" ', tet_file, sep = ''))
  #
  print(i)
  n_missing <- sum(tmp_dis$V2[c(15:nrow(tmp_dis))] != 
    tmp_tet$V2[c(1:(nrow(tmp_tet)-14))])
  print(n_missing)
}

[1] 1
[1] 0
[1] 2
[1] 0
[1] 3
[1] 0
[1] 4
[1] 91092
[1] 5
[1] 99986
[1] 6
[1] 99986
[1] 7
[1] 99986
[1] 8
[1] 99986
[1] 9
[1] 99986



