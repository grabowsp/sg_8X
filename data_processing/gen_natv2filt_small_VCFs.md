# Generate sub-VCFs for testing reproducibility

## Overview
* Make 1000 reps of 100k randomly chosen SNPs from the natv2filt SNP set
* Make 100 tetrasomic 100k VCFs
  * only need to make tetrasomic VCFs because these are used to generate
all downstream objects, whether disomic, tetrasomic, or polyploid

## Result Files
* Directory with sub-VCFs
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs`
* Directory with SNP subfiles
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/snp_subset_files`

## Important Files
* natv2filt SNPs
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/GW.disomic.CDS.natv2.filt.keptSNPs.txt`
* natv2filt full-VCF directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs`

## Generate Sets of SNPs
* in R
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

# in R

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
tot_snp_pos_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/GW.disomic.CDS.natv2.filt.keptSNPs.txt'

tot_snp_info <- fread(tot_snp_pos_file)

### SET OUTPUT ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/snp_subset_files/'

out_pre <- 'natv2filt.'

### SET VARIABLES ###
n_reps <- 1000

n_snps <- 1e5
n_snp_char <- paste(n_snps/1000, 'k.subsnps.', sep = '')
#############

for(i in seq(n_reps)){
  sub_inds <- sort(sample(seq(nrow(tot_snp_info)), size = n_snps))
  sub_snp_info <- tot_snp_info[sub_inds, ]
  n_rep_char <- sprintf("%04d", i)
  out_full <- paste(out_dir, out_pre, n_snp_char, n_rep_char, '.pos.txt',
    sep = '')
  fwrite(sub_snp_info, file = out_full, col.names = F, sep = '\t')
}

quit(save = 'no')

```

## Generate sub-VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs

sbatch gen_natv2filt_tet_subfiles_01_10.sh
sbatch gen_natv2filt_tet_subfiles_11_50.sh
sbatch gen_natv2filt_tet_subfiles_51_100.sh

sbatch gen_natv2filt_tet_subfiles_92_100.sh
```
### Example script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

IN_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/

IN_SUF=.tetrasomic.CDS.natv2filt.vcf.gz

SNP_POS_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/snp_subset_files/

SNP_POS_PRE=natv2filt.100k.subsnps.

NSNPS=100k

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/

OUT_PRE=GW.natv2filt.tet.

cd $OUT_DIR

for SNPSET in {0001..0010};
  do
  VCF_OUT=$OUT_PRE$NSNPS'.'$SNPSET'.vcf';
  SNP_SUB_FILE=$SNP_POS_DIR$SNP_POS_PRE$SNPSET'.pos.txt';
  #
  CHROM=Chr01K;
  echo $CHROM $SNPSET;
  VCF_IN=$IN_DIR$CHROM$IN_SUF;
  vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
    --recode --recode-INFO-all > $VCF_OUT;
  #
  CHROM=Chr01N;
  echo $CHROM $SNPSET;
  VCF_IN=$IN_DIR$CHROM$IN_SUF;
  vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
    --recode --recode-INFO-all | tail -n +635 >> $VCF_OUT;
  #
  for CHR_NUM in {02..09};
    do
    for CHR_LET in K N;
      do
      CHROM=Chr$CHR_NUM$CHR_LET;
      echo $CHROM;
      VCF_IN=$IN_DIR$CHROM$IN_SUF;
      vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
        --recode --recode-INFO-all | tail -n +635 >> $VCF_OUT;
      done;
    done;
  gzip $VCF_OUT;
  done;

```
