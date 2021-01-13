# Steps for generating genotypes to be used for HMM looking for introgression

## Overview
* Use Disomic genotypes
* Identify SNPs that have high Fst between the 3 source populations
* Generate VCF containing high-Fst SNPs and all the samples

## Location of input files
* Starting VCFs
  * Parent directory
    * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs`
* Sample files
  * Parent directory
    * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try`
  * Sample files
    * `Atlantic_ref_individuals.txt`
    * `Gulf_ref_individuals.txt`
    * `Midwest_ref_individuals.txt`
    * `Introgression_target_individuals.txt`

## Choosing SNPs
### Select filtering criteria
* MAF
  * at least in 5 samples = 5/(60*2) = 0.042
  * or for pairwise: 5/(40*2) + 0.0625
* Fst values
  * use 'weir-fst-pop' flag in vcf tools
  * want Fst above 0.75
### Get Fst values with VCFTools
#### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

sbatch get_fst_vals_01_03.sh
sbatch get_fst_vals_04_06.sh
sbatch get_fst_vals_07_09.sh
```
#### Example script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 24:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/
ATLANTIC_FILE=Atlantic_ref_names.txt
GULF_FILE=Gulf_ref_names.txt
MW_FILE=Midwest_ref_names.txt

for CHR_NUM in {01..04};
  do
  for CHR_LET in K N;
    do
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    #
    MAF=0.042;
    OUT_TXT=ALL3_Fst_Chr$CHR_NUM$CHR_LET;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE \
      --weir-fst-pop $SAMP_DIR$GULF_FILE \
      --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    MAF=0.062;
    OUT_TXT=GULFvMW_Fst_Chr$CHR_NUM$CHR_LET;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$GULF_FILE --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    OUT_TXT=ATLANTICvMW_Fst_Chr$CHR_NUM$CHR_LET;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    OUT_TXT=ATLANTICvGULF_Fst_Chr$CHR_NUM$CHR_LET;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE \
      --weir-fst-pop $SAMP_DIR$GULF_FILE;
    done;
  done;
```
### Consolidate SNP list
* Final file:
  * `global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/hi_fst_SNPS_v1.txt`
  * 85979 SNPs
* in R
  * `~/sg_8X/introgression/consolidate_v1_Fst_SNPs.r`

## Make Test Sample lists
* in R
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/'

samp_file <- paste(res_dir, 'Introgression_target_individuals.txt', sep = '')

samp_info <- fread(samp_file)

grp_vec <- unique(samp_info[, kmeans_grps_five])

grp_list <- list()

for(i in grp_vec){
  grp_list[[i]] <- samp_info[kmeans_grps_five == i, ID]
}

for(j in names(grp_list)){
  tmp_out_file <- paste(res_dir, j, '_sampnames.txt', sep = '')
  write.table(grp_list[[j]], tmp_out_file, quote = F, row.names = F, 
  col.names = F)
}
```

## Make VCFs
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

sbatch gen_popVCFs_01_02.sh
sbatch gen_popVCFs_03_06.sh
sbatch gen_popVCFs_07_09.sh

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
#SBATCH -t 24:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/
ATLANTIC_FILE=Atlantic_ref_names.txt
GULF_FILE=Gulf_ref_names.txt
MW_FILE=Midwest_ref_names.txt

SNP_POS_FILE=$SAMP_DIR'hi_fst_SNPS_v1.txt'

for CHR_NUM in {01..02};
  do
  for CHR_LET in K N;
    do
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET;
    #
    SAMP_FILE=$SAMP_DIR$ATLANTIC_FILE;
    OUT_TXT=$OUT_PRE'.Atlantic';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
      --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
    #
    SAMP_FILE=$SAMP_DIR$GULF_FILE;
    OUT_TXT=$OUT_PRE'.Gulf';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
      --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
    #
    SAMP_FILE=$SAMP_DIR$MW_FILE;
    OUT_TXT=$OUT_PRE'.Midwest';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
      --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
    #
    for GRP in {1..5};
      do
      SAMP_FILE=$SAMP_DIR'grp'$GRP'_sampnames.txt';
      OUT_TEXT=$OUT_PRE'.Group_'$GRP;
      vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
        --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
      done;
    done;
  done;

```

## Concatenate chromosome VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

cat Chr01K.Atlantic.recode.vcf > Atlantic.disomic.CDS.hiFst.vcf

awk '!/^#/ {print}' Chr01N.Atlantic.recode.vcf \
>> Atlantic.disomic.CDS.hiFst.vcf

OUT_FILE=Atlantic.disomic.CDS.hiFst.vcf

for CHR_NUM in {02..09};
  do
  for CHR_LET in K N;
    do
    VCF_NAME=Chr$CHR_NUM$CHR_LET'.Atlantic.recode.vcf';
    awk '!/^#/ {print}' $VCF_NAME >> $OUT_FILE;
    done;
  done;

#############

for SAMP_PRE in Gulf Midwest Group_1 Group_2 Group_3 Group_4 Group_5;
  do
  OUT_FILE=$SAMP_PRE'.disomic.CDS.hiFst.vcf'
  cat Chr01K.$SAMP_PRE'.recode.vcf' > $OUT_FILE
  awk '!/^#/ {print}' Chr01N.$SAMP_PRE'.recode.vcf' >> $OUT_FILE
  for CHR_NUM in {02..09};
    do
    for CHR_LET in K N;
      do
      VCF_NAME=Chr$CHR_NUM$CHR_LET'.'$SAMP_PRE'.recode.vcf';
      awk '!/^#/ {print}' $VCF_NAME >> $OUT_FILE;
      done;
    done;
  done;

```
### Compress final files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

for VCF_FILE in *hiFst.vcf;
  do
  gzip $VCF_FILE;
  done

```
### Transfer VCFs to HundsonAlpha
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try

scp *hiFst.vcf.gz grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_vcfs

```


