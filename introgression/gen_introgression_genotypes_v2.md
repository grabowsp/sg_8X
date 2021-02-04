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
    * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try`
  * Sample files
    * `Atlantic_ref_samps_40.txt`
    * `Gulf_ref_samps_40.txt`
    * `Midwest_ref_samps_40.txt`
    * Original candidate list: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/Introgression_target_individuals.txt`
    * Grp3 split: `Grp4_subgroups.csv`
    * 4X control file: `control_4X_names.txt`

## Generate file with just reference sampe names
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

tail -n +2 Atlantic_ref_samps_40.txt | cut -f 1 > Atlantic_ref_names_40.txt
tail -n +2 Gulf_ref_samps_40.txt | cut -f 1 > Gulf_ref_names_40.txt
tail -n +2 Midwest_ref_samps_40.txt | cut -f 1 > Midwest_ref_names_40.txt
```

## Choosing SNPs
### Select filtering criteria
* MAF
  * at least in 5 samples = 5/(120*2) = 0.021
  * or for pairwise: 5/(80*2) = 0.03125
* Fst values
  * use 'weir-fst-pop' flag in vcf tools
  * want Fst above 0.75
### Get Fst values with VCFTools
#### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

sbatch get_fst_vals_01_03_v2.sh
sbatch get_fst_vals_04_06_v2.sh
sbatch get_fst_vals_07_09_v2.sh
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/
ATLANTIC_FILE=Atlantic_ref_names_40.txt
GULF_FILE=Gulf_ref_names_40.txt
MW_FILE=Midwest_ref_names_40.txt

for CHR_NUM in {01..03};
  do
  for CHR_LET in K N;
    do
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    #
    MAF=0.021;
    OUT_TXT=ALL3_Fst_Chr$CHR_NUM$CHR_LET'_v2';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE \
      --weir-fst-pop $SAMP_DIR$GULF_FILE \
      --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    MAF=0.03125;
    OUT_TXT=GULFvMW_Fst_Chr$CHR_NUM$CHR_LET'_v2';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$GULF_FILE --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    OUT_TXT=ATLANTICvMW_Fst_Chr$CHR_NUM$CHR_LET'_v2';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    OUT_TXT=ATLANTICvGULF_Fst_Chr$CHR_NUM$CHR_LET'_v2';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE \
      --weir-fst-pop $SAMP_DIR$GULF_FILE;
    done;
  done;
```

### Consolidate SNP list
* Final file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/hi_fst_SNPS_v2.txt`
  * 87306 SNPs
* in R
  * `~/sg_8X/introgression/consolidate_Fst_SNPs_v2.r`

## Make Test Sample lists
* in R
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/'

samp_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/first_try/Introgression_target_individuals.txt'

samp_info <- fread(samp_file)

grp3_split_file <- paste(res_dir, 'Grp4_subgroups.csv', sep = '')
grp3_sub_info <- fread(grp3_split_file)

subgroup_names <- unique(grp3_sub_info$Sub_Group)

for(k in seq(length(subgroup_names))){
  tmp_sgn <- subgroup_names[k]
  tmp_totname <- paste('grp3_', k, sep = '')
  tmp_sinds <- which(grp3_sub_info$Sub_Group == tmp_sgn)
  for(m in tmp_sinds){
    tmp_i_ind <- which(samp_info$ID == grp3_sub_info$ID[m])
    samp_info$kmeans_grps_five[tmp_i_ind] <- tmp_totname
  }
}

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

## LD pruning
### Steps
1. Generate total list of samples
2. Use all_samp plink file to prune by LD 
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/GW_all_samps.bed`
### Generate list of samples
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cat *names_40.txt *sampnames.txt > all_v2_samp_names.txt
paste all_v2_samp_names.txt all_v2_samp_names.txt > all_v2_plink_names.txt
```
### Generate SNP list for plink
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cut -f 1 hi_fst_SNPS_v2.txt > tmp_chroms.txt
cut -f 2 hi_fst_SNPS_v2.txt > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > hi_fst_SNPS_v2.plink.txt
```
### Generate LD-prune list
* `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/hi_fst_SNPS_v2.LD.txt`
  * 78815 SNPs
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/all_v2_plink_names.txt

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/hi_fst_SNPS_v2.plink.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 50'kb' 100 0.3 --out introgress_v2_ld0.3

cp introgress_v2_ld0.3.prune.in /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cut -d '_' -f 1 introgress_v2_ld0.3.prune.in > tmp_chroms_in.txt 
cut -d '_' -f 2 introgress_v2_ld0.3.prune.in > tmp_pos_in.txt 
paste tmp_chroms_in.txt tmp_pos_in.txt > hi_fst_SNPS_v2.LD.txt
```

## Make VCFs
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

sbatch gen_v2_popVCFs_01_02.sh
sbatch gen_v2_popVCFs_03_04.sh
sbatch gen_v2_popVCFs_05_06.sh
sbatch gen_v2_popVCFs_07_08.sh
sbatch gen_v2_popVCFs_09.sh
# forgot to include the 4X controls
sbatch gen_v2_4xControlVCFs_all.sh
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/
ATLANTIC_FILE=Atlantic_ref_names_40.txt
GULF_FILE=Gulf_ref_names_40.txt
MW_FILE=Midwest_ref_names_40.txt

SNP_POS_FILE=$SAMP_DIR'hi_fst_SNPS_v2.LD.txt'

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
    for GRP in 1 2 4 5;
      do
      SAMP_FILE=$SAMP_DIR'grp'$GRP'_sampnames.txt';
      OUT_TXT=$OUT_PRE'.Group_'$GRP;
      vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
        --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
      done;
    for GRP in 3_1 3_2 3_3 3_4;
      do
      SAMP_FILE=$SAMP_DIR'grp'$GRP'_sampnames.txt';
      OUT_TXT=$OUT_PRE'.Group_'$GRP;
      vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
        --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
      done;
    done;
  done;

```

## Concatenate chromosome VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cat Chr01K.Atlantic.recode.vcf > Atlantic.disomic.CDS.hiFst.v2.vcf

awk '!/^#/ {print}' Chr01N.Atlantic.recode.vcf \
>> Atlantic.disomic.CDS.hiFst.v2.vcf

OUT_FILE=Atlantic.disomic.CDS.hiFst.v2.vcf

for CHR_NUM in {02..09};
  do
  for CHR_LET in K N;
    do
    VCF_NAME=Chr$CHR_NUM$CHR_LET'.Atlantic.recode.vcf';
    awk '!/^#/ {print}' $VCF_NAME >> $OUT_FILE;
    done;
  done;

#############

for SAMP_PRE in Gulf Midwest Group_1 Group_2 Group_4 Group_5 Group_3_1 Group_3_2 Group_3_3 Group_3_4;
  do
  OUT_FILE=$SAMP_PRE'.disomic.CDS.hiFst.v2.vcf'
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

SAMP_PRE=4Xcontrol

OUT_FILE=$SAMP_PRE'.disomic.CDS.hiFst.v2.vcf'
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
 
***# NEED TO CONCATENATE, COMPRESS AND TRANSFER THE 4X VCF #*** 

```
### Compress final files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

for VCF_FILE in *hiFst.v2.vcf;
  do
  gzip $VCF_FILE;
  done

gzip 4Xcontrol.disomic.CDS.hiFst.v2.vcf
```
### Transfer VCFs to HundsonAlpha
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

scp *hiFst.v2.vcf.gz grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v2

scp 4Xcontrol.disomic.CDS.hiFst.v2.vcf.gz grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v2
```
### Adjust Permissions on HA
```
cd /home/f2p1/work/grabowsk/data/switchgrass/introgression_v2

chmod 666 *vcf.gz
```

