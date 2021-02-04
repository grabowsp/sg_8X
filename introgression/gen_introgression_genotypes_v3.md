# Steps for generating genotypes to be used for HMM looking for introgression

## Overview
* Use Disomic genotypes
* Identify SNPs that have high Fst between pairwise comparisons of 3 source populations
  * Goal is 100k SNPs for each comparison
  * Then light LD pruning using all parent and test samples
* Generate 3 VCFs 
  * One each containing high-Fst SNPs for each pairwise comparions
  * Each containing all parent and all test samples
* Generate txt file
  * Library IDs
  * Group IDs

## Location of files
### Results
* `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3`
### Input files
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

## Choosing SNPs
### Select filtering criteria
* MAF
  * at least in 5 samples = 5/(120*2) = 0.021
  * or for pairwise: 5/(80*2) = 0.03125
* Fst values
  * use 'weir-fst-pop' flag in vcf tools
  * used Fst > 0.5 as cutoff
* Goal is 100k SNPs for each pairwise comparison
### Get Fst values with VCFTools
* Use Fst values from "try_2"
### Generate SNP lists
* Atlantic vs Gulf
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/AtlanticVsGulf_hiFst_SNPs_v3.txt`
  * 67057 SNPs
* Atlantic vs Midwest
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/AtlanticVsMW_hiFst_SNPs_v3.txt`
  * 100,000 SNPs (downsampled from ~156k)
* Gulf vs Midwest
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/GulfVsMW_hiFst_SNPs_v3.txt`
  * 100k SNPs (downsampled from ~150k)
* R script
  * `~/sg_8X/introgression/select_Fst_SNPs_v3.r`

## LD pruning
### Steps
1. Generate total list of samples
2. Use all_samp plink file to prune by LD 
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/GW_all_samps.bed`
### Generate list of samples
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try

cat *names_40.txt *sampnames.txt control_4X_names.txt > all_v3_samp_names.txt

paste all_v3_samp_names.txt all_v3_samp_names.txt > all_v3_plink_names.txt

cp all_v3_plink_names.txt /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

cp all_v3_samp_names.txt /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3
```
### Generate SNP list for plink
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

# Atlantic vs Gulf

SNP_FILE=AtlanticVsGulf_hiFst_SNPs_v3.txt
SNP_OUT=AtlanticVsGulf_hiFst_SNPs_v3.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

# Atlantic vs Midwest

SNP_FILE=AtlanticVsMW_hiFst_SNPs_v3.txt
SNP_OUT=AtlanticVsMW_hiFst_SNPs_v3.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

# Gulf vs Midwest

SNP_FILE=GulfVsMW_hiFst_SNPs_v3.txt
SNP_OUT=GulfVsMW_hiFst_SNPs_v3.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

```
### Generate LD-prune list
* `AtlanticVsGulf_hiFst_SNPs_v3.LD.txt`
  * 40768 SNPs
* `AtlanticVsMW_hiFst_SNPs_v3.LD.txt`
  * 63982 SNPs
* `GulfVsMW_hiFst_SNPs_v3.LD.txt`
  * 63829 SNPs
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/all_v3_plink_names.txt

## Atlantic vs Gulf

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/AtlanticVsGulf_hiFst_SNPs_v3.plink.txt

OUT_TXT_1=AvG_v3_LD_1
OUT_TXT_2=AvG_v3_LD_2
OUT_TXT_3=AtlanticVsGulf_hiFst_SNPs_v3.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR 

cd $OUT_DIR 

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt 
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt 
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

## Atlantic vs Midwest

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/AtlanticVsMW_hiFst_SNPs_v3.plink.txt

OUT_TXT_1=AvM_v3_LD_1
OUT_TXT_2=AvM_v3_LD_2
OUT_TXT_3=AtlanticVsMW_hiFst_SNPs_v3.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt       
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt       
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

## Gulf vs Midwest

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/GulfVsMW_hiFst_SNPs_v3.plink.txt

OUT_TXT_1=GvM_v3_LD_1
OUT_TXT_2=GvM_v3_LD_2
OUT_TXT_3=GulfVsMW_hiFst_SNPs_v3.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3
```

## Make VCFs
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

sbatch gen_v3_VCFs_01_02.sh

sbatch gen_v3_VCFs_03_04.sh
sbatch gen_v3_VCFs_05_06.sh
sbatch gen_v3_VCFs_07_09.sh

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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/all_v3_samp_names.txt

A_V_G_SNPS=$OUT_DIR'AtlanticVsGulf_hiFst_SNPs_v3.LD.txt'
A_V_M_SNPS=$OUT_DIR'AtlanticVsMW_hiFst_SNPs_v3.LD.txt'
G_V_M_SNPS=$OUT_DIR'GulfVsMW_hiFst_SNPs_v3.LD.txt'

for CHR_NUM in {01..02};
  do
  for CHR_LET in K N;
    do
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET;
    #
    OUT_TXT=$OUT_PRE'.AtlanticVsGulf';
    SNP_POS_FILE=A_V_G_SNPS;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
      --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
    #
    OUT_TXT=$OUT_PRE'.AtlanticVsMW';
    SNP_POS_FILE=A_V_M_SNPS;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
      --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
    #
    OUT_TXT=$OUT_PRE'.GulfVsMW';
    SNP_POS_FILE=G_V_M_SNPS;
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
      --keep $SAMP_FILE --recode --recode-INFO --recode-INFO-all;
    #
    done;
  done;

```

## Concatenate chromosome VCFs
* `AtlanticVsGulfhiFst.v3.disomic.CDS.vcf.gz`
  * 40,768 SNPs
* `AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz`
  * 63,982 SNPs
* `GulfVsMWhiFst.v3.disomic.CDS.vcf.gz`
  * 63,829 SNPs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

for SAMP_PRE in AtlanticVsGulf AtlanticVsMW GulfVsMW;
  do
  OUT_FILE=$SAMP_PRE'hiFst.v3.disomic.CDS.vcf'
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
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

for VCF_FILE in *v3.disomic.CDS.vcf;
  do
  gzip $VCF_FILE;
  done

```
### Transfer VCFs to HundsonAlpha
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

scp *vcf.gz grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3

```
### Adjust Permissions on HA
```
cd /home/f2p1/work/grabowsk/data/switchgrass/introgression_v3

chmod 666 *vcf.gz
```

## Make file with LIB_IDs and Groups
* File name:
  * `introgression_samps_and_groups.txt`
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

file_vec <- c(system('ls /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/*names_40.txt', intern = T), 
  system('ls /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/*sampnames.txt', intern = T),
  system('ls /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/*4X_names.txt', intern = T))

samp_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/second_try/'

samp_list <- list()

for(fv in file_vec){
  tmp_group <- sub(samp_dir, '', fv)
  samp_names <- fread(fv, header = F)
  samp_names[ , GROUP := tmp_group]
  samp_list[[tmp_group]] <- samp_names
}

tot_samp_tab <- rbindlist(samp_list)

tot_samp_tab$GROUP <- gsub('_ref_names_40.txt', '', tot_samp_tab$GROUP)
tot_samp_tab$GROUP <- gsub('_sampnames.txt', '', tot_samp_tab$GROUP)
tot_samp_tab$GROUP <- gsub('_names.txt', '', tot_samp_tab$GROUP)

colnames(tot_samp_tab) <- c('LIB', 'GROUP')

out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/introgression_samps_and_groups.txt'

fwrite(tot_samp_tab, out_file, sep = '\t')
```
### Transfor to HA
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

scp /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/introgression_samps_and_groups.txt grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3
```

