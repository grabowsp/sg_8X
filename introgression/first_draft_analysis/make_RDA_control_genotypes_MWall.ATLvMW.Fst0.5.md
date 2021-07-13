# Generating genotypes that can be used as controls to compare to RDA results
* For comparing to the MWall ATLvMW Fst 0.5 RDA results

## Overview
* Generate genotypes to be used for RDA to compare to hiFst results for MWall
samples using the ATLvMW Fst 0.5 SNPs
* Control SNPs have same chromosome and MAF distribution as the hiFst-RDA SNPs

## 1. Identify SNPs used for hiFst RDA
* This file has the information about the hiFst-RDA SNPs
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_ATL_into_MW_Fst_0.5.RDA_SNP_info.txt`
 * This includes the MAF for each SNP
* hiFst-RDA genotypes and SNP info  generated within this R script:
 * `~/sg_8X/introgression/first_draft_analysis/MWall_MWvATL_Fst0.5_indiv_introgression_firstdraft.r`

## 2. Get MAF for hiFst-RDA SNPs
* Contained in the hiFst-RDA SNP info file:
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_ATL_into_MW_Fst_0.5.RDA_SNP_info.txt`

## 3. Get allele frequencies of all CDS SNPs for MWall
### Location of Chromosome allele frequency files
* Look for: `Chr*MWall.frq`
* on NERSC
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft`
* on HAIB
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files`
### Generate Sample Name File to be used for VCFtools
* at NERSC
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

library(data.table)

res_file <- '/global/homes/g/grabowsp/data/switchgrass/results_tables_8X/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_file)

mw_names <- res_tab[grep('MW_', res_tab$subgrp_v2), list(samp_name)]

mw_name_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt'

fwrite(mw_names, mw_name_outfile, col.names = F, row.names = F)
```
### Get allele frequencies for MWall
* I ended up running interactively to make sure finished before Cori shutdown
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

sbatch get_MWall_CDS_freqs_01_03.sh
```
### Example Script
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt

for CHR_NUM in {01..03};
  do
  for CHR_LET in K N;
    do
    #
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET;
    #
    OUT_TXT=$OUT_PRE'.MWall';
    #
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --keep $SAMP_FILE --freq;
    done;
  done;
```
## Package and transfer to HAIB
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/
#ls *MWall.frq
tar -czvf MWall_CDS_freqs.tar.gz *MWall.frq
scp MWall_CDS_freqs.tar.gz grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files

# on HAIB

cd /home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files
tar -xzvf MWall_CDS_freqs.tar.gz
```

## 4. Select control sets for MWall
* Location of files on HAIB
 * `/home/f1p1/tmp/switchgrass_8X/MWall_control_snps/`
* For MW vs ATL, Fst 0.5
 * `ATLvMW_Fst_0.5_MWall_RDA_control_snps.REPNUMBER.txt`
 * generated with: `select_control_SNPs_MWall.ATLvMW.Fst0.5.r`

## 5. Generate MWall Control VCFs
### File Locations
* First control VCF
 * `/home/f1p1/tmp/switchgrass_8X/MWall_control_vcfs/ATLvMW_Fst_0.5_MWall_RDAcontrol.tot.1.vcf`
### Commands used to generate VCF
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f1p1/tmp/switchgrass_8X/MWall_control_vcfs/
cd $OUT_DIR

VCF_DIR=/home/f1p1/tmp/switchgrass_8X/MWall_big_vcfs/
VCF_SUF=.MWall.big.recode.vcf.gz

REP_NUM=1

SNP_FILE=/home/f1p1/tmp/switchgrass_8X/MWall_control_snps/ATLvMW_Fst_0.5_MWall_RDA_control_snps.$REP_NUM'.txt'

OUT_SUF=ATLvMW_Fst_0.5_MWall_RDAcontrol
TOT_OUT=$OUT_SUF'.tot.'$REP_NUM'.vcf'

for CHR_NUM in {01..09};
  do
  for CHR_LET in K N;
    do
    #
    VCF_FULL=$VCF_DIR'Chr'$CHR_NUM$CHR_LET$VCF_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET;
    #
    OUT_TXT=$OUT_PRE'.'$OUT_SUF'.'$REP_NUM;
    #
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_FILE \
    --recode --recode-INFO --recode-INFO-all;
    #
    done;
  done;

cat Chr01K.$OUT_SUF'.'$REP_NUM'.recode.vcf' > $TOT_OUT
awk '!/^#/ {print}' Chr01N.$OUT_SUF'.'$REP_NUM'.recode.vcf' >> $TOT_OUT
for CHR_NUM in {02..09};
 do
 for CHR_LET in K N;
  do
  VCF_NAME=Chr$CHR_NUM$CHR_LET'.'$OUT_SUF'.'$REP_NUM'.recode.vcf';
  awk '!/^#/ {print}' $VCF_NAME >> $TOT_OUT;
  done;
 done;

rm Chr*$OUT_SUF'.'$REP_NUM'.recode.vcf'
rm Chr*$OUT_SUF'.'$REP_NUM'.log'
```

## 6. Generate Control RDA genotypes
* Script used to generate genotype file
 * NOTE: Need to adjust to correct missing data
 * `~/sg_8X/introgression/first_draft_analysis/MWall_MWvATL_Fst0.5_make_RDAcontrol_genos.r`
* Location of control RDA genotypes
 * `/home/f1p1/tmp/switchgrass_8X/MWall_control_vcfs/ATLvMW_Fst_0.5_MWall_RDAcontrol.genos.1.rds`

## Extra stuff
### Get MAF for the hiFST SNPs
* Chose to run interactively to make sure finished before Cori shutdown
 * Technically didn't need to run this because could have extracted from
the full SNP set, but this made it easier
```
module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files/

for COMPS in AtlanticvsMidwest GulfvsMidwest;
 do
 for FST in 0.5 1.0;
  do
  VCF_NAME=$VCF_DIR$COMPS'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.vcf.gz';
  OUT_TXT=$COMPS'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.MWall';
  #
  vcftools --gzvcf $VCF_NAME --out $OUT_TXT --keep $SAMP_FILE --freq;
  done;
 done;

```


