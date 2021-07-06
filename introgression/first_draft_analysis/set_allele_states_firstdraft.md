# Steps for setting the allele states for first draft analysis

## Steps
1. Calculate allele freq and HW frequency in training sets
2. Use R script to set allele states

## Generate HW and allele freq files from vcfTools
```
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

VCF_FILE=

cd $OUT_DIR

TRAIN_SET=GULF
FST=0.5
COMP=AtlanticVsGulf

OUT_PRE=$TRAIN_SET'_train_'$COMP'_firstdraft_Fst'$FST
SAMP_FILE=$TRAIN_SET'_firstdraft_train_filt40.txt'
SNP_FILE=$COMP'_hiFst_SNPs_firstdraft_Fst'$FST'.LD.MISS.snps.txt'

vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq




```
