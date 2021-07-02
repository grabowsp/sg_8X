# Steps for Population-level analysis of introgression

## Overview
* use allele frequencies to identify SNPs introgressed at high frequency
* Look at:
  * All MW8X
  * MW8X without MW_01_hi
  * MW_01_hi
  * MW_01
  * MW_03
  * MW4X_introgression
* Compare them to both Gulf and Atlantic
* Steps for each comparison
  * Generate sample name files
  * Generate HW file for test samples using vcfTools
  * run pop_introgression_v1.r

## Generate sample name files
* R script used to generate name files:
  * `~/sg_8X/sample_sets/MW_introgression_sets.r`
* File names
  * Directory (on HAIB)
    * `/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/mw_introgression_sets/`
  * MW8X - all
    * `MW8X_all_names.txt`
  * MW8X - regular (no hi introgression)
    * `MW8X_reg_names.txt`
  * MW8X - hi introgression
    * `MW8X_hi_names.txt`
  * MW8X - West (MW_01; no hi introgression)
    * `MW_01_names.txt`
  * MW8X - East
    * `MW_03_names.txt`
  * MW4X - with high introgression levels
    * `MW4X_intro_names.txt`

## Generate HW files for test samples
* Perform for both GULF and ATL comparisons
```
bash
source activate bioinformatics_env

OUTDIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
VCF_SHORT=GulfVsMWhiFst.v3.disomic.CDS.vcf.gz
VCF_IN=$OUTDIR$VCF_SHORT

POS_SHORT=MW_v_GULF_keep_pos.txt
POS_FILE=$OUTDIR'MW_analysis/'$POS_SHORT

GROUP_DIR=/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/mw_introgression_sets/

cd $OUTDIR

for GN in MW8X_all MW8X_reg MW8X_hi MW_01 MW_03 MW4X_intro;
  do
  GROUP_SHORT=$GN'_names.txt';
  GROUP_IN=$GROUP_DIR$GROUP_SHORT;
  OUT_PRE=$GN'_MWvGULF';
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN \
--positions $POS_FILE --hardy;
  done;

# MW vs ATL
VCF_SHORT=AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz
VCF_IN=$OUTDIR$VCF_SHORT

POS_SHORT=MW_v_ATL_keep_pos.txt
POS_FILE=$OUTDIR'MW_analysis/'$POS_SHORT

for GN in MW8X_all MW8X_reg MW8X_hi MW_01 MW_03 MW4X_intro;
  do
  GROUP_SHORT=$GN'_names.txt';
  GROUP_IN=$GROUP_DIR$GROUP_SHORT;
  OUT_PRE=$GN'_MWvATL';
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN \
--positions $POS_FILE --hardy;
  done;
```

## Pull out outlier SNPs

### MW v GULF
* Directory: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/`
* Files explantaions:
  * `GROUPNAME_Fval_GULF_into_MW_10SNP_allwindows.txt`
    * 'F-value' results for all 10-SNP windows
  * `GROUPNAME_Fval_GULF_into_MW_10SNP_topwindows.txt`
    * 'F-value' results for 10-SNP windows above an SD-based threshold
  * `GROUPNAME_Fval_GULF_into_MW_allsnps.txt`
    * 'F-value' results for all individual SNPs
  * `GROUPNAME_Fval_GULF_into_MW_outliersnps.txt`
    * 'F-value' results for SNPs with F-values above a SD-based threshold
  * `GROUPNAME_all_Fval_GULF_into_MW_topwindowsnps.txt`
    * 'F-value' results for SNPs chosen from top 10-SNP windows
  * `GROUPNAME_Fval_GULF_into_MW_10SNP_windows.pdf`
    * Plot of F-value windows across the genome
```
bash
source activate R_analysis

OUTDIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/
cd $OUTDIR

R CMD BATCH /home/grabowsky/tools/workflows/sg_8X/introgression/MW8X_all_pop_MWvGULF.r

# NEXT: make files for other samples sets and for comparing to ATL

```
#### MW8X_reg
```
bash
source activate R_analysis

R_SCRIPT=/home/grabowsky/tools/workflows/sg_8X/introgression/pop_introgression_v2.r

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/

cd $OUT_DIR

# SUB_HW=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW8X_reg_MWvGULF.hwe
TRAIN_FREQ=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_ref_freq.txt
ALLELE_STATE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_allele_states.txt

WIN_SIZE=10
WIN_CUT=5
#OUT_DIR
#SUB_NAME=MW8X_reg
TRAIN_POP1_NAME=MW
TRAIN_POP2_NAME=GULF

# Rscript $R_SCRIPT $SUB_HW $TRAIN_FREQ $ALLELE_STATE $WIN_SIZE $WIN_CUT \
# $OUT_DIR $SUB_NAME $TRAIN_POP1_NAME $TRAIN_POP2_NAME

for GN in MW8X_all MW8X_reg MW4X_intro MW8X_hi MW_01 MW_03;
  do
  SUB_HW=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/$GN'_MWvGULF.hwe';
  SUB_NAME=$GN;
  Rscript $R_SCRIPT $SUB_HW $TRAIN_FREQ $ALLELE_STATE $WIN_SIZE $WIN_CUT \
$OUT_DIR $SUB_NAME $TRAIN_POP1_NAME $TRAIN_POP2_NAME;
  done;

TRAIN_FREQ=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_ref_freq.txt
ALLELE_STATE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_allele_states.txt

TRAIN_POP2_NAME=ATL

for GN in MW8X_all MW8X_reg MW4X_intro MW8X_hi MW_01 MW_03;
  do
  SUB_HW=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/$GN'_MWvATL.hwe';
  SUB_NAME=$GN;
  Rscript $R_SCRIPT $SUB_HW $TRAIN_FREQ $ALLELE_STATE $WIN_SIZE $WIN_CUT \
$OUT_DIR $SUB_NAME $TRAIN_POP1_NAME $TRAIN_POP2_NAME;
  done;

```

## Analysis of patterns of outlier SNPs
* R script
  * `~/sg_8X/introgression/MW_pop_intro_analysis.r`

#### Using a stricter SD cutoff
```
bash
source activate R_analysis

R_SCRIPT=/home/grabowsky/tools/workflows/sg_8X/introgression/pop_introgression_v3.r

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/

cd $OUT_DIR

# SUB_HW=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW8X_reg_MWvGULF.hwe
TRAIN_FREQ=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_ref_freq.txt
ALLELE_STATE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_allele_states.txt

WIN_SIZE=10
WIN_CUT=5
#OUT_DIR
#SUB_NAME=MW8X_reg
TRAIN_POP1_NAME=MW
TRAIN_POP2_NAME=GULF
SD_CUT=4

for GN in MW8X_all MW8X_reg MW4X_intro MW8X_hi MW_01 MW_03;
  do
  SUB_HW=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/$GN'_MWvGULF.hwe';
  SUB_NAME=$GN'_strict';
  Rscript $R_SCRIPT $SUB_HW $TRAIN_FREQ $ALLELE_STATE $WIN_SIZE $WIN_CUT \
$OUT_DIR $SUB_NAME $TRAIN_POP1_NAME $TRAIN_POP2_NAME $SD_CUT;
  done;

TRAIN_FREQ=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_ref_freq.txt
ALLELE_STATE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_allele_states.txt

TRAIN_POP2_NAME=ATL

for GN in MW8X_all MW8X_reg MW4X_intro MW8X_hi MW_01 MW_03;
  do
  SUB_HW=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/$GN'_MWvATL.hwe';
  SUB_NAME=$GN'_strict';
  Rscript $R_SCRIPT $SUB_HW $TRAIN_FREQ $ALLELE_STATE $WIN_SIZE $WIN_CUT \
$OUT_DIR $SUB_NAME $TRAIN_POP1_NAME $TRAIN_POP2_NAME $SD_CUT;
  done;



```

