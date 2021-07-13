# Steps for calculationg the F_val values for MW8X

## Steps
1. Calculate allele frequencies for training sets
 * This is usually done when calling allele states
2. Calculate HWE for MW8X
 * Used for calculating allele frequency in R script
3. Run R script to generate F_vals

## Location of files
* File with training set allele frequencies
 * Output from script used to call allele-states
  * `/home/f1p1/tmp/switchgrass_8X/introgression_files/MW_v_ATL_Fst0.5_ref_freq.txt`
* MW8X sample names
  * `/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/mw_introgression_sets/MW8X_all_names.txt`
* SNP file
 * Output from script used to call Allele-states
  * `/home/f1p1/tmp/switchgrass_8X/introgression_files/MW_v_ATL_Fst0.5_keep_pos.txt`
* VCF to use
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_vcfs/AtlanticvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`

## Location of Output Files
* MW8X MWvATL Fst 0.5 HWE file
 `/home/f1p1/tmp/switchgrass_8X/firstdraft_hwe_files/MWvATL_Fst0.5_files/MWvATL_Fst0.5_MW8X.hwe`
* MW8X MWvATL Fst 0.5 ATL into MW Fvals
 `/home/f1p1/tmp/switchgrass_8X/firstdraft_fval_files/MW8X_Fval_ATL_into_MW_allsnps_Fst0.5.txt`

## Calculate HWE for sample sets
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f1p1/tmp/switchgrass_8X/firstdraft_hwe_files/MWvATL_Fst0.5_files/

cd $OUT_DIR

VCF_IN=/home/f1p1/tmp/switchgrass_8X/firstdraft_vcfs/AtlanticvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz

POS_FILE=/home/f1p1/tmp/switchgrass_8X/introgression_files/MW_v_ATL_Fst0.5_keep_pos.txt

# MW8X samples

MW8X_IN=/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/mw_introgression_sets/MW8X_all_names.txt
MW8X_OUT=MWvATL_Fst0.5_MW8X

vcftools --gzvcf $VCF_IN --out $MW8X_OUT --keep $MW8X_IN --positions $POS_FILE \
--hardy
```

## Run R script to generate Fvals
```
bash
source activate R_analysis

RSCRIPT_FILE=/home/grabowsky/tools/workflows/sg_8X/introgression/pop_introgression_v3.r

MW8X_HWE=/home/f1p1/tmp/switchgrass_8X/firstdraft_hwe_files/MWvATL_Fst0.5_files/MWvATL_Fst0.5_MW8X.hwe

TRAIN_FREQ=/home/f1p1/tmp/switchgrass_8X/introgression_files/MW_v_ATL_Fst0.5_ref_freq.txt

ALLELE_STATES=/home/f1p1/tmp/switchgrass_8X/introgression_files/MW_v_ATL_Fst0.5_allele_states.txt

WINDOW_SIZE=10

WINDOW_CUT=5

OUT_DIR=/home/f1p1/tmp/switchgrass_8X/firstdraft_fval_files/

SUBGRP_NAME=MW8X
TRAIN_1_NAME=MW
TRAIN_2_NAME=ATL

SD_CUT=3

cd $OUT_DIR

Rscript $RSCRIPT_FILE $MW8X_HWE $TRAIN_FREQ $ALLELE_STATES $WINDOW_SIZE \
$WINDOW_CUT $OUT_DIR $SUBGRP_NAME $TRAIN_1_NAME $TRAIN_2_NAME $SD_CUT

mv MW8X_Fval_ATL_into_MW_allsnps.txt MW8X_Fval_ATL_into_MW_allsnps_Fst0.5.txt

```



