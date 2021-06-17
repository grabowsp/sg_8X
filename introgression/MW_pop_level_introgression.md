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



