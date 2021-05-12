# Steps for examining introgression in MW 8X samples

## Overview
* Select hi-Fst SNPs for MW_v_GULF and MW_v_ATL
  * done previously
* Assign allele-state for alleles at hi-Fst SNPs
* Individual-level analysis
  * Assign genotype-states to each MW-8X individual
  * Window analysis looking for elevated introgression levels
  * Plot window values
* Population-level analysis, for pops with sufficient population sizes
  * Look for SNPs with allele freqency shifts away from MW
  * Use window analysis and magintude of shift away from MW to select candidate SNPs
  * Plot window analysis results

## Location of files (at HAIB)
* Output directory
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis`
* Directory with many input files
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/`
* Hi-Fst VCFs
  * File generation outlined in:
    * `/home/grabowsky/tools/workflows/sg_8X/introgression/gen_introgression_genotypes_v3.md`
  * Atlantic vs Gulf hi-Fst VCF
    * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsGulfhiFst.v3.disomic.CDS.vcf.gz`
  * Gulf vs Midwest
    * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMWhiFst.v3.disomic.CDS.vcf.gz`
  * Atlantic vs Midwest
    * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz`
* Original Sample and Group Map
  * Contains training group samples and original 8X splits
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/introgression_samps_and_groups.txt`

## Get info to calculate genotype frequencies for training sets
### Generate sample files
* done previously in R
* Output files:
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/Atlantic_libnames.txt`
* R code:
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
samp_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/introgression_samps_and_groups.txt'

samp_info <- fread(samp_file)

### SET OUTPUT ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'

for(pop in unique(samp_info$GROUP)){
  tmp_names <- samp_info[GROUP == pop, list(LIB)]
  tmp_outname <- paste(out_dir, pop, '_libnames.txt', sep = '')
  fwrite(tmp_names, file = tmp_outname, col.names = F, sep = '\t')
}
```
### Extract genotype frequency info from training sets with VCFtools
#### MW vs GULF
* Done previously
```
bash
source activate bioinformatics_env

OUTDIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
VCF_SHORT=GulfVsMWhiFst.v3.disomic.CDS.vcf.gz
GROUP_SHORT=Midwest_libnames.txt
OUT_PRE=MW_training_MWvGULF

cd $OUTDIR

VCF_IN=$OUTDIR$VCF_SHORT
GROUP_IN=$OUTDIR$GROUP_SHORT
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --hardy

GROUP_SHORT=Gulf_libnames.txt
OUT_PRE=GULF_training_MWvGULF
VCF_IN=$OUTDIR$VCF_SHORT
GROUP_IN=$OUTDIR$GROUP_SHORT
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --hardy
```
#### MW vs ATL
```
bash
source activate bioinformatics_env

OUTDIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
VCF_SHORT=AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz
GROUP_SHORT=Midwest_libnames.txt
OUT_PRE=MW_training_MWvATL

cd $OUTDIR

VCF_IN=$OUTDIR$VCF_SHORT
GROUP_IN=$OUTDIR$GROUP_SHORT
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --hardy

GROUP_SHORT=Atlantic_libnames.txt
OUT_PRE=ATL_training_MWvATL
VCF_IN=$OUTDIR$VCF_SHORT
GROUP_IN=$OUTDIR$GROUP_SHORT
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --hardy
```

## Assign allele states
* scripts adapted from  `~/sg_8X/introgression/set_introgression_allele_states.r`
  * decided to hard-code because wanted to test different cutoffs so that all hi-Fst SNPs get included
### MW vs GULF
* R script
  * `~/sg_8X/introgression/MW_v_GULF_allele_states_v2.r`

* Output files:
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_allele_states.txt`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_keep_pos.txt`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_ref_freq.txt`

* Important info:
  * 58,829 SNPs pass filtering (missing data < 20%)
  * Chose 2.5 as probability allele ratio
    * all hi-Fst SNPs have at least 1 informative allele
  * 38,217 SNPs with alleles informative about MW ancestry
  * 40,885 SNPs with alleles informative about GULF ancestry
  * 20,273 SNPs have both alleles as informative

### MW vs ATL
* R script
  * `~/sg_8X/introgression/MW_v_ATL_allele_states.r`

* Output files:
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_allele_states.txt`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_keep_pos.txt`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_ref_freq.txt`

* Important info:
  * 59,113 SNPs pass filtering (missing data < 20%)
  * Chose 2.5 as probability allele ratio
    * all hi-Fst SNPs have at least 1 informative allele
  * 42,195 SNPs with alleles informative about MW ancestry
  * 37,1443 SNPs with alleles informative about ATL ancestry
  * 20,525 SNPs have both alleles as informative

## Generate VCFs for MW-8X samples
### Make sample name list for all MW-8X
* Output
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_names.txt`
* File with sample subgroups and ploidy
  * `/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.0.txt`
* R script:
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
info_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.0.txt'
samp_info <- fread(info_file)

meta_file <- '/home/grabowsky/tools/workflows/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(meta_file)

### SET OUTPUT ###
mw_8X_name_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/MW_8X_names.txt', sep = '')

mw_all_name_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/MW_all_names.txt', sep = '')

#####
mw_8X_inds <- intersect(which(samp_info$ploidy == '8X'), 
  grep('MW_', samp_info$sub_grp))

mw_8X_names <- samp_info[mw_8X_inds, list(samp_name)]
fwrite(mw_8X_names, file = mw_8X_name_out, col.names = F, sep = '\t')

mw_inds <- grep('MW_', samp_info$sub_grp)
mw_names <- samp_info[mw_inds, list(samp_name)]

fwrite(mw_names, file = mw_all_name_out, col.names = F, sep = '\t')
```
### VCF for MW8X - MW vs GULF
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/
VCF_IN=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMWhiFst.v3.disomic.CDS.vcf.gz
SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_names.txt
POS_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_keep_pos.txt
OUTFILE=MW8X_MWvGULF.vcf

cd $OUT_DIR

vcftools --gzvcf $VCF_IN -c  --keep $SAMP_FILE \
--positions $POS_FILE --recode> $OUTFILE
```
### VCF for MW8X - MW vs ATL
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/
VCF_IN=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz
SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_names.txt
POS_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_keep_pos.txt
OUTFILE=MW8X_MWvATL.vcf

cd $OUT_DIR

vcftools --gzvcf $VCF_IN -c  --keep $SAMP_FILE \
--positions $POS_FILE --recode> $OUTFILE
```
### VCFs for all MW
* generate VCFs of all samples for MWvsGULF and MWvsATL SNPs on NERSC
#### Transfer necessary info to NERSC 
* Only for MWvsGULF and MWvsATL at this point
  * will have to do GULFvsATL at another time, if choose to do that
```
# in NERSC
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3
scp grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_ATL_keep_pos.txt .
scp grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_v_GULF_keep_pos.txt .
```
#### Generate subVCFs on NERSC
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

sbatch gen_v3_allsamps_VCFs_01_02.sh
sbatch gen_v3_allsamps_VCFs_03_06.sh
sbatch gen_v3_allsamps_VCFs_07_09.sh
```
#### Concatenate VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

for SAMP_PRE in AtlanticVsMW_allsamps GulfVsMW_allsamps;
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
#### Compress final files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3
for VCF_FILE in *allsampshiFst.v3.disomic.CDS.vcf;
  do
  gzip $VCF_FILE;
  done
```
#### Transfer to HAIB
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

scp *_allsampshiFst.v3.disomic.CDS.vcf.gz grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3
```
#### Generate VCFs of MW_all for MWvsGULF
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/
VCF_IN=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMW_allsampshiFst.v3.disomic.CDS.vcf.gz
SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_all_names.txt
OUTFILE=MWall_MWvGULF.vcf

cd $OUT_DIR

vcftools --gzvcf $VCF_IN -c  --keep $SAMP_FILE \
--recode> $OUTFILE
```
### VCF for all MW - MW vs ATL
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/
VCF_IN=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMW_allsampshiFst.v3.disomic.CDS.vcf.gz
SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_all_names.txt
OUTFILE=MWall_MWvATL.vcf

cd $OUT_DIR

vcftools --gzvcf $VCF_IN -c  --keep $SAMP_FILE \
--recode> $OUTFILE
```

## MW-8X individual introgression analysis
* all MW-8X
* based on `~`
* had to hard-code because my Anaconda environment doesn't have 'Rscript' (why????!!!!) and I didn't want to mess it up by installing the 'base' package

### MW-8X Gulf into MW background
* R script
  * `~/sg_8X/introgression/MW8X_MWvGULF_indiv_introgression.r`

* Outputs
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_GULF_into_MW_10SNP_window_results.rds`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_GULF_into_MW_10SNP_window.adj_n_pop2.pdf`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_GULF_into_MW_10SNP_window.pop2_score.pdf`

### MW-8X ATL into MW background
* R script
  * `~/sg_8X/introgression/MW8X_MWvATL_indiv_introgression.r`
* Outputs
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_ATL_into_MW_10SNP_window_results.rds`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_ATL_into_MW_10SNP_window.adj_n_pop2.pdf`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MW_8X_ATL_into_MW_10SNP_window.pop2_score.pdf`

### MW-all Gulf into MW background
* R script
  * `~/sg_8X/introgression/MWall_MWvGULF_indiv_introgression.r`
* Outputs
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_GULF_into_MW_10SNP_window_results.rds`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_GULF_into_MW_10SNP_window.adj_n_pop2.pdf`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_GULF_into_MW_10SNP_window.pop2_score.pdf`

### MW-all ATL into MW background
* R script
  * `MWall_MWvATL_indiv_introgression.r`
* Outputs
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_ATL_into_MW_10SNP_window_results.rds`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_ATL_into_MW_10SNP_window.adj_n_pop2.pdf`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_ATL_into_MW_10SNP_window.pop2_score.pdf`


## Generate matrices formatted for RDA analysis
* R script used to generate matrices for all MW-8X, MW_01 (MW8X-West), and MW_03 (MW8X-East)
  * `~/sg_8X/introgression/MW8X_gen_RDA_files.r`
### Generate RDA genotypes for all MW samples
* R Script
  * `~/sg_8X/introgression/MWall_gen_RDA_files.r`
* Outputs
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_GULF_into_MW.RDA_genos.rds`
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_analysis/MWall_ATL_into_MW.RDA_genos.rds`


### NEXT
* Plots ordered by subgroups



