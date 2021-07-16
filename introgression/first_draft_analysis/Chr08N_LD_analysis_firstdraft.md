# LD analysis of Fval/RDA overlap on Chr08N

## Steps
* r^2 calculations using plink
* LDheatmap visualizations in R
* Overall analysis in R

## Calculating r^2 in plink
* LDheatmap needs a square matrix
```
--r2 square gz
``` 
* I don't think I need the `--ldwindow-r2 0.05` flag if requesting a square matrix
* I should include a MAF filter of 0.05
```
--maf 0.05
```

* I will do a progression:
 * RDA snps
 * hiFst snps
 * random set of SNPs

## Analysis plan
* RDA snps and hiFst snps
 * all natv2 samples
 * all MW samples
 * all Gulf samples
 * all Atlantic samples
* Random set of SNPs
 * MW all
 * MW 8X
 * MW 4X
 * MW introgression
 * MW no introgression

## File Locations
### SNPs in plink format
* Chr08N RDA-only SNPs
 * NERSC:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5_RDASNPs.plink.txt`
* Chr08N all Fval SNPs
 * NERSC:
  `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5_FvalSNPs.plink.txt`
* RDA-only SNPs - chromosome-wide
 * HAIB:
  * `/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvsMW_Fst0.5_RDASNPs.plink.txt`
 * NERSC
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/ATLvsMW_Fst0.5_RDASNPs.plink.txt`
* All Fval SNPS 0 chromosome-wide
 * HAIB
  * `/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvsMW_Fst0.5_FvalSNPs.plink.txt`
 * NERSC
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/ATLvsMW_Fst0.5_FvalSNPs.plink.txt`

### Samples in plink format
* all natv2 samples
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/natv2_samps.plink.txt`
* MW samples
* Gulf samples
* Atlantic samples
* MW8X samples
* MW low introgression samples
* MW4X samples

## Make SNP files for plink
* RDA-only SNPs
 * HAIB:
  * `/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvsMW_Fst0.5_RDASNPs.plink.txt`
 * NERSC
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/ATLvsMW_Fst0.5_RDASNPs.plink.txt`
* All Fval SNPS
 * HAIB
  * `/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvsMW_Fst0.5_FvalSNPs.plink.txt`
 * NERSC
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/ATLvsMW_Fst0.5_FvalSNPs.plink.txt`
* Chr08N RDA-only SNPs
 * NERSC:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5_RDASNPs.plink.txt`
* Chr08N all Fval SNPs
 * NERSC:
  `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5_FvalSNPs.plink.txt`
### Generate Files
```
# on HAIB
# bash
# source activate

library(data.table)

file_dir <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_fval_files/'

fval_file <- paste(file_dir, 'MW8X_Fval_ATL_into_MW_allsnps_Fst0.5.txt',
  sep = '')
fvals <- fread(fval_file)

rda_snp_file <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_ATL_into_MW_Fst_0.5.RDA_SNP_info.txt'
rda_snps <- fread(rda_snp_file)

out_dir <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/'

allsnp_out <- paste(out_dir, 'ATLvsMW_Fst0.5_FvalSNPs.plink.txt', sep = '')

rda_snp_out <- paste(out_dir, 'ATLvsMW_Fst0.5_RDASNPs.plink.txt', sep = '')

allsnp_tab <- data.table(snpname = paste(fvals$CHR, fvals$POS, sep = '_'))

rdasnp_tab <- data.table(snpname = rda_snps$ID)

fwrite(allsnp_tab, file = allsnp_out, col.names = F, sep = '\t')
fwrite(rdasnp_tab, file = rda_snp_out, col.names = F, sep = '\t')
```
### Transfer Files to NERSC
```
# on NERSC
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files

scp grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvsMW_Fst0.5_RDASNPs.plink.txt .

scp grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvsMW_Fst0.5_FvalSNPs.plink.txt . 
```
### Generate Chr08N SNP file
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files

grep Chr08N ATLvsMW_Fst0.5_FvalSNPs.plink.txt > \
Chr08N_ATLvsMW_Fst0.5_FvalSNPs.plink.txt

grep Chr08N ATLvsMW_Fst0.5_RDASNPs.plink.txt > \
Chr08N_ATLvsMW_Fst0.5_RDASNPs.plink.txt
```

## Make RDA control SNP file
* Chr08N RDA control
 * NERSC
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt`
* RDA control - genome-wide
 * On HAIB:
  * `/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt`
 * On NERSC:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt`
### Make SNP name file
```
SNP_IN=/home/f1p1/tmp/switchgrass_8X/MWall_control_snps/ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.txt

OUT_DIR=/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/

cd $OUT_DIR

OUT_FILE=$OUT_DIR'ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt'

cut -f 1 $SNP_IN > tmp_chroms.txt
cut -f 2 $SNP_IN > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $OUT_FILE
rm tmp_chroms.txt
rm tmp_pos.txt
```
### Transfer to NERSC
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files

scp grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt .
```
### Make Chr08N file
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files

grep Chr08N ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt > \
Chr08N_ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt
```

## Make Sample Files for plink
* Need to generate:
 * MW samples
 * Gulf samples
 * Atlantic samples
 * MW8X samples
 * MW low introgression samples
 * MW4X samples

### Make MW4X names for vcftools
```
bash
source activate R_analysis

library(data.table)

samp_info_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.1.txt'
samp_info <- fread(samp_info_file)

mw4X_inds <- intersect(grep('MW_', samp_info$subgrp_v2), 
  which(samp_info$ploidy == '4X'))

mw4X_names <- samp_info[mw4X_inds, list(samp_name)]

mw4X_out <- '/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_4X_names.txt'

fwrite(mw4X_names, mw4X_out, col.names = F)
```
### Make name files for plink
```
bash

IN_DIR=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/

OUT_DIR=/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/

cd $OUT_DIR

MW_ALL_IN=$IN_DIR'MW_all_names.txt'
MW_ALL_OUT=$OUT_DIR'MW_all_names.plink.txt'
paste -d '\t' $MW_ALL_IN $MW_ALL_IN > $MW_ALL_OUT

GULF_IN=$IN_DIR'GULF_all_names.txt'
GULF_OUT=$OUT_DIR'GULF_all_names.plink.txt'
paste -d '\t' $GULF_IN $GULF_IN > $GULF_OUT

ATL_IN=$IN_DIR'ATL_all_names.txt'
ATL_OUT=$OUT_DIR'ATL_all_names.plink.txt'
paste -d '\t' $ATL_IN $ATL_IN > $ATL_OUT

MW8X_IN=$IN_DIR'MW8X_introgressed_names.txt'
MW8X_OUT=$OUT_DIR'MW8X_introgressed_names.plink.txt'
paste -d '\t' $MW8X_IN $MW8X_IN > $MW8X_OUT

MW_NO_IN=$IN_DIR'MW_no_introgress_names.txt'
MW_NO_OUT=$OUT_DIR'MW_no_introgress_names.plink.txt'
paste -d '\t' $MW_NO_IN $MW_NO_IN > $MW_NO_OUT

MW4X_IN=$IN_DIR'MW_4X_names.txt'
MW4X_OUT=$OUT_DIR'MW_4X_names.plink.txt'
paste -d '\t' $MW4X_IN $MW4X_IN > $MW4X_OUT
```
### Transfer files to NERSC
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files

scp grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/*names.plink.txt .
```

## Calculate r^2 of RDA SNPs in plink
```
# on Cori
module load python/3.7-anaconda-2019.10
source activate plink_1_env

START_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

cd $START_DIR

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5/

FILE_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/

SNP_FILE=$FILE_DIR'Chr08N_ATLvsMW_Fst0.5_RDASNPs.plink.txt'

MAF_CUT=0.05

# all samples

ALL_SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/natv2_samps.plink.txt

ALL_SAMP_OUT=natv2_Chr08N_ATLvsMW_Fst0.5_RDAonly

plink --bfile GW_all_samps --keep $ALL_SAMP_FILE --extract $SNP_FILE \
--maf $MAF_CUT --r2 square --out $ALL_SAMP_OUT

mv $ALL_SAMP_OUT'.ld' $OUT_DIR 

# all Midwest samples

cd $START_DIR

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/

MW_ALL_SAMP_FILE=$SAMP_DIR'MW_all_names.plink.txt'

MW_ALL_OUT=MWall_Chr08N_ATLvsMW_Fst0.5_RDAonly

plink --bfile GW_all_samps --keep $MW_ALL_SAMP_FILE --extract $SNP_FILE \
--maf $MAF_CUT --r2 square --out $MW_ALL_OUT

mv $MW_ALL_OUT'.ld' $OUT_DIR

# Gulf samples

cd $START_DIR

GULF_SAMP_FILE=$SAMP_DIR'GULF_all_names.plink.txt'
GULF_OUT=GULF_Chr08N_ATLvsMW_Fst0.5_RDAonly

plink --bfile GW_all_samps --keep $GULF_SAMP_FILE --extract $SNP_FILE \
--r2 square --out $GULF_OUT

mv $GULF_OUT'.ld' $OUT_DIR

# Atlantic samples

cd $START_DIR

ATL_SAMP_FILE=$SAMP_DIR'ATL_all_names.plink.txt'
ATL_OUT=ATL_Chr08N_ATLvsMW_Fst0.5_RDAonly

plink --bfile GW_all_samps --keep $ATL_SAMP_FILE --extract $SNP_FILE \
--r2 square --out $ATL_OUT

mv $ATL_OUT'.ld' $OUT_DIR

# MW8X

cd $START_DIR

MW8X_SAMP_FILE=$SAMP_DIR'MW8X_introgressed_names.plink.txt'
MW8X_OUT=MW8X_Chr08N_ATLvsMW_Fst0.5_RDAonly

plink --bfile GW_all_samps --keep $MW8X_SAMP_FILE --extract $SNP_FILE \
--r2 square --out $MW8X_OUT

mv $MW8X_OUT'.ld' $OUT_DIR

# MW4X No introgression

cd $START_DIR

MW4X_NO_SAMP_FILE=$SAMP_DIR'MW_no_introgress_names.plink.txt'
MW4X_NO_OUT=MW4X_NoIntrogression_Chr08N_ATLvsMW_Fst0.5_RDAonly

plink --bfile GW_all_samps --keep $MW4X_NO_SAMP_FILE --extract $SNP_FILE \
--r2 square --out $MW4X_NO_OUT

mv $MW4X_NO_OUT'.ld' $OUT_DIR
```

## Calculate r^2 for control SNPs
```
# on Cori
module load python/3.7-anaconda-2019.10
source activate plink_1_env

START_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

cd $START_DIR

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5/

FILE_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/

SNP_FILE=$FILE_DIR'Chr08N_ATLvMW_Fst_0.5_MWall_RDA_control_snps.1.plink.txt'

MAF_CUT=0.05

# all Midwest samples, control SNPs

cd $START_DIR

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/

MW_ALL_SAMP_FILE=$SAMP_DIR'MW_all_names.plink.txt'

MW_ALL_OUT=MWall_Chr08N_ATLvsMW_ControlSNPs

plink --bfile GW_all_samps --keep $MW_ALL_SAMP_FILE --extract $SNP_FILE \
--maf $MAF_CUT --r2 square --out $MW_ALL_OUT

mv $MW_ALL_OUT'.ld' $OUT_DIR
```

## Transfer results to HAIB
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/plink_LD_files/Chr08N_ATLvsMW_Fst0.5

scp *ld grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_plink_files/Chr08N_ATLvsMW_Fst0.5_results
```

## Evaluate LD Decay
* see `~/sg_8X/introgression/first_draft_analysis/Chr08N_firstdraft_LD_analsysis.r`



