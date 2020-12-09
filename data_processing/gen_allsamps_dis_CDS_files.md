# Steps for making `all_samp` VCFs for CDS SNPs

## Generate `all_samp` CDS VCFs
* Result directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs`
* Data directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files`
* Sample file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt`

## Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs

sbatch generate_Chr01_Chr05_CDS_allsamps_vcf.sh
sbatch generate_Chr06_Chr09_CDS_allsamps_vcf.sh

```
### Example script used for submitting jobs
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/
BED_FILE=/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt
SAMPSET=allsamps

cd $OUT_DIR

for CHROM in 01K 01N 02K 02N 03K 03N 04K 04N 05K 05N;
  do
  vcftools --gzvcf \
    $DATA_DIR'Pvirgatum_1070g_Chr'$CHROM'.snp.sort.norepeats.vcf.gz' \
    --max-missing 0.8 \
    --stdout --bed $BED_FILE --keep $SAMP_FILE --recode --recode-INFO-all | \
    gzip -c > $OUT_DIR'Chr'$CHROM'.disomic.CDS.'$SAMPSET'.vcf.gz';
  gunzip -kc $OUT_DIR'Chr'$CHROM'.disomic.CDS.'$SAMPSET'.vcf.gz' | \
    split -l 100000 -d - $OUT_DIR'Chr'$CHROM'.disomic.CDS.'$SAMPSET'.vcf_';  
  done

```

## Generate sample header file for VCFs
* File name:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/CDS.disomic.allsamps.vcf.sampheader.txt`
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs
 
SAMPSET=allsamps

gunzip -kc Chr01K.disomic.CDS.$SAMPSET'.vcf.gz' | head -648 | \
tail -1 > CDS.disomic.$SAMPSET'.vcf.sampheader.txt'
```

## Generate `genlight` object for each chromosome
* If want 6+ copies of the MAF, then use 0.0015 as MAF
  * this is pretty permissive because assumes all samps are 8X
### Submit scripts
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs

sbatch gen_allsamp_dip_Chr01_05_genlight.sh
sbatch gen_allsamp_dip_Chr06_09_genlight.sh

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
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/

cd $DATA_DIR

HEADER_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/CDS.disomic.allsamps.vcf.sampheader.txt

SAMPSET=allsamps

MAF_CUT=0.0015

for CHR_N in 01 02 03 04 05;
  do
  for CHR_T in K N;
    do
    TEST_CHR=$CHR_N$CHR_T
    SEARCH_STRING=Chr$TEST_CHR'.disomic.CDS.'$SAMPSET'.vcf_'
    Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/make_dip_only_Chr_genlight_objs.r \
    $DATA_DIR $SEARCH_STRING'*' $HEADER_FILE $MAF_CUT;
    done;
  done;

```

## Generate genome-wide subsampled `genlight` object 
* Subsample SNPs from all chromosomes to get desired number of genome-wide SNPs
* Final file:
  `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/GW.50kSNPs.disomic.CDS.allsamps.genlight.rds`
### Get number of SNPs in each chromosome `genlight` object
* SNP count file
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/allsamps_diploid.SNPcount.txt`
```
module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
FILE_SUB_SHORT=allsamps.genlight.rds
OUT_SHORT=allsamps_diploid

cd $DATA_DIR

Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/get_tot_nSNPs.r \
$DATA_DIR '*'$FILE_SUB_SHORT $OUT_SHORT
```
### Calculate sub-sampling rate
* in R
```
# module load python/3.7-anaconda-2019.10
# source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

library(data.table)
res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/allsamps_diploid.SNPcount.txt'
res <- fread(res_file)

goal_n <- 1e5

goal_n / sum(res$nSNPs)
# [1] 0.005905933
```
### Generate subsampled `genlight` object
* Is faster to just run in interactive session
```
module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs
FILE_SUB_SHORT=allsamps.genlight.rds
OUT_SHORT='GW.50kSNPs.disomic.CDS.allsamps.genlight.rds'
PER_SUBSAMP=0.0065
TOT_SNP=1e5

cd $DATA_DIR

Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/subsample_genlight.r \
$DATA_DIR '*'$FILE_SUB_SHORT $OUT_SHORT $PER_SUBSAMP $TOT_SNP
```






