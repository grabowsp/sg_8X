# Steps for making `all_samp` VCFs for CDS SNPs

## Generate `all_samp` CDS VCFs
* Result directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs`
* Data directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs`
* Sample file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt`
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs

sbatch generate_Chr01K_01N_tet_CDS_allsamps_vcf.sh
sbatch generate_Chr02_Chr05_tet_CDS_allsamps_vcf.sh

# NEED TO RUN THESE ONCE DONE TRANSFERRING DATA
#sbatch generate_Chr06_Chr09_tet_CDS_allsamps_vcf.sh
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

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs/
BED_FILE=/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/
SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt
SAMPSET=allsamps

cd $OUT_DIR

for CHROM in 01K 01N;
  do
  vcftools --gzvcf \
    $DATA_DIR'Pvirgatum_1070g_Chr'$CHROM'.5alleles.snp.sort.norepeats.vcf.gz' \
    --max-missing 0.8
    --stdout --bed $BED_FILE --keep $SAMP_FILE --recode --recode-INFO-all | \
    gzip -c > $OUT_DIR'Chr'$CHROM'.tetrasomic.CDS.'$SAMPSET'.vcf.gz';
  gunzip -kc $OUT_DIR'Chr'$CHROM'.tetrasomic.CDS.'$SAMPSET'.vcf.gz' | \
    split -l 100000 -d - $OUT_DIR'Chr'$CHROM'.tetrasomic.CDS.'$SAMPSET'.vcf_';  
  done

```

## Generate sample header file for VCFs
* File name: 
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/CDS.tetrasomic.allsamps.vcf.sampheader.txt`
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs

SAMPSET=allsamps

gunzip -kc Chr01K.tetrasomic.CDS.$SAMPSET'.vcf.gz' | head -634 | \
tail -1 > CDS.tetrasomic.$SAMPSET'.vcf.sampheader.txt'
```

## Generate `genlight` object for each chromosome
* If want 6+ copies of the MAF, then use 0.0014 as MAF
  * this is pretty permissive because assumes all samps are 8X

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

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/r_adegenet_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/

cd $DATA_DIR

HEADER_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/CDS.tetrasomic.allsamps.vcf.sampheader.txt

SAMPSET=allsamps

MAF_CUT=0.0014

for CHR_N in 01;
  do
  for CHR_T in K N;
    do
    TEST_CHR=$CHR_N$CHR_T
    SEARCH_STRING=Chr$TEST_CHR'.tetrasomic.CDS.'$SAMPSET'.vcf_'
    Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/make_Chr_genlight_objs.r \
    $DATA_DIR $SEARCH_STRING'*' $HEADER_FILE $MAF_CUT;
    done;
  done;



```
