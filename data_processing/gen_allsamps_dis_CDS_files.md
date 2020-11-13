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

sbatch generate_Chr01K_01N_CDS_allsamps_vcf.sh
sbatch generate_Chr02_Chr05_CDS_allsamps_vcf.sh
sbatch generate_Chr06_Chr09_CDS_allsamps_vcf.sh

```

## Example script used for submitting jobs
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

for CHROM in 01K 01N;
  do
  vcftools --gzvcf \
    $DATA_DIR'Pvirgatum_1070g_Chr'$CHROM'.snp.sort.norepeats.vcf.gz' \
    --stdout --bed $BED_FILE --keep $SAMP_FILE --recode --recode-INFO-all | \
    gzip -c > $OUT_DIR'Chr'$CHROM'.disomic.CDS.'SAMPSET'.vcf.gz';
  gunzip -kc $OUT_DIR'Chr'$CHROM'.disomic.CDS.'SAMPSET'.vcf.gz' | \
    split -l 100000 -d - $OUT_DIR'Chr'$CHROM'.disomic.CDS.'SAMPSET'.vcf_';  
  done

```

