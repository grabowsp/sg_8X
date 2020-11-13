# Steps for generating CDS genotype files for `geoperm` sample set

## Generate `geoperm` CDS VCFs
### Location of files
* Result directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/geoperm_tet_vcfs`
* Data directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs`
* Sample file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_perm_names.txt`
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/geoperm_tet_vcfs

sbatch generate_Chr01_05_tet_CDS_geoperm_vcf.sh
sbatch generate_Chr06_09_tet_CDS_geoperm_vcf.sh
```
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
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs/
BED_FILE=/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/geoperm_tet_vcfs/
SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_perm_names.txt
SAMPSET=geoperm

cd $OUT_DIR

for CHROM in 01K 01N 02K 02N 03K 03N 04K 04N 05K 05N;
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