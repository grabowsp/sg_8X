# Steps for generating allsamps VCF that includes ancestral genotypes

## Steps
1. Make new directory
1. Make diploid VCF (from original VCF) using 'allsamps' and SNP positions of the ancestral.VCF files
2. Copy ancestral VCFs to working directory
2. Sort, bgzip, and tabix the sample and ancestral VCF files prior to merging
  * see: `~/sg_ha_effort/polyploid_genos/genotype_file_generation/combine_sim8X_geo.md`
3. Merge VCFs

## Location of Files
* Result directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs`
* Original VCF directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files`
* allsamps samp name file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt`
* ancestral VCF snp positions
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/$CHR_NAME'.ancestral_AP13_kept_SNP_positions.txt'`
* ancestral diploid VCF
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/$CHR_NAME'Chr09N.ancestral_AP13_dip_format.vcf'`

## Generate, sort, bgzip, and tabix sample VCFs
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs

sbatch gen_tmp_allsamp_ancest_pos_vcfs_Chr01.sh
sbatch gen_tmp_allsamp_ancest_pos_vcfs_Chr02_05.sh
sbatch gen_tmp_allsamp_ancest_pos_vcfs_Chr06_09.sh

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

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/
ANCESTRAL_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/
ANCESTRAL_SNP_SUF=.ancestral_AP13_kept_SNP_positions.txt
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs/

SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt
SAMPSET=allsamps

SORT_SUF='.'$SAMPSET'.ancest_snp_pos.dip.sort.vcf.gz'

cd $OUT_DIR

for CHROM in 01K 01N;
  do
  SORT_OUT=Chr$CHROM$SORT_SUF;
  vcftools --gzvcf \
    $DATA_DIR'Pvirgatum_1070g_'$CHROM'.snp.sort.norepeats.vcf.gz' \
    --stdout --positions $ANCESTRAL_DIR'Chr'$CHROM$ANCESTRAL_SNP_SUF \
    --keep $SAMP_FILE --recode --recode-INFO-all | \
  bcftools sort - -Oz -o $SORT_OUT;
  tabix -p vcf $SORT_OUT;
  done;
```

## Sort, bgzip, and tabix ancestral VCFs
* saved into working directory
### Submit script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs

sbatch gen_ancest_sort_vcf_Chr01.sh
sbatch gen_ancest_sort_vcf_Chr02_09.sh

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

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs/
ANCEST_FILE_SUF=.ancestral_AP13_dip_format.vcf

cd $OUT_DIR

for CHROM in 01K 01N;
  do
  ANCEST_FILE=$DATA_DIR'Chr'$CHROM$ANCEST_FILE_SUF;
  SORT_OUT=Chr$CHROM'.ancestral.dip.sort.vcf.gz';
  bcftools sort $ANCEST_FILE -Oz -o $SORT_OUT;
  tabix -p vcf $SORT_OUT;
  done;

```

## Merge sample and ancestral VCFs
### Run scripts
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs

sbatch merge_allsamps_ancest_vcfs_Chr01.sh
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
#SBATCH -t 24:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs/
SAMP_VCF_SUF=.allsamps.ancest_snp_pos.dip.sort.vcf.gz
ANC_VCF_SUF=.ancestral.dip.sort.vcf.gz
OUT_SUF=.allsamps.ancestral.sort.vcf.gz

cd $DATA_DIR

for CHROM in 01K 01N;
  do
  SAMP_VCF=Chr$CHROM$SAMP_VCF_SUF;
  ANC_VCF=Chr$CHROM$ANC_VCF_SUF;
  OUT_FILE=Chr$CHROM$OUT_SUF;
  bcftools merge $SAMP_VCF $ANC_VCF -Oz -o $OUT_FILE;
  done;




```
