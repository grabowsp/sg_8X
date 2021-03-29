# Steps for making distance matrixes for 100 (for now) repeats of 100k subsampled SNPs

## Important info
* Location of results
 `/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist`
* Location of input VCFs
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs`
* R script used for generating distance matrices
  * `/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r`


## Generate disomic distance matrices
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist 

sbatch gen_dip_dist_01_10.sh
sbatch gen_dip_dist_11_30.sh
sbatch gen_dip_dist_31_50.sh
sbatch gen_dip_dist_51_70.sh
sbatch gen_dip_dist_71_90.sh
sbatch gen_dip_dist_91_100.sh
# need to wait until remaining vcfs are done to do next rounds

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
module swap PrgEnv-intel PrgEnv-gnu
source activate /global/homes/g/grabowsp/.conda/envs/R_tidy

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist

cd $OUT_DIR

RSC=/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r

COMP_TYPE=diploid

VCF_PRE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/GW.natv2filt.tet.100k.

for SETNUM in {0001..0010};
  do
  echo $SETNUM;
  IN_VCF=$VCF_PRE$SETNUM'.vcf.gz';
  Rscript $RSC $IN_VCF $COMP_TYPE $OUT_DIR;
  done;
```

## Generate tetrasomic distance matrices
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist

sbatch gen_tet_dist_01_10.sh
sbatch gen_tet_dist_11_30.sh
sbatch gen_tet_dist_31_50.sh
sbatch gen_tet_dist_51_70.sh
sbatch gen_tet_dist_71_90.sh
sbatch gen_tet_dist_91_100.sh
# need to wait until remaining vcfs are done to do next rounds
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
module swap PrgEnv-intel PrgEnv-gnu
source activate /global/homes/g/grabowsp/.conda/envs/R_tidy

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist

cd $OUT_DIR

RSC=/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r

COMP_TYPE=tetraploid

VCF_PRE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/GW.natv2filt.tet.100k.

for SETNUM in {0001..0010};
  do
  echo $SETNUM;
  IN_VCF=$VCF_PRE$SETNUM'.vcf.gz';
  Rscript $RSC $IN_VCF $COMP_TYPE $OUT_DIR;
  done;

```

## Generate polyploid distance matrices
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist

sbatch gen_poly_dist_01_10.sh
sbatch gen_poly_dist_11_30.sh
sbatch gen_poly_dist_31_50.sh
sbatch gen_poly_dist_51_70.sh
sbatch gen_poly_dist_71_90.sh
sbatch gen_poly_dist_91_100.sh

# need to wait until remaining vcfs are done to do next rounds
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
module swap PrgEnv-intel PrgEnv-gnu
source activate /global/homes/g/grabowsp/.conda/envs/R_tidy

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2filt_sub_dist

cd $OUT_DIR

RSC=/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r

COMP_TYPE=polyploid

VCF_PRE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/replicated_sub_vcfs/GW.natv2filt.tet.100k.

for SETNUM in {0001..0010};
  do
  echo $SETNUM;
  IN_VCF=$VCF_PRE$SETNUM'.vcf.gz';
  Rscript $RSC $IN_VCF $COMP_TYPE $OUT_DIR;
  done;

```



