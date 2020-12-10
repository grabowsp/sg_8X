# Steps for measuring diversity measures in `geo_big` samples

## Overview
* NERSC Resulst directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity`

## Calculate homozygosity in CDS across chromosomes
### Run script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity

sbatch calc_geobig_div.sh
```

```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
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
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity
DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/
DATA_FILE_SUF=.tetrasomic.CDS.geobig.genlight.rds

OUT_PRE_SUB=geobig.CDS

BLOCKSIZE=1e5

cd $OUT_DIR

for CHR_NUM in {01..09};
  do
  for CHR_LET in K N;
    do
    DATA_FILE=$DATA_DIR'Chr'$CHR_NUM$CHR_LET$DATA_FILE_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET'.'$OUT_PRE_SUB;
    Rscript /global/homes/g/grabowsp/tools/sg_8X/diversity_measures/calc_hom_het_vals.r \
    $DATA_FILE $OUT_DIR $OUT_PRE $BLOCKSIZE;
    done;
  done;




```
