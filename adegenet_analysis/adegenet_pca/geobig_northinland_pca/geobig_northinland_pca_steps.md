# Steps for running PCA on 'geobig_northinland' sample set

## Run PCA
### Run submit script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_northinland_pca

sbatch run_geobig_northinland_PCA.sh
```
### Example of script
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

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_northinland_pca/

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_northinland_tet_vcfs/GW.50kSNPs.tetrasomic.CDS.geobig_northinland.genlight.rds

cd $OUT_DIR

Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_pca/adegenet_pca.r $DATA_FILE $OUT_DIR


```


