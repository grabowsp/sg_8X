# Steps for running PCA on `natv2filt` genepool sample sets

## Run PCA
### Run submit script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/natv2filt_gp_pca

sbatch run_natv2filt_MW_tot_PCA.sh
sbatch run_natv2filt_GULF_tot_PCA.sh
sbatch run_natv2filt_ATL_tot_PCA.sh
sbatch run_natv2filt_MISFIT_tot_PCA.sh 
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/natv2filt_gp_pca

DATA_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2filt_tet_vcfs/GW.100kSNPs.tetrasomic.CDS.natv2filt.MW_tot.genlight.rds

cd $OUT_DIR

Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_pca/adegenet_pca.r $DATA_FILE $OUT_DIR

```


