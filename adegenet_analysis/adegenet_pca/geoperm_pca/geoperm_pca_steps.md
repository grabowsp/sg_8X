# Steps for running PCA on 'geoperm' sample set

## Run PCA
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geoperm_pca

sbatch run_geoperm_full_PCA.sh

ln -s /global/cscratch1/sd/grabowsp/sg_8X_scratch/geoperm_tet_vcfs/GW.50kSNPs.tetrasomic.CDS.geoperm.genlight.PCAresults.rds ./GW.50kSNPs.tetrasomic.CDS.geoperm.genlight.PCAresults.rds
```
