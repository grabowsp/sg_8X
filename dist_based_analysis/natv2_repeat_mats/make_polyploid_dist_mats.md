# Steps for making distance matrixes for 100 (for now) repeats of 100k subsampled SNPs

## Important info
* Location of results
 `/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist`
* Location of input VCFs
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs`
* R script used for generating distance matrices
  * `/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r`

## Test script
```
module load python/3.7-anaconda-2019.10
module swap PrgEnv-intel PrgEnv-gnu
source activate R_tidy

RSC=/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r

IN_VCF=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/GW.natv2.ancestral.tet.100k.0001.vcf.gz

COMP_TYPE=polyploid

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist

Rscript $RSC $IN_VCF $COMP_TYPE $OUT_DIR

```
## Run for 2-10
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist

sbatch gen_poly_dist_02_10.sh
sbatch gen_poly_dist_11_20.sh

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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist

cd $OUT_DIR

RSC=/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dist_mat_tet_in.r

COMP_TYPE=polyploid

VCF_PRE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/GW.natv2.ancestral.tet.100k.

for SETNUM in {0002..0010};
  do
  echo $SETNUM;
  IN_VCF=$VCF_PRE$SETNUM'.vcf.gz';
  Rscript $RSC $IN_VCF $COMP_TYPE $OUT_DIR;
  done;
```

### Generate list of NJ trees to later analyze with `ape` package
```
### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
library(ape, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')
library(phangorn, lib.loc = '/global/homes/g/grabowsp/tools/r_libs')

### INPUT DATA ###
dist_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist/'

dist_files <- system(paste('ls ', dist_dir, '*rds', sep = ''), intern = T)

dist_list <- list()
for(df in seq(length(dist_files))){
  tmp_data <- readRDS(dist_files[df])
  dist_list[[df]] <- tmp_data[['euclidean_dist']]
}

### SET OUTPUT ###
out_file <- paste(dist_dir, 'natv2.100k.NJtree.list.v1.rds', sep = '')

################

tree_list <- lapply(dist_list, function(x) root(nj(x), 
                                                outgroup = 'ANCESTRAL_TET_GENO',
                                                resolve.root = T))

saveRDS(tree_list, file = out_file)

```





