# Steps for generating pseudohaploid genotype files for 'geobig'

## Generate 'genlight' object for each chromosome
* If want 6+ copies of the MAF, then use 0.0075 as MAF
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/

sbatch gen_geobig_pseudohap_genlight_Chr01.sh

sbatch gen_geobig_pseudohap_genlight_Chr02_03.sh
sbatch gen_geobig_pseudohap_genlight_Chr04_05.sh
sbatch gen_geobig_pseudohap_genlight_Chr06_07.sh
sbatch gen_geobig_pseudohap_genlight_Chr08_09.sh


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
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/

cd $DATA_DIR

HEADER_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/CDS.tetrasomic.geobig.vcf.sampheader.txt

SAMPSET=geobig

MAF_CUT=0.0075

for CHR_N in 01;
  do
  for CHR_T in K N;
    do
    TEST_CHR=$CHR_N$CHR_T
    SEARCH_STRING=Chr$TEST_CHR'.tetrasomic.CDS.'$SAMPSET'.vcf_'
    Rscript /global/homes/g/grabowsp/tools/sg_8X/adegenet_analysis/adegenet_genotype_generation/make_pseudohap_Chr_genlight_objs.r \
    $DATA_DIR $SEARCH_STRING'*' $HEADER_FILE $MAF_CUT;
    done;
  done;

```


