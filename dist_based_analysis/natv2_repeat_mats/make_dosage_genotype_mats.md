# Steps for making dosage-based genotype matrices 
* Use genotypes to look for correlations with odd PCoA results that change depending on disomic or tetrasomic genotypes
  * Variation in PCo4 in Midwest seen with disomic genotypes disappears with tetrasomic genotypes
  * Variation in PCo5 in lowland seen with disomic genotypes disappears with tetrasomic genotypes

## Generate Genotype matrices
```
module load python/3.7-anaconda-2019.10
module swap PrgEnv-intel PrgEnv-gnu
source activate /global/homes/g/grabowsp/.conda/envs/R_tidy

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/dist_analysis/natv2_ancestral_sub_dist

cd $OUT_DIR

RSC=/global/homes/g/grabowsp/tools/sg_8X/dist_based_analysis/generate_dosage_mat_tet_in.r

COMP_TYPE=diploid

VCF_PRE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/GW.natv2.ancestral.tet.100k.

for SETNUM in {0001..0010};
  do
  echo $SETNUM;
  IN_VCF=$VCF_PRE$SETNUM'.vcf.gz';
  Rscript $RSC $IN_VCF $COMP_TYPE $OUT_DIR;
  done;


```
