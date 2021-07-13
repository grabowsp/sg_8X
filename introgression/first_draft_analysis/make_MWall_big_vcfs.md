# Generate VCFs for MWall sample set

## Location of final files
* NERSC:
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_big_vcfs`
* HAIB:
 * `/home/f1p1/tmp/switchgrass_8X/MWall_big_vcfs`

## Generate VCFs
```
module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_big_vcfs/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt

MAF_CUT=0.04

for CHR_NUM in {07..08};
 do
 for CHR_LET in K N;
  do
  VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
  OUT_PRE=Chr$CHR_NUM$CHR_LET;
  OUT_TXT=$OUT_PRE'.MWall.big'
  #
  vcftools --gzvcf $VCF_FULL --out $OUT_TXT --keep $SAMP_FILE --maf $MAF_CUT \
  --recode --recode-INFO --recode-INFO-all;
  done;
 done;
```

## Compress Files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_big_vcfs

for CHR_NUM in {07..09};
 do
 for CHR_LET in K N;
  do
  gzip Chr$CHR_NUM$CHR_LET'.MWall.big.recode.vcf';
  done;
 done;
```

## Transfer to HAIB
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_big_vcfs

ls *vcf.gz

scp *vcf.gz grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/MWall_big_vcfs
```

