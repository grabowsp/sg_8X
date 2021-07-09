# Selecting SNPs that can be used as a control to compare to RDA results

## Overview
* RDA results show a lot of significant associtaions
* Want a set of random SNPs in MWall that has similar allele frequency spectrum

## Plan
* Use VCFtools to get MAF from each chromosome
  * use CDS snps
  * make sure to select MWall samples
* Use VCFtools to extract the MAF for just the RDA SNPs
* Divide RDA SNPs by Chr, then sort by MAF, then split by deciles (1/10)
* Select number of control SNPs to select by Chr, then random sample within decile allele frequency ranges

## Location of Files
### Allele Frequency Files
* Look for: `Chr*MWall.frq`
* on NERSC
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft`
* on HAIB
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files`

## Extract allele frequency info from VCFs
### Generate Sample Name File
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

library(data.table)

res_file <- '/global/homes/g/grabowsp/data/switchgrass/results_tables_8X/natv2filt_res_tab_v4.1.txt'
res_tab <- fread(res_file)

mw_names <- res_tab[grep('MW_', res_tab$subgrp_v2), list(samp_name)]

mw_name_outfile <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt'

fwrite(mw_names, mw_name_outfile, col.names = F, row.names = F)
```
### Get allele frequencies for MWall
* I ended up running interactively to make sure finished before Cori shutdown
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

sbatch get_MWall_CDS_freqs_01_03.sh
```
#### Example Script
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

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt

for CHR_NUM in {07..09};
  do
  for CHR_LET in K N;
    do
    #
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET;
    #
    OUT_TXT=$OUT_PRE'.MWall';
    #
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --keep $SAMP_FILE --freq;
    done;
  done;

```

## Get MAF for the hiFST SNPs
* Chose to run interactively to make sure finished before Cori shutdown
```
module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MWall_sampnames.txt

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files/

for COMPS in AtlanticvsMidwest GulfvsMidwest;
 do
 for FST in 0.5 1.0;
  do
  VCF_NAME=$VCF_DIR$COMPS'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.vcf.gz';
  OUT_TXT=$COMPS'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.MWall';
  #
  vcftools --gzvcf $VCF_NAME --out $OUT_TXT --keep $SAMP_FILE --freq;
  done;
 done;

```

## Package and transfer to HAIB
```
# cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files

# mv Chr*frq ..

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

ls *MWall.frq

tar -czvf MWall_CDS_freqs.tar.gz *MWall.frq

scp MWall_CDS_freqs.tar.gz grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files

# on HAIB

cd /home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files

tar -xzvf MWall_CDS_freqs.tar.gz
```





