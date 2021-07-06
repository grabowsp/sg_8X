# Generate VCFs containing the SNPs for first_draft analysis

## Steps
1. Final all files
2. Generate VCF with all samples

## Location of files
* Original VCF location on NERSC:
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs`
* Output directory
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft`

## Generate Chromosome VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files 

sbatch gen_firstdraft_VCFs_01_03.sh

sbatch gen_firstdraft_VCFs_04_06.sh
sbatch gen_firstdraft_VCFs_07_09.sh
```
## Example Script
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files

SNP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

for CHR_NUM in {01..03};
  do
  for CHR_LET in K N;
    do
    for COMP in AtlanticVsGulf AtlanticVsMidwest GulfVsMidwest;
      do
      for FST in 0.5 1.0;
        do
        #
        VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
        OUT_PRE=Chr$CHR_NUM$CHR_LET;
        #
        OUT_TXT=$OUT_PRE'.'$COMP'.hiFst_Fst'$FST'.firstdraft';
        #
        SNP_POS_FILE=$SNP_DIR$COMP'_hiFst_SNPs_firstdraft_Fst'$FST'.LD.MISS.snps.txt'
        #
        vcftools --gzvcf $VCF_FULL --out $OUT_TXT --positions $SNP_POS_FILE \
        --recode --recode-INFO --recode-INFO-all;
        done;
      done;
    done;
  done;
```
