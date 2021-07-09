# Generate VCFs containing the SNPs for first_draft analysis

## Overview
* Generate VCF for hiFst SNPs chosen for each comparison and each Fst cutoff
 * Will end up being 6 VCFs total
* Uses SNPs that have been LD-pruned and filtered by missing data in training sets
* Workflow is to generate Chromosome VCFs and then concatenate them

## Location of final files
* HAIB Directory:
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_vcfs/`
* NERSC Directory:
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files/`
* Midwest vs Atlantic
 * `AtlanticvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`
  * 90,151 SNPs with Fst > 0.5
 * `AtlanticvsMidwest.hiFst.firstdraft.Fst1.0.disomic.CDS.vcf.gz`
  * 1,714 SNPs with Fst = 1
* Midwest vs Gulf
 * `GulfvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`
  * 82,345 SNPs with Fst > 0.5
 * `GulfvsMidwest.hiFst.firstdraft.Fst1.0.disomic.CDS.vcf.gz`
  * 740 SNPs with Fst = 1
* Atlantic vs Gulf
 * `AtlanticvsGulf.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`
  * 38,762 SNPs with Fst > 0.5
 * `AtlanticvsGulf.hiFst.firstdraft.Fst1.0.disomic.CDS.vcf.gz`
  * 104 SNPs with Fst = 1

## Generate Chromosome VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files 

sbatch gen_firstdraft_VCFs_01_03.sh
sbatch gen_firstdraft_VCFs_04_06.sh
sbatch gen_firstdraft_VCFs_07_09.sh
```
### Example Script
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
    for COMP in AtlanticvsGulf AtlanticvsMidwest GulfvsMidwest;
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

## Concatenate chromosome VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files

for SAMP_PRE in AtlanticvsGulf AtlanticvsMidwest GulfvsMidwest;
 do
 for FST in 0.5 1.0;
  do
  OUT_FILE=$SAMP_PRE'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.vcf';
  cat Chr01K.$SAMP_PRE'.hiFst_Fst'$FST'.firstdraft.recode.vcf' > $OUT_FILE;
  awk '!/^#/ {print}' \
    Chr01N.$SAMP_PRE'.hiFst_Fst'$FST'.firstdraft.recode.vcf' >> $OUT_FILE;
  for CHR_NUM in {02..09};
   do
   for CHR_LET in K N;
    do
    VCF_NAME=Chr$CHR_NUM$CHR_LET'.'$SAMP_PRE'.hiFst_Fst'$FST'.firstdraft.recode.vcf';
    awk '!/^#/ {print}' $VCF_NAME >> $OUT_FILE;
    done;
   done;
  done;
 done;
```

## Compress final files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files

for VCF_FILE in *firstdraft.Fst0.5.disomic.CDS.vcf;
  do
  gzip $VCF_FILE;
  done

for VCF_FILE in *firstdraft.Fst1.0.disomic.CDS.vcf;
  do
  gzip $VCF_FILE;
  done
```

## Transfer VCFs to HudsonAlpha
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files

scp *vcf.gz grabowsk@pants.hagsc.org:/home/f1p1/tmp/switchgrass_8X/firstdraft_vcfs
```

