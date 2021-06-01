# Steps for analyzing Chromsome-level introgression patterns

## Introgression levels across chromosomes
* R script
  `~/sg_8X/introgression/MW_chr_amount.r`

## Analysis of Chr08N
* R script
  `~/sg_8X/introgression/MW_Chr08N_analysis.r`

## Calculate Heterozygosity in MW samples in each chromsome
### Location of files on Cori
* Starting VCFs
  * Parent directory
    * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs`
* Sample name file
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/MW_all_names.txt`
### Run script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

sbatch calc_MW_chr_het.sh
```
### Script
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3/MW_all_names.txt

for CHR_NUM in {01..04};
  do
  for CHR_LET in K N;
    do
    CHR_TXT=Chr$CHR_NUM$CHR_LET
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET;
    OUT_TXT=$OUT_PRE'.MWall';
    #
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --chr $CHR_TXT \
      --keep $SAMP_FILE --max-missing 0.8 --maf 0.01 --het;
    #
    done;
  done;

```
### Combine files and transfer to HAIB
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/v3

tar -czvf MWall_het.tar.gz *het

scp MWall_het.tar.gz grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
```
### Untar files at HAIB
```
cd /home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/

tar -xzvf MWall_het.tar.gz
```

## Process Het file data
* Final table
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MWall_Chr_HOM_vals.txt`
```
bash
source activate R_analysis

library(data.table)

mw_het_files <- system('ls /home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/*het', intern = T)

file_name_pre <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'
file_name_suf <- '.MWall.het'

chr_names <- gsub(file_name_pre, '', 
  gsub(file_name_suf, '', mw_het_files, fixed = T), 
  fixed = T)

res_list <- list()
for(i in seq(length(mw_het_files))){
  res_list[[chr_names[i]]] <- fread(mw_het_files[i])
}

hom_list <- lapply(res_list, function(x) x$'O(HOM)'/x$N_SITES)

hom_mat <- matrix(unlist(hom_list), ncol = length(hom_list), byrow = F)
colnames(hom_mat) <- names(hom_list)

hom_tab <- data.table(samp_name = res_list[[1]]$INDV, hom_mat)

hom_file_out <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MWall_Chr_HOM_vals.txt'

fwrite(hom_tab, file = hom_file_out, sep = '\t')


```
