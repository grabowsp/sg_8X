# Steps for running ADMIXTURE on allsamp samples

## Steps
* generate tped files with VCF tools
* generate bed, .bim, and .fam (?) files using plink
* run admixture
* Look at CV error
* Make barplot(s)
* Make result file(s)

## Overview
* Directory with files and results
  * '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix'

## Generate tped file with VCFtools
### What we need
* list of SNPs that we want to keep - make using genlight object
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix/allsamps_200k_snps.txt`
* "Chromosome map" - chromosome name to integer name map
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/plink_chr_name_map.txt`

### Generate list of SNPs
```
# module load python/3.7-anaconda-2019.10
# source activate adegenet_2_env

### LOAD PACKAGES ###
library(adegenet)

### LOAD INPUTS ###
geno_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/GW.200kSNPs.tetrasomic.CDS.allsamps.genlight.rds'

genos <- readRDS(geno_file)

### SET OUTPUTS ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix/'

out_short <- 'allsamps_200k_snps.txt'

out_full <- paste(out_dir, out_short, sep = '')

#########
loc_list <- strsplit(locNames(genos), split = '_')

chr_vec <- unlist(lapply(loc_list, function(x) x[1]))
pos_vec <- unlist(lapply(loc_list, function(x) x[2]))

snp_df <- data.frame(chr = chr_vec, pos = pos_vec, stringsAsFactors = F)

write.table(snp_df, file = out_full, quote = F, sep = '\t', row.names = F,
  col.names = F)

```
### Generate tped file with VCFtools
#### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix

sbatch gen_allsamp_200k_tped_Chr01_05.sh
sbatch gen_allsamp_200k_tped_Chr06_09.sh

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

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix

VCF_FILE_PRE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/Pvirgatum_1070g_

VCF_FILE_SUF=.snp.sort.norepeats.vcf.gz

SUB_PREFIX=_allsamps_200k

POS_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix/allsamps_200k_snps.txt

SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt

CHROM_MAP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/plink_chr_name_map.txt

cd $OUT_DIR

for CHRNAME in 01K 01N 02K 02N 03K 03N 04K 04N 05K 05N;
  do
  VCF_FILE=$VCF_FILE_PRE'Chr'$CHRNAME$VCF_FILE_SUF;
  OUT_PREFIX=Chr$CHRNAME$SUB_PREFIX;
  vcftools --gzvcf $VCF_FILE --out $OUT_PREFIX --positions $POS_FILE \
    --keep $SAMP_FILE --plink-tped --chrom-map $CHROM_MAP_FILE;
  done
```

# CONTINUE FROM HERE

#### Concatenate .tped files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix

cat *tped > GW_allsamps_200k.tped
cp Chr01K_allsamps_200k.tfam GW_allsamps_200k.tfam
```

## Generte .bed
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/allsamps_admix 

plink --tfile GW_allsamps_200k --maf 0.0001 --make-bed --out GW_200k_allsamps
```
## Run ADMIXTURE
* Example submit script
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
source activate /global/homes/g/grabowsp/.conda/envs/admixture_env

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix

cd $OUT_DIR

IN_FILE=GW_50k_geobigSouthCoastal.bed

K_NUM=1

CV_NUM=10

admixture --cv=$CV_NUM $IN_FILE $K_NUM | tee log${K_NUM}.out

```
### Make additional submit files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix
 
for KN in {2..10};
do
sed 's/K_NUM=1/'K_NUM="$KN"'/g' run_geobigSCoastal_admixture_K1.sh > \
run_geobigSCoastal_admixture_K$KN'.sh';
done

```
### Run jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix 

for KN in {1..10};
do
sbatch run_geobigSCoastal_admixture_K$KN'.sh';
done
```

## Analyze CV error
### Generate CV error file
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix

for K in {1..10};
  do
  grep -h CV log$K'.out' >> cv_error_vals.txt;
  done
```
### Generate line plot
* in R
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

library(data.table)
library(tidyverse)

data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/'

cv_val_file <- paste(data_dir, 'cv_error_vals.txt', sep = '')
cv_vals <- fread(cv_val_file, header = F)

CV_vals <- cv_vals$V4
K_num <- as.numeric(gsub('):', '',
  gsub('(K=', '', cv_vals$V3, fixed = T), fixed = T))

CV_dt <- data.table(K_num, CV_vals)

gg_cv <- ggplot(CV_dt, aes(x = K_num, y = CV_vals)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(10)) +
  xlab('K') + ylab('CV error') +
  ggtitle('ADMIXTURE CV error for geobig_NorthInland\nsamps and 50k SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/geobig_SouthCoastal_50k_CV_error.pdf'

pdf(out_file, width = 4.5, height = 4.5)
gg_cv
dev.off()
```

## Generate barplots
### K = 2
* R script for K=2
  * `/home/grabowsky/tools/workflows/sg_8X/admixture_analysis/geobig_SCoast_admix/geobig_SCoast_admix_K2_barplot.r`
* Figure
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.2.ADMIXTURE.memb.pdf`
### K = 4
* R script for K=4
  * `/home/grabowsky/tools/workflows/sg_8X/admixture_analysis/geobig_SCoast_admix/geobig_SCoast_admix_K4_barplot.r`
* Figure
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.4.ADMIXTURE.memb.pdf`
### K = 5
* R script for K=5
  * `/home/grabowsky/tools/workflows/sg_8X/admixture_analysis/geobig_SCoast_admix/geobig_SCoast_admix_K5_barplot.r`
* Figure
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.5.ADMIXTURE.memb.pdf`


## Generate Results File
* File paths
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.2.results.txt`
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.4.results.txt`
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/GW_50k_geobigSouthCoastal.5.results.txt`
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_SouthCoastal_admix/

cut -d " " -f 1 GW_50k_geobigSouthCoastal.fam | \
paste -d " " - GW_50k_geobigSouthCoastal.2.Q > \
GW_50k_geobigSouthCoastal.2.results.txt

cut -d " " -f 1 GW_50k_geobigSouthCoastal.fam | \
paste -d " " - GW_50k_geobigSouthCoastal.4.Q > \
GW_50k_geobigSouthCoastal.4.results.txt

cut -d " " -f 1 GW_50k_geobigSouthCoastal.fam | \
paste -d " " - GW_50k_geobigSouthCoastal.5.Q > \
GW_50k_geobigSouthCoastal.5.results.txt


```


