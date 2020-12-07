# Steps for running ADMIXTURE on geobig samples

## Steps
* generate tped files with VCF tools
* generate bed, .bim, and .fam (?) files using plink
* run admixture
* Look at CV error
* Make barplot(s)
* Make result file(s)

## Overview
* Directory with files and results
  * '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix'

## Generate tped file with VCFtools
### What we need
* list of SNPs that we want to keep - make using genlight object
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/geobig_50k_snps.txt`
* "Chromosome map" - chromosome name to integer name map
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/plink_chr_name_map.txt`

### Generate list of SNPs
```
# module load python/3.7-anaconda-2019.10
# source activate adegenet_2_env

### LOAD PACKAGES ###
library(adegenet)

### LOAD INPUTS ###
geno_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.rds'

genos <- readRDS(geno_file)

### SET OUTPUTS ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/'

out_short <- 'geobig_50k_snps.txt'

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
#### Test for 1 chromosome
```
# module load python/3.7-anaconda-2019.07
# source activate gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix

VCF_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/Pvirgatum_1070g_Chr01K.snp.sort.norepeats.vcf.gz

OUT_PREFIX=Chr01K_geobig_50k

POS_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/geobig_50k_snps.txt

SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/geo_big_names.txt

CHROM_MAP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/plink_chr_name_map.txt

cd $OUT_DIR

vcftools --gzvcf $VCF_FILE --out $OUT_PREFIX --positions $POS_FILE \
--keep $SAMP_FILE --plink-tped --chrom-map $CHROM_MAP_FILE

```
#### Submit jobs for remaining chromosomes
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix
sbatch gen_geobig_tped_Chr01N_03.sh
sbatch gen_geobig_tped_Chr04_06.sh
sbatch gen_geobig_tped_Chr07_09.sh
```
#### Concatenate .tped files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix

cat *tped > GW_geobig_50k.tped
cp Chr01K_geobig_50k.tfam GW_geobig_50k.tfam
```

## Generte .bed
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix

plink --tfile GW_geobig_50k --maf 0.0001 --make-bed --out GW_50k_geobig
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

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix

cd $OUT_DIR

IN_FILE=GW_50k_geobig.bed

K_NUM=1

CV_NUM=10

admixture --cv=$CV_NUM $IN_FILE $K_NUM | tee log${K_NUM}.out

```
### Make additional submit files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix

for KN in {2..10};
do
sed 's/K_NUM=1/'K_NUM="$KN"'/g' run_geobig_admixture_K01.sh > \
run_geobig_admixture_K$KN'.sh';
done

```
### Run jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix

sbatch run_geobig_admixture_K01.sh

for KN in {2..10};
do
sbatch run_geobig_admixture_K$KN'.sh';
done
```

## Analyze CV error
* in R
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

library(data.table)
library(tidyverse)

# manually inputted this from the cori output files; should save log files
#  in the future
CV_vals <- c(0.13281, 0.10978, 0.09949, 0.09734, 0.09603, 0.09483, 0.09410, 
  0.09377, 0.09328, 0.09336)

K_num <- c(1:10)

CV_dt <- data.table(K_num, CV_vals)

gg_cv <- ggplot(CV_dt, aes(x = K_num, y = CV_vals)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(10)) +
  xlab('K') + ylab('CV error') +
  ggtitle('ADMIXTURE CV error for geobig samps\nand 50k SNPs')

out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/geobig_50k_CV_error.pdf'

pdf(out_file, width = 4.5, height = 4.5)
gg_cv
dev.off()
```

## Generate barplots
* R script
  * `/home/grabowsky/tools/workflows/sg_8X/admixture_analysis/geobig_admix/admix_barplot_template.r`
* Figure
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.ADMIXTURE.memb.pdf`

## Generate Results File
* File path
  `/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/GW_50k_geobig.3.results.txt`
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/

cut -d " " -f 1 GW_50k_geobig.fam | paste -d " " - GW_50k_geobig.3.Q > \
GW_50k_geobig.3.results.txt

cut -d " " -f 1 GW_50k_geobig.fam | paste -d " " - GW_50k_geobig.2.Q > \
GW_50k_geobig.2.results.txt

cut -d " " -f 1 GW_50k_geobig.fam | paste -d " " - GW_50k_geobig.4.Q > \
GW_50k_geobig.4.results.txt
```




