# Scratch space for figuring out how to best use tools for filtering samples

## Overview
* I think I will use plink 2.0
* I will start by using disomic genotypes
* Try calculating a few distance metrics
  * plink 1.9
    * `--distance` is the command for Identity by state calculation
      * may want to use `square0` flag to make square matrix for downstream analysis and so that can use `--parallel` flag
      * may want to use `1-ibs` to get distance equivalent
    * `--make-rel` makes 'realized relationship matrix' 
      * use `square0` flage to be able to use `--parallel` 
  * plink 2
    * `--sample-diff` looks like Identify-by-state calculation
      * Looks like this works with dosage, too, with `dosage` flag

## Workflow
* start with Chr01K
* set up plink 1.9 environment on Cori
* Need to generate chromosome name map - look at my files
* Generate tped files using VCF tools so that chromosome names get converted
to format that works for plink
* Input VCF into plink using all_samp names, CDS positions, 
  * double check that name file is properly formatted
  * double check that CDS position file is proper, I think need `--extract range` flag and modifier
* Prune SNPs by LD
* Save genotypes as binary plink file
* Calculate distance matrix

## Location of files on Cori
* Directory with all_samp tped files
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds`
* Chromosome map file
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/sg_chr_map.txt`
* CDS BED file
  * `/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed`
* all_samps samp name file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt`
* Location of VCFs
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files`

## Generate tped files
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds

sbatch gen_Chr01_tpeds.sh
sbatch gen_Chr02_05_tpeds.sh
sbatch gen_Chr06_09_tpeds.sh
```
### Example Script for generating tped file using VCFtools
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
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/
VCF_FILE_PRE=Pvirgatum_1070g_Chr
VCF_FILE_SUF=.snp.sort.norepeats.vcf.gz

CDS_BED=/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed

SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/all_samp_names.txt

CHROM_MAP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/sg_chr_map.txt

MAF_LOW_CUT=0.003

MAF_HI_CUT=0.997

MAX_MISS=0.8

cd $OUT_DIR

for CHR_NUM in 01;
  do
  for CHR_LET in K N;
  do
  OUT_PREFIX=all_samps_Chr$CHR_NUM$CHR_LET;
  VCF_FILE=$VCF_DIR$VCF_FILE_PRE$CHR_NUM$CHR_LET$VCF_FILE_SUF;
  #
  vcftools --gzvcf $VCF_FILE --out $OUT_PREFIX --bed $CDS_BED \
    --keep $SAMP_FILE --maf $MAF_LOW_CUT --max-maf $MAF_HI_CUT \
    --max-missing $MAX_MISS \
    --plink-tped --chrom-map $CHROM_MAP_FILE;
  done
done
```
# CONTINUE FROM HERE
### Concatenat .tped files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds

cat *tped > GW_all_samps.tped
cp all_samps_Chr01K.tfam GW_all_samps.tfam
```

## Generate genome-wide plink .bed file
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds

plink -tfile GW_all_samps --maf 0.0001 --make-bed --out GW_all_samps
```

## Select LD-pruned SNPs
* SNP list
  * should be: `GW_all_samps_ld0.3.prune.in`
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds

plink --bfile GW_all_samps --indep-pairwise 50'kb' 100 0.3 --out GW_all_samps_ld0.3
```

## Generate distance matrix using LD-pruned SNPs
```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds

plink --bfile GW_all_samps --distance square '1-ibs' --extract GW_all_samps_ld0.3.prune.in --out GW_all_samps_ld0.3_symmetric
```





* Generate tped files
* Generate plink file using all chromosomes
* Prune by LD
* Calculate distance matrix




## Test on HA while NERSC is down
## Overview
* Test process for Chr01K on HA, then do entire process on NERSC
  * I had trouble consistently accessing the original VCFs from Sujan's scratch space, so decided to just test using one chromosome to save time and preserve space
  * I needed to copy the vcf to my directory, so doing that for all 18 chromosomes would take up too much space
* Sujan's scratch directory
  * `/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/`
* My directory where I transfered the original VCFs
  * `/home/f2p1/work/grabowsk/data/switchgrass/sujan_vcf_copies/`
## Chromosome map file for VCFtools
* `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/sg_chr_map.txt`
```
out_file <- '/home/f2p1/work/grabowsk/data/switchgrass/plink_files/sg_chr_map.txt'

chr_names <- paste('Chr0', rep(seq(9), each = 2), c('K', 'N'), sep = '')
chr_nums <- seq(length(chr_names))

chr_df <- data.frame(name = chr_names, num = chr_nums, stringsAsFactors = F)

write.table(chr_df, out_file, quote = F, sep = '\t', row.names = F, 
  col.names = F)
```
### Script for generating tped file using VCFtools
* the meat of the script can be adapted for use with Cori submissions
  * need to add more chromosome numbers and the N subgenome
```
#!/bin/bash

#$ -q "all.q"
#$ -pe smp 1
#$ -M pgrabowski@hudsonalpha.org
#$ -m ae
#$ -V
#$ -cwd
#$ -N make_tpeds_02

source /home/raid2/LINUXOPT/miniconda2/bin/activate \
/home/grabowsky/.conda/envs/bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps

VCF_DIR=/home/f2p1/work/grabowsk/data/switchgrass/sujan_vcf_copies/
VCF_FILE_PRE=Pvirgatum_1070g_Chr
VCF_FILE_SUF=.snp.poly.sort.norepeats.vcf.gz

CDS_BED=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/sg_v5_CDS.bed

SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samp_names.txt

CHROM_MAP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/sg_chr_map.txt

MAF_LOW_CUT=0.003

MAF_HI_CUT=0.997

MAX_MISS=0.8

cd $OUT_DIR

for CHR_NUM in 01;
  do
  for CHR_LET in K;
  do
  OUT_PREFIX=all_samps_Chr$CHR_NUM$CHR_LET;
  VCF_FILE=$VCF_DIR$VCF_FILE_PRE$CHR_NUM$CHR_LET$VCF_FILE_SUF;
  #
  vcftools --gzvcf $VCF_FILE --out $OUT_PREFIX --bed $CDS_BED \
    --keep $SAMP_FILE --maf $MAF_LOW_CUT --max-maf $MAF_HI_CUT \
    --max-missing $MAX_MISS \
    --plink-tped --chrom-map $CHROM_MAP_FILE;
  done
done

```
### Submit scripts
```
cd /home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps

qsub gen_tped_01K.sh
qsub gen_tped_01N.sh
```
### Generate plink .bed
* for all chromosomes, can concatenate the tped, and copy/rename a .fam file, then use this same type of command
```
bash
source activate plink_1.9_env

cd /home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps

plink -tfile all_samps_Chr01K --maf 0.0001 --make-bed --out all_samps_Chr01K
```
### Generate distance matrix
```
bash
source activate plink_1.9_env

cd /home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps

plink --bfile all_samps_Chr01K --distance square0 '1-ibs' --out all_samps_Chr01K

plink --bfile all_samps_Chr01K --distance square '1-ibs' --out all_samps_Chr01K_symmetric

plink --bfile all_samps_Chr01K --make-rel square0 --out all_samps_Chr01K


```
#### Try pruning by LD
```
bash
source activate plink_1.9_env

cd /home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps

plink --bfile all_samps_Chr01K --indep-pairwise 50'kb' 100 0.3 --out all_samps_Chr01K_ld0.3

plink --bfile all_samps_Chr01K --distance square '1-ibs' --extract all_samps_Chr01K_ld0.3.prune.in --out all_samps_Chr01K_ld0.3_symmetric 
```

### Analysis in R
#### Important files
* Metadata file
  * `/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_metadata_v1.0.csv`
* square0 distance files
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K.mdist`
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K.mdist.id`
* symmetric distance files
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_symmetric.mdist`
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_symmetric.mdist.id`
* LD-pruned symmetric distance files
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_ld0.3_symmetric.mdist`
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_ld0.3_symmetric.mdist.id`

```
bash
source activate R_analysis

library(data.table)

meta_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

weird_locality_inds <- c(53,797,921,924,936)

dist_res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K.mdist'
dist_res <- fread(dist_res_file, header = F)
dist_ids <- fread(paste(dist_res_file, 'id', sep = '.'), header = F)
dist_mat <- as.matrix(dist_res)
colnames(dist_mat) <- rownames(dist_mat) <- dist_ids$V1

sym_dist_res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_symmetric.mdist'
sym_dist_res <- fread(sym_dist_res_file, header = F)
sym_dist_res_ids <- fread(paste(sym_dist_res_file, 'id', sep = '.'), header = F)
sym_dist_mat <- as.matrix(sym_dist_res)
colnames(sym_dist_mat) <- rownames(sym_dist_mat) <- sym_dist_res_ids$V1

ld_sym_dist_res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_ld0.3_symmetric.mdist'
ld_sym_dist_res <- fread(ld_sym_dist_res_file, header = F)
ld_sym_dist_ids <- fread(paste(ld_sym_dist_res_file, 'id', sep = '.'), 
  header = F)
ld_sym_dist_mat <- as.matrix(ld_sym_dist_res)
colnames(ld_sym_dist_mat) <- rownames(ld_sym_dist_mat) <- ld_sym_dist_ids$V1


cor(c(sym_dist_mat), c(ld_sym_dist_mat))
# [1] 0.9426157

samp_cor_vec <- c()
for(i in seq(ncol(sym_dist_mat))){
  tmp_val <- cor(sym_dist_mat[i,], ld_sym_dist_mat[i,])
  samp_cor_vec <- c(samp_cor_vec, tmp_val)
}

summary(samp_cor_vec)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8330  0.9416  0.9519  0.9557  0.9866  0.9908

sum(samp_cor_vec < 0.9)
# 47

summary(setdiff(c(ld_sym_dist_mat), 0))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.003156 0.046690 0.053973 0.053334 0.060379 0.080834

low_dist_inds <- which(ld_sym_dist_mat > 0 & ld_sym_dist_mat < 0.02, 
  arr.ind = T)

# Continue from here

218, 481, 522, 924

tmp_test_inds <- which(ld_sym_dist_ids$V1 %in% 
  samp_meta$VCF_NAME[c(218, 481, 522, 924)])

ld_sym_dist_mat[tmp_test_inds, tmp_test_inds]

ld_sym_dist_mat[ld_sym_dist_ids$V1 %in% samp_meta$VCF_NAME[c(218, 481, 522, 924)], ]

sort(ld_sym_dist_mat[ld_sym_dist_ids$V1 == 'J515.B', ])[1:10]
```

# Ways to choose problematic samples
# 1) Set arbitrary distance cutoff
# 2) Search by "accession" and look for samples that don't cluster within acceession
# 3) Search by cultivar - find samples that cluster within a cultivar
# 4) Search by ploidy - search for samples that 
# Also try adding a "unified" accession so that all samples that should be
#  from the same accession but different sources can be combined





## Generate all_samps plink objects
* Directory on HA for results
  * `/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps`
* Directory on HA with original VCF files
  * `/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/`

```
bash

source activate plink_1.9_env

# start with single chromosome
VCF_DIR=/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/
VCF_FILE_SHORT=Pvirgatum_1070g_Chr01K.snp.poly.sort.norepeats.vcf.gz
VCF_FILE=$VCF_DIR$VCF_FILE_SHORT

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps

SAMP_FILE=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samp_names_forPlink.txt
CDS_POS_FILE=/home/f2p1/work/grabowsk/data/switchgrass/plink_files/sg_v5_CDS.bed





cd $OUT_DIR

plink --vcf $VCF_FILE --double-id --keep $SAMP_FILE \
--extract 'range' $CDS_POS_FILE



```




