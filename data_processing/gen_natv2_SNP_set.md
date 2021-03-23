# Generate List of CDS SNPs to be used for natv2 analysis

## Steps/Overview
* Start with disomic genotypes
* Generate disomic VCF for natv2 samples and MAF of 5+ copies across all samples
  * natv2_dis_tmp vcfs
  * 5/(706*2) = 0.00354
* Generate lists of samples that fall into the following categories
  * MW (>0.9) 4X
  * GULF (>0.9) 4X
  * ATLANTIC (>0.9) 4X
  * MW (>0.9) 8X
* Use VCF tools and --hardy flag to output the number of HETs at each SNP 
for each group
* Identify SNPs with f(HET) > HW maximum
  * 0.5 for 4X
  * 0.875 for 8X
* Use VCF tools and --kept-sites and --exclude-postions to generate list of SNPs

## Location of files
* Result directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering`
* Chromosome-specific SNPs
  * ex: `Chr01K.disomic.CDS.natv2.filt.kept.sites`
* Genome-wide file of kept SNPs
  * `GW.disomic.CDS.natv2.filt.keptSNPs.txt`
* Data directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/`
* Sample file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt`

## Generate natv2 diploid CDS temporary VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering

sbatch gen_tmp_natv2_dip_vcfs_01_05.sh
sbatch gen_tmp_natv2_dip_vcfs_06_09.sh
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
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/
BED_FILE=/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/
SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt
SAMPSET=natv2

cd $OUT_DIR

for CHROM_NUM in {01..05};
do
for CHROM_LET in K N;
  do
  CHROM=Chr$CHROM_NUM$CHROM_LET;
  vcftools --gzvcf \
    $DATA_DIR'Pvirgatum_1070g_'$CHROM'.snp.sort.norepeats.vcf.gz' \
    --max-missing 0.8 \
    --stdout --bed $BED_FILE --keep $SAMP_FILE --recode --recode-INFO-all | \
    gzip -c > $OUT_DIR$CHROM'.disomic.CDS.'$SAMPSET'.vcf.gz';
  done;
done

```

## Generate lists of samples in each gene pool
```
### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
admix_k3_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/natv2_admix/GW_50k_natv2.3.results.txt'
admix_k3_res <- fread(admix_k3_file)
# $V2 = Atlantic, $V3 = MW, $V4 = Gulf

samp_meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

### SET OUTPUTS ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/'

atl_4X_out <- paste(out_dir, 'atlantic_4X_names.txt', sep = '')
gulf_4X_out <- paste(out_dir, 'gulf_4X_names.txt', sep = '')
mw_4X_out <- paste(out_dir, 'midwest_4X_names.txt', sep = '')
mw_8X_out <- paste(out_dir, 'midwest_8X_names.txt', sep = '')

#############

atl_tot_names <- admix_k3_res$V1[which(admix_k3_res$V2 > 0.9)]
mw_tot_names <- admix_k3_res$V1[which(admix_k3_res$V3 > 0.9)]
gulf_tot_names <- admix_k3_res$V1[which(admix_k3_res$V4 > 0.9)]

atl_4X <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '4X'],
  atl_tot_names)
mw_4X <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '4X'],
  mw_tot_names)
gulf_4X <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '4X'],
  gulf_tot_names)
mw_8X <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '8X'],
  mw_tot_names)

fwrite(data.table(name = atl_4X), file = atl_4X_out, col.names = F)
fwrite(data.table(name = mw_4X), file = mw_4X_out, col.names = F)
fwrite(data.table(name = gulf_4X), file = gulf_4X_out, col.names = F)
fwrite(data.table(name = mw_8X), file = mw_8X_out, col.names = F)

```

## Calculate number of HET genotypes at each SNP
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/

sbatch ATL4X_hwe_sub.sh

sbatch GULF4X_hwe_sub.sh
sbatch MW4X_hwe_sub.sh
sbatch MW8X_hwe_sub.sh
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
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/
SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/atlantic_4X_names.txt
SAMPSET=natv2
POPNAME=ATL4X

cd $OUT_DIR

for CHR_NUM in 02..09;
  do
  for CHR_LET in K N;
    do
    CHROM=Chr$CHR_NUM$CHR_LET;
    OUT_PRE=$CHROM'.'$POPNAME;
    vcftools --gzvcf \
      $OUT_DIR$CHROM'.disomic.CDS.natv2.vcf.gz' \
      --out $OUT_PRE --keep $SAMP_FILE --hardy;
    done;
  done

```

## Identify SNPs with excessive HET genotypes
```
### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
res_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/'

### SET OUTPUTS ###
out_full <- paste(res_dir, 'tot_hiHET_SNPs.txt', sep = '')

### SET VARIABLES ###
pop_names <- c('ATL4X', 'GULF4X', 'MW4X', 'MW8X')

max_f_HET_vec <- c(0.5, 0.5, 0.5, 0.875)
names(max_f_HET_vec) <- pop_names

#########

pop_hiHET_list <- list()
for(popname in pop_names){
  print(popname)
  tmp_files <- system(paste('ls ', res_dir, '*', popname, '.hwe', sep = ''), 
    intern = T)
  tmp_hiHet_list <- list()
  for(i in seq(length(tmp_files))){
    print(i)
    tmp_res <- fread(tmp_files[i])
    tmp_counts <- strsplit(unlist(tmp_res[, 3]), split = '/', fixed = T)
    tmp_counts <- lapply(tmp_counts, function(x) as.numeric(x))
    tmp_freq <- unlist(lapply(tmp_counts, function(x) x[2]/sum(x)))
    tmp_hi_inds <- which(tmp_freq > max_f_HET_vec[popname])
    tmp_het_SNPs <- tmp_res[tmp_hi_inds, c(1,2)]
    tmp_hiHet_list[[i]] <- tmp_het_SNPs
  }
  tmp_pop_tab <- rbindlist(tmp_hiHet_list)
  pop_hiHET_list[[popname]] <- tmp_pop_tab
}

lapply(pop_hiHET_list, nrow)
# $ATL4X
# [1] 36526

# $GULF4X
# [1] 24993

# $MW4X
# [1] 48264

# $MW8X
# [1] 6481

pop_hiHET_tab <- rbindlist(pop_hiHET_list)
pop_hiHET_tab[, NAME := paste(pop_hiHET_tab$CHR, pop_hiHET_tab$POS, sep = '_')]

table(table(pop_hiHET_tab$NAME))
#     1     2     3     4 
# 41217 13002 10029  4739 

dup_inds <- which(duplicated(pop_hiHET_tab$NAME))

tmp_tab <- pop_hiHET_tab[-dup_inds, ]
table(table(tmp_tab$NAME))
#     1 
# 68987 

unique_hiHET_SNPs <- pop_hiHET_tab[-dup_inds, c('CHR', 'POS')]

unique_hiHET_SNPs_2 <- unique_hiHET_SNPs[order(unique_hiHET_SNPs$POS)]
unique_hiHET_SNPs_3 <- unique_hiHET_SNPs_2[order(unique_hiHET_SNPs_2$CHR)]

fwrite(unique_hiHET_SNPs_3, file = out_full, row.names = F, col.names = F, 
  sep = '\t')
```

## Get list of SNPs to keep
* Ended up just running it interactively
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
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/
HI_HET_SNP_FILE=$OUT_DIR'tot_hiHET_SNPs.txt'

cd $OUT_DIR

for CHROM_NUM in {03..09};
do
for CHROM_LET in K N;
  do
  CHROM=Chr$CHROM_NUM$CHROM_LET;
  OUT_PRE=$CHROM'.disomic.CDS.natv2.filt'
  echo $CHROM;
  vcftools --gzvcf \
    $DATA_DIR$CHROM'.disomic.CDS.natv2.vcf.gz' \
    --maf 0.00354 --exclude-positions $HI_HET_SNP_FILE \
    --kept-sites --out $OUT_PRE;
  done;
done

```
### Generate GW list
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_snp_filtering/

OUT_FILE=GW.disomic.CDS.natv2.filt.keptSNPs.txt

for KEPT in *kept.sites;
  do
  tail -n +2 $KEPT >> $OUT_FILE;
  done;

```


