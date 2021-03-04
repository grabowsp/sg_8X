# Steps for generating allsamps VCF that includes ancestral genotypes

## Steps
1. Make new directory
1. Make tetraploid VCF (from original VCF) using 'natv2' and SNP positions of the ancestral.VCF files
2. Copy ancestral VCFs to working directory
2. Sort, bgzip, and tabix the sample and ancestral VCF files prior to merging
3. Merge VCFs

## Location of Files
* Result directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs`
* Original VCF directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs`
* natv2 samp name file
  * `/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt`
* ancestral VCF snp positions
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/$CHR_NAME'.ancestral_AP13_kept_SNP_positions.txt'`
* ancestral tetrasomic VCF
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/$CHR_NAME'Chr09N.ancestral_AP13_tet_format.vcf'`

## Generate, sort, bgzip, and tabix sample VCFs
### Submit jobs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs

sbatch gen_tmp_natv2_ancest_pos_vcfs_Chr01.sh
sbatch gen_tmp_natv2_ancest_pos_vcfs_Chr02_05.sh
sbatch gen_tmp_natv2_ancest_pos_vcfs_Chr06_09.sh
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

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs/
ANCESTRAL_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/
ANCESTRAL_SNP_SUF=.ancestral_AP13_kept_SNP_positions.txt
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/

SAMP_FILE=/global/homes/g/grabowsp/data/switchgrass/metadata_8X/nat_filt_v2_names.txt
SAMPSET=natv2

SORT_SUF='.'$SAMPSET'.ancest_snp_pos.tet.sort.vcf.gz'

cd $OUT_DIR

for CHROM in 01K 01N;
  do
  SORT_OUT=Chr$CHROM$SORT_SUF;
  vcftools --gzvcf \
    $DATA_DIR'Pvirgatum_1070g_Chr'$CHROM'.5alleles.snp.sort.norepeats.vcf.gz' \
    --stdout --positions $ANCESTRAL_DIR'Chr'$CHROM$ANCESTRAL_SNP_SUF \
    --keep $SAMP_FILE --recode --recode-INFO-all | \
  bcftools sort - -Oz -o $SORT_OUT;
  tabix -p vcf $SORT_OUT;
  done;
```

## Sort, bgzip, and tabix ancestral VCFs
* saved into working directory
### Submit script
* ended up just running interactively - is pretty fast
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs 

sbatch gen_ancest_sort_vcf_Chr01_09.sh
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

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/
OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/
ANCEST_FILE_SUF=.ancestral_AP13_tet_format.vcf

cd $OUT_DIR

for CHROM_NUM in {01..09};
do
for CHROM_LET in K N;
  do
  CHROM=$CHROM_NUM$CHROM_LET;
  ANCEST_FILE=$DATA_DIR'Chr'$CHROM$ANCEST_FILE_SUF;
  SORT_OUT=Chr$CHROM'.ancestral.tet.sort.vcf.gz';
  bcftools sort $ANCEST_FILE -Oz -o $SORT_OUT;
  tabix -p vcf $SORT_OUT;
  done;
done

```

## Manually merge sample and ancestral VCFs
* Can't merge using bcftools because not using standard genotype notation
* Need to manually paste columns together
### Extract header without sample name row
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs

gunzip -c Chr01K.natv2.ancest_snp_pos.tet.sort.vcf.gz | head -633 \
> tet.format.vcf.nosamp.header.txt
```
### Merge sample and ancestral files
* Elements of script
  * check that sample and ancetral VCFs contain same SNPs
  * remove headers
  * paste columns, including sample name row
  * add header
  * gzip
* Submit job
  * manually ran commands for Chr01K and Chr01N during testing
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs

sbatch merge_natv2_tet_vcfs_02_09.sh
```
* Example of script
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

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs

cd $DATA_DIR

ANCEST_SUF=.ancestral.tet.sort.vcf.gz
SAMP_VCF_SUF=.natv2.ancest_snp_pos.tet.sort.vcf.gz
OUT_SUF=.natv2.ancestral.tet.vcf.gz

for CHR_NUM in {02..09};
do
for CHR_LET in K N;
do
  CHROM=Chr$CHR_NUM$CHR_LET;
  echo 'Start' $CHROM;
  #
  gunzip -c $CHROM$ANCEST_SUF | awk '/^Chr/ {print $1,$2}' > tmp_check_1.txt;
  gunzip -c $CHROM$SAMP_VCF_SUF | awk '/^Chr/ {print $1,$2}' > tmp_check_2.txt;
  cmp --silent tmp_check_1.txt tmp_check_2.txt && \
  echo 'SUCCESS: Files are identical' || echo 'WARNING: files are different';
  rm tmp_check_1.txt;
  rm tmp_check_2.txt;
  #
  gunzip -c $CHROM$ANCEST_SUF | tail -n +634 | cut -f 10-11 > tmp_merge_1.txt;
  gunzip -c $CHROM$SAMP_VCF_SUF | tail -n +634 > tmp_merge_2.txt;
  paste -d '\t' tmp_merge_2.txt tmp_merge_1.txt > tmp_merge_3.txt;
  cat tet.format.vcf.nosamp.header.txt tmp_merge_3.txt | gzip \
  > $CHROM$OUT_SUF;
  #
  rm tmp_merge_1.txt;
  rm tmp_merge_2.txt;
  rm tmp_merge_3.txt;
  #
  echo 'Finished' $CHROM;
  done;
done

```

## Make VCF with 100k SNPs
* Final file
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/GW.100k.natv2.ancestral.vcf`
    * 98,482 SNPs
### Select subset of SNPs
* `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/100k_sub_snp_pos.txt`
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

# in R

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
tot_snp_pos_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/sg_8x_ancestral_snp_pos.txt'

tot_snp_info <- fread(tot_snp_pos_file)

### SET OUTPUT ###
out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/100k_sub_snp_pos.txt'

#############

sub_inds <- sort(sample(seq(nrow(tot_snp_info)), size = 1e5))

sub_snp_info <- tot_snp_info[sub_inds, ]

fwrite(sub_snp_info, file = out_file, col.names = F, sep = '\t')

quit(save = 'no')
```

### Generate subset VCFs
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
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.10
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/
IN_VCF_SUF=.natv2.ancestral.tet.vcf.gz

SNP_SUB_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/100k_sub_snp_pos.txt

OUT_SUF=.natv2.ancestral.sub100k.vcf

cd $DATA_DIR

for CHR_NUM in {01..09};
do
for CHR_LET in K N;
do
  CHROM=$CHR_NUM$CHR_LET;
  echo $CHROM;
  VCF_OUT=Chr$CHROM$OUT_SUF;
  VCF_IN=$DATA_DIR'Chr'$CHROM$IN_VCF_SUF;
  vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
    --recode --recode-INFO-all > $VCF_OUT;
  done;
done

```
### Concatenate SNPs
* Need to do manually since using non-typical genotype notation
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs

TOT_VCF=GW.100k.natv2.ancestral.vcf

cat Chr01K.natv2.ancestral.sub100k.vcf > $TOT_VCF

cat Chr01N.natv2.ancestral.sub100k.vcf | tail -n +635 >> $TOT_VCF

for CHR_NUM in {02..09};
do
for CHR_LET in K N;
do
  CHROM=Chr$CHR_NUM$CHR_LET;
  echo $CHROM;
  cat $CHROM.natv2.ancestral.sub100k.vcf | tail -n +635 >> $TOT_VCF;
done;
done

```

## Generate Concatenated Sub-SNP VCFs to test reproducibility of tree
### Steps
* Generate 100 SNP positions (100k randomly chosen positions)
* In loop:
  * Make sub-vcf (merged chromosome-wide VCFs already exist)
  * If Chr01K, then start new sub-VCF
  * If other chromosome, then remove header and concatenate to sub-VCF
### Make SNP position files
* Directory with files:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/snp_subset_files`
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

# in R

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
tot_snp_pos_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/sg_8x_ancestral_snp_pos.txt'

tot_snp_info <- fread(tot_snp_pos_file)

### SET OUTPUT ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/snp_subset_files/'

out_pre <- 'natv2.anc.'

### SET VARIABLES ###
n_reps <- 1000

n_snps <- 1e5
n_snp_char <- paste(n_snps/1000, 'k.subsnps.', sep = '')
#############

for(i in seq(n_reps)){
  sub_inds <- sort(sample(seq(nrow(tot_snp_info)), size = n_snps))
  sub_snp_info <- tot_snp_info[sub_inds, ]
  n_rep_char <- sprintf("%04d", i)
  out_full <- paste(out_dir, out_pre, n_snp_char, n_rep_char, '.pos.txt', 
    sep = '')
  fwrite(sub_snp_info, file = out_full, col.names = F, sep = '\t')
}

quit(save = 'no')

```
### Generate sub-VCFs
#### Run script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs

sbatch gen_natv2_ancest_subfiles_01_10.sh
sbatch gen_natv2_ancest_subfiles_11_50.sh
sbatch gen_natv2_ancest_subfiles_51_100.sh

```
#### Example script
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

IN_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/

IN_SUF=.natv2.ancestral.tet.vcf.gz

SNP_POS_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/snp_subset_files/

SNP_POS_PRE=natv2.anc.100k.subsnps.

NSNPS=100k

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/natv2_ancestral_tet_vcfs/sub100k_vcfs/

OUT_PRE=GW.natv2.ancestral.tet.

cd $OUT_DIR

for SNPSET in {0001..0010};
  do
  VCF_OUT=$OUT_PRE$NSNPS'.'$SNPSET'.vcf';
  SNP_SUB_FILE=$SNP_POS_DIR$SNP_POS_PRE$SNPSET'.pos.txt';
  #
  CHROM=Chr01K;
  echo $CHROM $SNPSET;
  VCF_IN=$IN_DIR$CHROM$IN_SUF;
  vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
    --recode --recode-INFO-all > $VCF_OUT;
  #
  CHROM=Chr01N;
  echo $CHROM $SNPSET;
  VCF_IN=$IN_DIR$CHROM$IN_SUF;
  vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
    --recode --recode-INFO-all | tail -n +635 >> $VCF_OUT;
  #
  for CHR_NUM in {02..09};
    do
    for CHR_LET in K N;
      do
      CHROM=Chr$CHR_NUM$CHR_LET;
      echo $CHROM;
      VCF_IN=$IN_DIR$CHROM$IN_SUF;
      vcftools --gzvcf $VCF_IN --stdout --positions $SNP_SUB_FILE \
        --recode --recode-INFO-all | tail -n +635 >> $VCF_OUT;
      done;
    done;
  gzip $VCF_OUT;
  done;
```




## Make NEXUS file - this is from testing before - not sure will adapt...
* `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs/GW.allsamps.ancestral.sub.min4.nexus`
```
module load python/3.7-anaconda-2019.10

DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_ancestral_dis_vcfs/

cd $DATA_DIR

VCF_IN=GW.allsamps.ancestral.sub.vcf.gz

/global/homes/g/grabowsp/tools/vcf2phylip/vcf2phylip.py -i $VCF_IN \
-o ANCESTRAL_DIP_GENO -p -n
```
