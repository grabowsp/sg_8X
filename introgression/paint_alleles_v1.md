# First attempt at painting chromosomes by alleles

## Steps
1) Find hi-Fst files
2) For each comparison VCF, pull out genepool training-set genotype frequencies
3) Compare training-set genotype frequencies between genepools
4) Assign genepool-diagnistic and/or non-informative states to each allele/genotype
5) Code sample genotypes by allele/genotype designation

## Location of files at HA
* Directory with files
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/`
* Atlantic vs Gulf
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsGulfhiFst.v3.disomic.CDS.vcf.gz`
* Gulf vs Midwest
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMWhiFst.v3.disomic.CDS.vcf.gz`
* Atlantic vs Midwest
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMWhiFst.v3.disomic.CDS.vcf.gz`
* Sample and Group Map
  * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/introgression_samps_and_groups.txt`


## Make genepool-control sample files for VCFtools
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
samp_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/introgression_samps_and_groups.txt'

samp_info <- fread(samp_file)

### SET OUTPUT ###
out_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'


for(pop in unique(samp_info$GROUP)){
  tmp_names <- samp_info[GROUP == pop, list(LIB)]
  tmp_outname <- paste(out_dir, pop, '_libnames.txt', sep = '')
#  print(tmp_outname)
  fwrite(tmp_names, file = tmp_outname, col.names = F, sep = '\t')
}

```

## Extract genotype frequencies from training sets
### MW vs GULF
```
bash
source activate bioinformatics_env

OUTDIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
VCF_SHORT=GulfVsMWhiFst.v3.disomic.CDS.vcf.gz
GROUP_SHORT=Midwest_libnames.txt
OUT_PRE=MW_training_MWvGULF

cd $OUTDIR

VCF_IN=$OUTDIR$VCF_SHORT
GROUP_IN=$OUTDIR$GROUP_SHORT
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --hardy

GROUP_SHORT=Gulf_libnames.txt
OUT_PRE=GULF_training_MWvGULF
VCF_IN=$OUTDIR$VCF_SHORT
GROUP_IN=$OUTDIR$GROUP_SHORT
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN --hardy

```

## Assign allele states
* Allele states
  * genepool-only
  * genepool-likely
  * not informative
### MW vs GULF
* Used R script: `~sg_8X/introgression/mw_v_gulf_states.r`
* Output (on HA):
  * REF allele freq in MW and GULF at each hiFst SNP
    * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_v_GULF_ref_freq.txt`
  * Allele states at each hiFst SNP 
    * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_v_GULF_allele_states.txt`
  * Position file of filtered hiFst SNPs for further VCFtools work
    * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/MW_v_GULF_keep_pos.txt`
#### Notes on allele and genotype coding
```
# Allele coding
# 1 = MW_only_REF
# 2 = MW_mainly_REF
# 3 = MW_only_ALT
# 4 = MW_mainly_ALT

# A = GULF_only_REF
# B = GULF_mainly_REF
# C = GULF_only_ALT
# D = GULF_mainly_ALT

# Y = non-informative_REF
# Z = non-informative_ALT

# Genotype Categories and Color options
# cat_a = 1:1 = 3:3 = dark blue (HOM MW_ONLY)
# cat_b = 2:2 = 4:4 = medium blue (HOM MW_MAINLY)
# cat_c = 1:Z = 3:Y = dark blue, less saturation (HET MW_ONLY:NO_INFO)
# cat_d = 2:Z = 4:Y = medium blue, less saturation (HET MW_MAINLY:NO_INFO)

# cat_e = A:A = C:C = dark red (HOM GULF_ONLY)
# cat_f = B:B = D:D = medium red (HOM GULF_MAINLY)
# cat_g = A:Z = C:Y = dark red, less saturaion (HET GULF_ONLY:NO_INFO)
# cat_h = B:Z = D:Y = medium red, less saturation (HET GULF_MAINLY:NO_INFO)

# cat_i = 1:C = 3:A =  magenta (HET MW_ONLY:GULF_ONLY)
# cat_j = 1:D = 3:B =  dark magenta (more blue) (HET MW_ONLY:GULF_MAINLY)
# cat_k = 2:C = 4:A =  lighter magenta (more red) (HET MW_MAINLY:GULF_ONLY)
# cat_l = 2:D = 4:B =  magenta with less saturation (HET MW_MAINLY:GULF_MAINLY)

# cat_m = Z:Z = Y:Y = grey (HOM NO_INFO)
```

## Generate grp2 MW vs GULF VCF
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
VCF_SHORT=GulfVsMWhiFst.v3.disomic.CDS.vcf.gz
SAMP_SHORT=grp2_libnames.txt
POS_SHORT=MW_v_GULF_keep_pos.txt
OUTFILE=grp2_MWvGULF.vcf
#OUT_PRE=grp2_MWvGULF

cd $OUT_DIR

vcftools --gzvcf $VCF_SHORT -c  --keep $SAMP_SHORT \
--positions $POS_SHORT --recode> $OUTFILE

```

### Try painting grp 2
* `/home/grabowsky/tools/workflows/sg_8X/introgression/mw_v_gulf_grp2_painting.r`

### Calculate allele frequency in grp2
* Steps
  * Get small-scale population splits
  * Find grp2 samples that are part of same small-group
    * old grp2 = MW_01 (35 samples)
  * Generate HW file for grp2 subset
  * Compare grp2 freq to training freqs
* R scratch for checking grp2
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

res_tab_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_result_tabs/natv2filt_res_tab_v4.0.txt'
res_tab <- fread(res_tab_file)

res_tab[grep('MW_', res_tab$sub_grp), .N, by = sub_grp]


vcf_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/', 'grp2_MWvGULF.vcf', sep = '')
vcf_in <- read.table(vcf_file, header = F, stringsAsFactors = F)

vcf_head_tmp <- system(paste('grep CHR ', vcf_file, sep = ''), intern = T)
vcf_head_2 <- sub('#CHROM', 'CHROM', vcf_head_tmp)
vcf_head_3 <- unlist(strsplit(vcf_head_2, split = '\t'))

colnames(vcf_in) <- vcf_head_3
vcf_in <- data.table(vcf_in)







