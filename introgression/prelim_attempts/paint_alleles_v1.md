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

## Calculate allele frequency in MW_03 (most of old grp2)
* Steps
  * Get small-scale population splits
  * Find grp2 samples that are part of same small-group
    * old grp2 = MW_03 (52 samples), MW_02 (1 weird 8X), MW_04 (7), and MW_05 (6)
  * Generate HW file for grp2 subset
  * Compare grp2 freq to training freqs

### Generate HW file for MW_03 for MW vs GULF
```
### MW vs GULF
```
bash
source activate bioinformatics_env

OUTDIR=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/
VCF_SHORT=GulfVsMWhiFst.v3.disomic.CDS.vcf.gz
VCF_IN=$OUTDIR$VCF_SHORT

POS_SHORT=MW_v_GULF_keep_pos.txt
POS_FILE=$OUTDIR$POS_SHORT

GROUP_DIR=/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/natv2filt_subgrps/

cd $OUTDIR

GROUP_SHORT=MW_03_names.txt
GROUP_IN=$GROUP_DIR$GROUP_SHORT
OUT_PRE=MW_03_MWvGULF
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN \
--positions $POS_FILE--freq
vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $GROUP_IN \
--positions $POS_FILE --hardy

```
#### Look at MW_03 allele frequencies
```
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
data_dir <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/'

subgrp_hw_file <- paste(data_dir, 'MW_03_MWvGULF.hwe', sep = '')
subgrp_hw <- fread(subgrp_hw_file)

train_freq_file <- paste(data_dir, 'MW_v_GULF_ref_freq.txt', sep = '')
train_freq <- fread(train_freq_file)

allele_state_file <- paste(data_dir, 'MW_v_GULF_allele_states.txt', sep = '')
allele_states <- fread(allele_state_file)



##########

subgrp_geno_count <- lapply(strsplit(
  unlist(subgrp_hw[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

subgrp_n_genos <- unlist(lapply(subgrp_geno_count, function(x) sum(x)))
subgrp_low_inds <- which(subgrp_n_genos < (max(subgrp_n_genos)*0.8))
# [1] 667

subgrp_ref_freq <- unlist(lapply(subgrp_geno_count, function(x)
  (x[1] * 2 + x[2])/(sum(x)*2)))
subgrp_alt_freq <- unlist(lapply(subgrp_geno_count, function(x)
  (x[3] * 2 + x[2])/(sum(x)*2)))

###

# calculate the REF freq difference between pop2 and pop1 and convert to 
#  their range
ref_freq_dif <- unlist(train_freq[,4] - train_freq[,3])
ref_freq_range <- abs(ref_freq_dif)

# find SNPs where pop1 has higher REF freq
pop1_ref_hi <- which(ref_freq_dif < 0)

# calc diff in REF freq between subgrp and pop1
sub_v_pop1_ref <- subgrp_ref_freq - unlist(train_freq[,3])

# adjust diff for SNPs where pop1 has higher REF freq than pop2 because
#  want this to represent difference in direction of pop2
sub_v_pop1_ref_2 <- sub_v_pop1_ref
sub_v_pop1_ref_2[pop1_ref_hi] <- sub_v_pop1_ref[pop1_ref_hi]*-1

F_sub_v_pop1 <- sub_v_pop1_ref_2 / ref_freq_range
F_sub_v_pop1[subgrp_low_inds] <- NA

# CONTINUE FROM HERE WITH WINDOW ANALYSIS
### decide on window cutoff of F_sub_v_pop1

```

## Notes from R when trying to actually paint chromosomes
```
# My feeling about individual colors - it's going to be too hard to make
#  a meaningful plot where colors show good resolution. I think, for now,
#  the line plots will show the best pattern. Perhaps we can show color
#  patterns for close-ups of certain regions...

MW_G_palette <- colorRampPalette(brewer.pal(9, 'Oranges'))
mwg_cont <- scale_colour_gradientn(colours = MW_G_palette(100),
  limits = c(0.5,20))

score_color_names <- as.character(seq(from=0, to = 20, by = 0.5))
gray_palette <- rev(gray.colors(n = length(score_color_names), start = 0.2,
  end = 1.0))
names(gray_palette) <- score_color_names

mwg_cont <- scale_colour_gradientn(colours = MW_G_palette(100),
  limits = c(0,20))

scale_color_grey(start = 0.2, end = 1.0, limits = c(0,20))

tmp_tab <- tot_pop2_list[[1]][[4]]
tmp_tab[, PLOT_COL := as.character(NA)]
for(i in unique(tmp_tab$POP2_SCORE)){
  tmp_inds <- which(tmp_tab$POP2_SCORE == i)
  tmp_tab[tmp_inds, PLOT_COL := gray_palette[as.character(i)]]
}
gg_col_test_2 <- ggplot(tmp_tab) +
  theme_void() +
#  geom_point() + 
  xlim(-10000, (max(tmp_tab$POS_CUM)+10000)) +
  ylim(0, 1) +
  geom_vline(xintercept = tmp_tab$POS_CUM, color = tmp_tab$PLOT_COL)
#  mwg_cont
#  scale_color_grey(start = 0.2, end = 1.0, limits = c(0,20))
#  scale_colour_gradient(low = 'grey20', high = 'black', limits = c(0,20))

test_color_plot_2 <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/J567.A_Chr01K_introgression_color_2.pdf'

pdf(test_color_plot_2, width = 18, height = 4)
gg_col_test_2
dev.off()

tmp_geno_cat_vec <- tot_pop2_list[[1]][['geno_cat_vec']]

geno_color_vec <- rep(NA, times = length(tmp_geno_cat_vec))
for(cn in unique(tmp_geno_cat_vec)){
  geno_color_vec[which(tmp_geno_cat_vec == cn)] <- intro_col_vec[[cn]]
}

tmp_tab <- data.table(CHROM = vcf_in$CHROM, POS_CUM = vcf_in$POS_CUM,
  PLOT_COLOR = geno_color_vec)

test_color_plot_1 <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/J567.A_Chr01K_introgression_color_1.pdf'

```
### Notes about colors
```
# Set color variables
intro_col_vec <- c()
intro_col_vec$cat_a <- 'blue4'
intro_col_vec$cat_b <- 'blue2'
intro_col_vec$cat_c <- 'steelblue3'
intro_col_vec$cat_d <- 'steelblue1'
intro_col_vec$cat_e <- 'red3'
intro_col_vec$cat_f <- 'red1'
intro_col_vec$cat_g <- 'indianred3'
intro_col_vec$cat_h <- 'indianred1'
intro_col_vec$cat_i <- 'mediumorchid3'
intro_col_vec$cat_j <- 'purple3'
intro_col_vec$cat_k <- 'magenta3'
intro_col_vec$cat_l <- 'mediumorchid1'
intro_col_vec$cat_m <- 'gray75'
intro_col_vec$cat_n <- 'white'


```


