# Assign allele states to the alleles in introgression comparisons

### LOAD ENVIRONMENT ###
# bash
# source activate R_analysis

# args <- commandArgs(trailingOnly=T)

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
#pop1_hw_file <- args[1]
pop1_hw_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/', 'MW_training_MWvATL.hwe', sep = '')
pop1_hw <- fread(pop1_hw_file)

#pop2_hw_file <- args[2]
pop2_hw_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/', 'ATL_training_MWvATL.hwe', sep = '')
pop2_hw <- fread(pop2_hw_file)

### SET OUTPUT ###
#pop1_name <- args[3]
pop1_name <- 'MW'
#pop2_name <- args[4]
pop2_name <- 'ATL'

#out_dir <- args[5]
out_dir <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/MW_analysis/', sep = '')
tmp_string <- rev(unlist(strsplit(out_dir, split = '')))
if(tmp_string[1] != '/'){
  out_dir <- paste(out_dir, '/', sep = '')
}

allele_freq_file <- paste(out_dir, pop1_name, '_v_', pop2_name, 
  '_ref_freq.txt', sep = '')

keep_position_file <- paste(out_dir, pop1_name, '_v_', pop2_name, 
  '_keep_pos.txt', sep = '')

allele_state_file <- paste(out_dir, pop1_name, '_v_', pop2_name, 
  '_allele_states.txt', sep = '')

### SET VARIABLES ###
#missing_cut <- as.numeric(args[6])
missing_cut <- 0.8
# minimum % of samples with genotypes; max missing = 1-missing_cut

#n_training <- as.numeric(args[7])
n_training <- 40
# the number of samples in the training set

#min_keep_freq <- as.numeric(args[8])
min_keep_freq <- 0.5
# The minimum frequency for an allele to be considered informative in a
#  population

prob_allele_ratio <- 2.5
# The probability ratio of seeing the allele in popA vs popB for
#   the allele to be considered informative for popA ancestry

#####################

length(setdiff(pop1_hw$POS, pop2_hw$POS))
# [1] 0

### Filter out sites with high missing data

pop1_geno_count_pre <- lapply(strsplit(
  unlist(pop1_hw[, c('OBS(HOM1/HET/HOM2)')]), 
  split = '/', fixed = T), function(x) as.numeric(x))

pop1_n_genos <- unlist(lapply(pop1_geno_count_pre, function(x) sum(x)))
pop1_lowcount_inds <- which(pop1_n_genos < (n_training * missing_cut))

pop2_geno_count_pre <- lapply(strsplit(
  unlist(pop2_hw[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

pop2_n_genos <- unlist(lapply(pop2_geno_count_pre, function(x) sum(x)))
pop2_lowcount_inds <- which(pop2_n_genos < (n_training * missing_cut))

lowcount_inds <- sort(union(pop1_lowcount_inds, pop2_lowcount_inds))

pop1_hw_filt <- pop1_hw[-lowcount_inds]
pop2_hw_filt <- pop2_hw[-lowcount_inds]

#######
# tally the number of samples with each genotype at each SNP and the
#   REF and ALT frequency at each SNP in each pop
pop1_geno_count <- lapply(strsplit(
  unlist(pop1_hw_filt[, c('OBS(HOM1/HET/HOM2)')]),    
  split = '/', fixed = T), function(x) as.numeric(x))

pop1_ref_freq <- unlist(lapply(pop1_geno_count, function(x) 
  (x[1] * 2 + x[2])/(sum(x)*2)))
pop1_alt_freq <- unlist(lapply(pop1_geno_count, function(x) 
  (x[3] * 2 + x[2])/(sum(x)*2)))

pop2_geno_count <- lapply(strsplit(
  unlist(pop2_hw_filt[, c('OBS(HOM1/HET/HOM2)')]),
  split = '/', fixed = T), function(x) as.numeric(x))

pop2_ref_freq <- unlist(lapply(pop2_geno_count, function(x) 
  (x[1] * 2 + x[2])/(sum(x)*2)))
pop2_alt_freq <- unlist(lapply(pop2_geno_count, function(x) 
  (x[3] * 2 + x[2])/(sum(x)*2)))

summary(abs(pop2_ref_freq - pop1_ref_freq))

freq_table <- data.table(CHR = pop1_hw_filt$CHR, POS = pop1_hw_filt$POS,
  pop1_ref_freq = pop1_ref_freq,
  pop2_ref_freq = pop2_ref_freq)

pop_freq_names <- paste(c(pop1_name, pop2_name), '_ref_freq', sep = '')
colnames(freq_table)[c(3,4)] <- pop_freq_names

fwrite(freq_table, file = allele_freq_file)

###
# calculate the probabilty of detecting the REF of ALT allele at least
#  once in a sample given the training sets allele frequencies
pop1_prob_ref_allele <- sapply(pop1_ref_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))
pop1_prob_alt_allele <- sapply(pop1_alt_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))

pop2_prob_ref_allele <- sapply(pop2_ref_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))
pop2_prob_alt_allele <- sapply(pop2_alt_freq, function(x)
  pbinom(1, size = 2, prob = (1-x)))

pop2_only_REF <- which(pop1_ref_freq == 0)
pop1_only_REF <- which(pop2_ref_freq == 0)
pop1_mainly_REF <- which(
  pop1_ref_freq > min_keep_freq &
  (pop1_prob_ref_allele/pop2_prob_ref_allele) > prob_allele_ratio &
  pop2_ref_freq > 0)
pop2_mainly_REF <- which(
  pop2_ref_freq > min_keep_freq &
  (pop2_prob_ref_allele / pop1_prob_ref_allele) > prob_allele_ratio &
  pop1_ref_freq > 0)
tmp_REF_tot <- sort(c(pop2_only_REF, pop1_only_REF, 
  pop1_mainly_REF, pop2_mainly_REF))
noinfo_REF <- setdiff(seq(length(pop1_ref_freq)), tmp_REF_tot)

length(tmp_REF_tot)
# 28655
sum(duplicated(tmp_REF_tot))
# 0
length(noinfo_REF)
# [1] 30458

pop2_only_ALT <- which(pop1_alt_freq == 0)
pop1_only_ALT <- which(pop2_alt_freq == 0)
pop1_mainly_ALT <- which(
  pop1_alt_freq > min_keep_freq &
  (pop1_prob_alt_allele/pop2_prob_alt_allele) > prob_allele_ratio &
  pop2_alt_freq > 0)
pop2_mainly_ALT <- which(
  pop2_alt_freq > min_keep_freq &
  (pop2_prob_alt_allele / pop1_prob_alt_allele) > prob_allele_ratio &
  pop1_alt_freq > 0)
tmp_ALT_tot <- sort(c(pop2_only_ALT, pop1_only_ALT,
  pop1_mainly_ALT, pop2_mainly_ALT))
noinfo_ALT <- setdiff(seq(length(pop1_alt_freq)), tmp_ALT_tot)

length(tmp_ALT_tot)
# 50983
sum(duplicated(tmp_ALT_tot))
# 0
length(noinfo_ALT)
# 8130 

tmp_BOTH <- c(tmp_REF_tot, tmp_ALT_tot)
length(setdiff(seq(length(pop1_alt_freq)), tmp_BOTH))
# 0 not included in any category if use prob_allele_ratio = 2.5
# 42 not included in any category if use prob_allele_ratio = 2.67

sum(duplicated(tmp_BOTH))
# 20525 SNPs informative for both alleles

length(c(pop1_only_REF, pop1_mainly_REF, pop1_only_ALT, pop1_mainly_ALT))
# 42195 SNPs informative about MW
length(c(pop2_only_REF, pop2_mainly_REF, pop2_only_ALT, pop2_mainly_ALT))
# 37443 SNPs informative about ATL

allele_assign_tab <- data.table(CHR = pop1_hw_filt$CHR,
  POS = pop1_hw_filt$POS,
  REF_STATE = as.character(NA),
  ALT_STATE = as.character(NA))

allele_assign_tab[pop2_only_REF, REF_STATE := paste(pop2_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop1_only_REF, REF_STATE := paste(pop1_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop2_mainly_REF, REF_STATE := paste(pop2_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[pop1_mainly_REF, REF_STATE := paste(pop1_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[noinfo_REF, REF_STATE := 'NO_INFO']

allele_assign_tab[pop2_only_ALT, ALT_STATE := paste(pop2_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop1_only_ALT, ALT_STATE := paste(pop1_name, '_ONLY', 
  sep = '')]
allele_assign_tab[pop2_mainly_ALT, ALT_STATE := paste(pop2_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[pop1_mainly_ALT, ALT_STATE := paste(pop1_name, '_MAINLY', 
  sep = '')]
allele_assign_tab[noinfo_ALT, ALT_STATE := 'NO_INFO']

fwrite(allele_assign_tab, file = allele_state_file,
  sep = '\t')

fwrite(allele_assign_tab[, c('CHR', 'POS')], 
  file = keep_position_file, sep = '\t')

quit(save = 'no')

