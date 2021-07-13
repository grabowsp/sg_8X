# Select SNPs to be used as control for RDA

# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
rda_snp_info_file <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_ATL_into_MW_Fst_0.5.RDA_SNP_info.txt'
rda_snp_info <- fread(rda_snp_info_file)

mwall_rda_snp_freq_file <- paste('/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files/', 
  'AtlanticvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.MWall.frq', sep = '')
mwall_rda_snp_freq <- fread(mwall_rda_snp_freq_file)
colnames(mwall_rda_snp_freq) <- c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 
  'ALLELE1_FREQ', 'ALLELE2_FREQ')


mwall_tot_freq_dir <- '/home/f1p1/tmp/switchgrass_8X/firstdraft_introg_results/MWall_freq_files/'
mwall_tot_freq_files <- system(paste('ls ', mwall_tot_freq_dir, 
  'Chr*MWall.frq', sep = ''), intern = T)

### SET OUTPUTS ###
out_dir <- '/home/f1p1/tmp/switchgrass_8X/MWall_control_snps/'

out_file_pre <- 'ATLvMW_Fst_0.5_MWall_RDA_control_snps'

### SET VARIABLES ###

maf_cut <- 0.05

n_control_sets <- 10

###########

# Generate MAF list for all SNPs
mwall_tot_freq_list <- list()
for(i in seq(length(mwall_tot_freq_files))){
  mwall_tot_freq_list[[i]] <- fread(mwall_tot_freq_files[i])
}

tot_freq <- rbindlist(mwall_tot_freq_list)
colnames(tot_freq) <- c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 
  'ALLELE1_FREQ', 'ALLELE2_FREQ')

af_1 <- as.numeric(unlist(
  lapply(strsplit(tot_freq$ALLELE1_FREQ, split = ':'), function(x) x[2])))
af_2 <- as.numeric(unlist(
  lapply(strsplit(tot_freq$ALLELE2_FREQ, split = ':'), function(x) x[2])))
af_tab <- data.table(af_1 = af_1, af_2 = af_2)
maf_vec <- apply(af_tab, 1, function(x) x[which.min(x)])

tot_freq$MAF <- maf_vec

tot_freq_1 <- tot_freq[MAF >= maf_cut, ]

######

choose_sub_quantile_snps <- function(test_freq_tab, tot_freq_tab, goal_n, 
  n_quantiles){
  # Generate random set of SNPs with same MAF distribution as a test set
  # INPUTS
  # test_freq_tab = SNP frequency data.table for test SNPs, must include CHROM,
  #                  POS, and MAF columns
  # tot_freq_tab = SNP frequency data.table for all SNPs to be sampled from;
  #                  must include CHROM, POS, and MAF columns
  # goal_n = the goal number of SNPs to sample; the total number will
  #           usually be a little more because equal numbers of SNPs are
  #           sampled from each quantile
  # n_quantiles = the number of quantiles to sample MAF's from;
  # OUTPUT
  # data.table of sampled SNPs in same format as tot_freq_tab
  #########
  # number of SNPs to choose for each quantile
  n_quant <- ceiling(goal_n/n_quantiles)
  # quantile vector from test SNPs
  tmp_qnt <- quantile(test_freq_tab[, MAF], prob = seq(0,1,(1/n_quantiles)))
  #
  # remove test SNPs from total SNP list, if they are present
  test_snp_names <- paste(test_freq_tab$CHROM, test_freq_tab$POS, sep = '_')
  tot_snp_names <- paste(tot_freq_tab$CHROM, tot_freq_tab$POS, sep = '_')
  overlap_inds <- which(tot_snp_names %in% test_snp_names)
  if(length(overlap_inds) > 0){
    tot_freq_tab <- tot_freq_tab[-overlap_inds, ]
  }
  #
  subsamp_list <- list()
  for(i in seq(n_quantiles)){
    tmp_tot_qnt_inds <- sample(which(tot_freq_tab$MAF >= tmp_qnt[i] &
      tot_freq_tab$MAF <= tmp_qnt[(i+1)]), size = n_quant)
    subsamp_list[[i]] <- tmp_tot_qnt_inds
  }
  tot_sub_inds <- sort(unlist(subsamp_list))
  sub_tab <- tot_freq_tab[tot_sub_inds, ]
  return(sub_tab)
}

###

# Generate random SNP list

tot_tab_list <- list()
for(TOT_REP in seq(n_control_sets)){
  test_snp_tab <- table(rda_snp_info$CHROM)
  tot_sub_list <- list()
  for(j in seq(length(test_snp_tab))){
    chr_name <- names(test_snp_tab)[j]
    n_goal <- test_snp_tab[j]
    tot_sub_list[[chr_name]] <- choose_sub_quantile_snps(
      test_freq_tab = rda_snp_info[CHROM == chr_name],
      tot_freq_tab = tot_freq_1[CHROM == chr_name], goal_n = n_goal, 
      n_quantiles = 10)
  }
  tot_sub_tab <- rbindlist(tot_sub_list)
  tot_tab_list[[TOT_REP]] <- tot_sub_tab
}

for(outf in seq(length(tot_tab_list))){
  tmp_out <- paste(out_dir, out_file_pre, '.', outf, '.txt', sep = '')
  fwrite(tot_tab_list[[outf]][, c('CHROM', 'POS')], tmp_out, col.names = F,
    sep = '\t')
}

