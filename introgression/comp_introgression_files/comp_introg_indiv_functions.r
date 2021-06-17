# Functions for introgression analysis of individuals based on training
#  populations comparisons

gen_geno_state_vec <- function(vcf, samp_name, allele_states){
  ############
  # generate the "genotype state" for each SNP in a sample by pasting together
  #   the allele-state for each allele in the genotype
  ### INPUTS ###
  # vcf = vcf with genotypes
  # samp_name = sample name for analysis; must match up with name used in
  #    vcf
  # allele_states = table with the REF and ALT allele states at each SNP, as
  #    well as the chromosome and position info for each SNP
  ### OUTPUT ###
  # vector with "genotype-state" containing allele state of both alleles 
  #  at each SNP for a sample  
  ############
  geno_vec_1 <- unlist(vcf[, ..samp_name])
  geno_vec_2 <- unlist(lapply(strsplit(geno_vec_1, split = ':'),
    function(x) x[1]))
  #
  state_vec <- rep(NA, times = nrow(vcf))
  tmp_miss_inds <- which(geno_vec_2 == './.')
  state_vec[tmp_miss_inds] <- 'MISSING'
  #
  tmp_homREF_inds <- which(geno_vec_2 == '0/0')
  state_vec[tmp_homREF_inds] <- paste(
    allele_states$REF_STATE[tmp_homREF_inds],
    allele_states$REF_STATE[tmp_homREF_inds], sep = '_')
  #
  tmp_homALT_inds <- which(geno_vec_2 == '1/1')
  state_vec[tmp_homALT_inds] <- paste(
    allele_states$ALT_STATE[tmp_homALT_inds],
    allele_states$ALT_STATE[tmp_homALT_inds], sep = '_')
  #
  tmp_HET_inds <- which(geno_vec_2 == '0/1')
  state_vec[tmp_HET_inds] <- paste(
    allele_states$REF_STATE[tmp_HET_inds],
    allele_states$ALT_STATE[tmp_HET_inds], sep = '_')
  return(state_vec)
}

gen_cat_state_vec <- function(geno_state_vec, pop1_name, pop2_name){
  ################
  # generate geneotype category vector for a sample; categories are defined
  #   elsewhere out of this function
  ### INPUTS ###
  # geno_state_vec = vector generated with 'gen_geno_state_vec' function;
  #      contains allele state for each allele of the genotype at all SNPs
  # pop1_name = the name of pop1 used for labeling columns in introgression 
  #    analysis; pop1 is the majority background population
  # pop2_name = the name of pop2 used for labeling columns in introgression 
  #    analysis; pop2 is the candidate introgression population
  ### OUTPUT ###
  # vector of the category of each genotype for a sample
  ##################
  # find any SNPs that contain at least one allele of the following categories:
  pop2_only_inds <- grep(paste(pop2_name, '_ONLY', sep = ''), tmp_state_vec)
  pop2_mainly_inds <- grep(paste(pop2_name, '_MAINLY', sep = ''), tmp_state_vec)
  pop1_only_inds <- grep(paste(pop1_name, '_ONLY', sep = ''), tmp_state_vec)
  pop1_mainly_inds <- grep(paste(pop1_name, '_MAINLY', sep = ''), tmp_state_vec)
  na_inds <- grep('NO_INFO', tmp_state_vec)
  # assign each SNP to appropriate category
  tmp_cat_vec <- rep(NA, times = nrow(vcf_in))
  tmp_cat_vec[grep(paste(pop1_name, 'ONLY', pop1_name, 'ONLY', sep = '_'),
    tmp_state_vec)] <- 'cat_a'
  tmp_cat_vec[grep(paste(pop1_name, 'MAINLY', pop1_name, 'MAINLY', sep = '_'),
    tmp_state_vec)] <- 'cat_b'
  tmp_cat_vec[intersect(pop1_only_inds, na_inds)] <- 'cat_c'
  tmp_cat_vec[intersect(pop1_mainly_inds, na_inds)] <- 'cat_d'
  #
  tmp_cat_vec[grep(paste(pop2_name, 'ONLY', pop2_name, 'ONLY', sep = '_'),
    tmp_state_vec)] <- 'cat_e'
  tmp_cat_vec[grep(paste(pop2_name, 'MAINLY', pop2_name, 'MAINLY', sep = '_'),
    tmp_state_vec)] <- 'cat_f'
  tmp_cat_vec[intersect(pop2_only_inds, na_inds)] <- 'cat_g'
  tmp_cat_vec[intersect(pop2_mainly_inds, na_inds)] <- 'cat_h'
  #
  tmp_cat_vec[intersect(pop1_only_inds, pop2_only_inds)] <- 'cat_i'
  tmp_cat_vec[intersect(pop1_only_inds, pop2_mainly_inds)] <- 'cat_j'
  tmp_cat_vec[intersect(pop1_mainly_inds, pop2_only_inds)] <- 'cat_k'
  tmp_cat_vec[intersect(pop1_mainly_inds, pop2_mainly_inds)] <- 'cat_l'
  #
  tmp_cat_vec[grep('NO_INFO_NO_INFO', tmp_state_vec)] <- 'cat_m'
  tmp_cat_vec[grep('MISSING', tmp_state_vec)] <- 'cat_n'
  #
  return(tmp_cat_vec)
}

gen_pop2_presence_nonoverlap_window_tab <- function(samp_state_vec, vcf, 
  chrom, window_size){
  ###########
  # within NON-overlapping windows, calculate 1) the number SNPs with any 
  #   pop2 allele; 2) the number of SNPs with pop2_only alleles; and 
  #   3) the number of SNPs with pop2 alleles * number of SNPs with pop2-only
  #   alleles
  ### INPUTS ###
  # samp_state_vec = vector generated with 'gen_geno_state_vec' function;
  #      contains allele state for each allele of the genotype at all SNPs
  # vcf = vcf with genotypes
  # chrom = the chromosome getting analyzed
  # window_size = the number of SNPs in each window
  ### OUTPUT ###
  # data.table with chromosome, position of beginning of each window, 
  #  N_ANY_POP2 = the count of SNPs with any pop2 allele in each window, 
  #  N_POP2_ONLY = the count of SNPs with pop2-only alleles in each window,
  #  ADJ_N_POP2 = N_ANY_POP2 * N_POP2_ONLY (adjusted count, gives more weight
  #    to windows with more pop2_only alleles)
  ################
  chr_inds <- grep(chrom, vcf$CHROM)
  window_start_inds <- chr_inds[seq(1, length(chr_inds), by = window_size)]
  window_start_inds <- window_start_inds[-length(window_start_inds)]
  window_nonoverlap_any_pop2_vec <- c()
  window_nonoverlap_pop2_only_vec <- c()
  for(i in window_start_inds){
    test_inds <- c(i:(i+(window_size - 1)))
    any_pop2_count <- length(grep(pop2_name, samp_state_vec[test_inds]))
    window_nonoverlap_any_pop2_vec <- c(window_nonoverlap_any_pop2_vec,
      any_pop2_count)
    pop2only_count <- length(grep(paste(pop2_name, '_ONLY', sep = ''), 
      samp_state_vec[test_inds]))
    window_nonoverlap_pop2_only_vec <- c(window_nonoverlap_pop2_only_vec,
      pop2only_count)
  }
  adj_n_pop2 <- window_nonoverlap_any_pop2_vec * 
    window_nonoverlap_pop2_only_vec
  nonoverlap_plot_tab <- data.table(CHROM = vcf$CHROM[window_start_inds],
    POS = vcf$POS[window_start_inds], 
    POS_CUM = vcf$POS_CUM[window_start_inds],
    N_ANY_POP2 = window_nonoverlap_any_pop2_vec,
    N_POP2_ONLY = window_nonoverlap_pop2_only_vec,
    ADJ_N_POP2 = adj_n_pop2)
  return(nonoverlap_plot_tab)
}

#test_nonoverlap <- gen_pop2_presence_nonoverlap_window_vec(
#  samp_state_vec = tmp_state_vec, vcf = vcf_in, chrom = 'Chr01K',
#  window_size = 10)

gen_pop2_presence_overlap_window_tab <- function(samp_state_vec, vcf,
  chrom, window_size){
  ###########
  # within OVERLAPPING windows, calculate 1) the number SNPs with any 
  #   pop2 allele and 2) the number of SNPs with pop2_only alleles
  ### INPUTS ###
  # samp_state_vec = vector generated with 'gen_geno_state_vec' function;
  #      contains allele state for each allele of the genotype at all SNPs
  # vcf = vcf with genotypes
  # chrom = the chromosome getting analyzed
  # window_size = the number of SNPs in each window
  ### OUTPUT ###
  # data.table with chromosome, position of beginning of each window, 
  #  N_ANY_POP2 = the count of SNPs with any pop2 allele in each window, 
  #  N_POP2_ONLY = the count of SNPs with pop2-only alleles in each window,
  #  ADJ_N_POP2 = N_ANY_POP2 * N_POP2_ONLY (adjusted count, gives more weight
  #    to windows with more pop2_only alleles)
  ################
  chr_inds <- grep(chrom, vcf$CHROM)
  window_inds <- chr_inds[1:(length(chr_inds)-(window_size - 1))]
  window_overlap_any_pop2_vec <- c()
  window_overlap_pop2_only_vec <- c()
  for(i in window_inds){
    test_inds <- c(i:(i+(window_size - 1)))
    any_pop2_count <- length(grep(pop2_name, samp_state_vec[test_inds]))
    window_overlap_any_pop2_vec <- c(window_overlap_any_pop2_vec,
      any_pop2_count)
    pop2_only_count <- length(grep(paste(pop2_name, '_ONLY', sep = ''), 
      samp_state_vec[test_inds]))
    window_overlap_pop2_only_vec <- c(window_overlap_pop2_only_vec,
      pop2_only_count)
  }
  adj_n_pop2 <- window_overlap_any_pop2_vec * window_overlap_pop2_only_vec
  overlap_plot_tab <- data.table(CHROM = vcf$CHROM[window_inds],
    POS = vcf$POS[window_inds],
    POS_CUM = vcf$POS_CUM[window_inds],
    N_ANY_POP2 = window_overlap_any_pop2_vec,
    N_POP2_ONLY = window_overlap_pop2_only_vec,
    ADJ_N_POP2 = adj_n_pop2)
  return(overlap_plot_tab)
}

#test_overlap <- gen_pop2_presence_overlap_window_tab(
#  samp_state_vec = tmp_state_vec, vcf = vcf_in, chrom = 'Chr01K',
#  window_size = 10)

gen_pop2_score_nonoverlap_window_tab <- function(samp_score_vec, vcf, 
  chrom, window_size){
  ###########
  # within NON-overlapping windows, calculate cumulative "pop2 score";
  #   scores for each genotype class are defined elsewhere and used to make
  #   samp_score_vec
  ### INPUTS ###
  # samp_score_vec = vector of pop2 genotypes scores for each SNP based on 
  #      the genotype states; scores are defined and calculated elsewhere
  # vcf = vcf with genotypes
  # chrom = the chromosome getting analyzed
  # window_size = the number of SNPs in each window
  ### OUTPUT ###
  # data.table with chromosome, position of beginning of each window, 
  #  and "pop2 score" for each window
  ################
  chr_inds <- grep(chrom, vcf$CHROM)
  window_start_inds <- chr_inds[seq(1, length(chr_inds), by = window_size)]
  window_start_inds <- window_start_inds[-length(window_start_inds)]
  window_score_vec <- c()
  for(i in window_start_inds){
    test_inds <- c(i:(i+(window_size - 1)))
    score_count <- sum(samp_score_vec[test_inds])
    window_score_vec <- c(window_score_vec, score_count)
  }
  nonoverlap_score_tab <- data.table(CHROM = vcf$CHROM[window_start_inds],
    POS = vcf$POS[window_start_inds],
    POS_CUM = vcf$POS_CUM[window_start_inds],
    POP2_SCORE = window_score_vec)
  return(nonoverlap_score_tab)
}

gen_pop2_score_overlap_window_tab <- function(samp_score_vec, vcf,
  chrom, window_size){
  ###########
  # within OVERLAPPING windows, calculate cumulative "pop2 score";
  #   scores for each genotype class are defined elsewhere and used to make
  #   samp_score_vec
  ### INPUTS ###
  # samp_score_vec = vector of pop2 genotypes scores for each SNP based on 
  #      the genotype states; scores are defined and calculated elsewhere
  # vcf = vcf with genotypes
  # chrom = the chromosome getting analyzed
  # window_size = the number of SNPs in each window
  ### OUTPUT ###
  # data.table with chromosome, position of beginning of each window, 
  #  and "pop2 score" for each window
  ################
  chr_inds <- grep(chrom, vcf$CHROM)
  window_inds <- chr_inds[1:(length(chr_inds)-(window_size - 1))]
  window_score_vec <- c()
  for(i in window_inds){
    test_inds <- c(i:(i+(window_size - 1)))
    score_count <- sum(samp_score_vec[test_inds])
    window_score_vec <- c(window_score_vec, score_count)
  }
  overlap_score_tab <- data.table(CHROM = vcf$CHROM[window_inds],
    POS = vcf$POS[window_inds],
    POS_CUM = vcf$POS_CUM[window_inds],
    POP2_SCORE = window_score_vec)
  return(overlap_score_tab)
}



