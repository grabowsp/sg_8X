# Functions to be used for calculated numbers of private/shared alleles, etc.

gen_admix_popind_list <- function(admix_res, group_cutoff){
  # Get sample indices for samples with ancestry above "group" cutoff for
  #   each population in admix_res
  # INPUTS
  # admix_res = data.table or data.frame with first column = sample names and
  #              remaining columns from Q matrix outputted by ADMIXTURE
  # group_cutoff = the ancestry value above which a sample is considered to be 
  #                  a "pure" member of a group
  # OUTPUT
  # list with as many elements as K used for ADMIXTURE, each element includes
  #    the row indices for "pure" samples in the respective populations
  #############################3
  pop_ind_list <- list()
  for(i in c(2:ncol(admix_res))){
    tmp_name_num <- i-1
    tmp_name <- paste('g', tmp_name_num, sep = '')
    pop_ind_list[[tmp_name]] <- which(admix_res[[paste('V', i, sep = '')]] >
      group_cutoff)
  }
  return(pop_ind_list)
}

###########
gen_admix_pop_samp_list <- function(admix_res, group_cutoff){
  # Get names of samples with ancestry above "group" cutoff for
  #   each population in admix_res
  # INPUTS
  # admix_res = data.table or data.frame with first column = sample names and
  #              remaining columns from Q matrix outputted by ADMIXTURE
  # group_cutoff = the ancestry value above which a sample is considered to be 
  #                  a "pure" member of a group
  # OUTPUT
  # list with as many elements as K used for ADMIXTURE, each element includes
  #    the sample names of "pure" samples in the respective populations
  #############################3
  pop_ind_list <- list()
  for(i in c(2:ncol(admix_res))){
    tmp_name_num <- i-1
    tmp_name <- paste('g', tmp_name_num, sep = '')
    pop_ind_list[[tmp_name]] <- admix_res$V1[which(
      admix_res[[paste('V', i, sep = '')]] > group_cutoff)]
  }
  return(pop_ind_list)
}


###########

find_share_unique_inds <- function(missing_allele_ind_list, n_snps){
  # Find indices of alleles that are private to individual populations 
  #  or combinations of population
  # INPUTS
  # missing_allele_ind_list = list containing indices of SNPs missing the 
  #   allele of interest in each population; can be generated with 
  #   `lapply(tmp_genos, function(x) which(glMean(x) == 1))` or `== 0`
  # OUTPUT
  # each element contains indices of SNPs where allele is present in a 
  #   population or group of populations but missing from remaining populations
  #################
  allele_share_list <- list()
  for(i in c(1:(length(missing_allele_ind_list)-1))){
  # set the number of pops in the "test" category to look for shared alleles
    n_test <- i
    ncomb_mat <- combn(length(missing_allele_ind_list), m = i)
    for(j in seq(ncol(ncomb_mat))){
      test_inds <- ncomb_mat[, j]
      test_name <- paste(names(missing_allele_ind_list)[test_inds],
        collapse = '_')
      comp_inds <- setdiff(seq(length(missing_allele_ind_list)), test_inds)
      # indices of snps/alleles that are present in all 'test' pops
      test_yes <- setdiff(seq(n_snps), sort(unique(
        unlist(missing_allele_ind_list[test_inds]))))
      # get indices of alleles missing in all 'comp' pops
      comp_miss_inds <- unlist(missing_allele_ind_list[comp_inds])
      comp_tot_miss_inds <- as.numeric(names(which(table(comp_miss_inds) ==
        length(comp_inds))))
      test_priv <- intersect(test_yes, comp_tot_miss_inds)
      allele_share_list[[test_name]] <- test_priv
    }
  }
  return(allele_share_list)
}

#############

gen_share_unique_list <- function(group_ind_list, genos, n_test_reps){
  # Function to run replicate tests identifying unique alleles and alleles
  #  shared by combinations of groups and all groups
  # INPUTS
  # group_ind_list = list containing the sample indices for samples
  #                   that belong to different groups; the indices
  #                   relate to the sample indices in 'genos'
  # genos = genlight genotype object
  # n_test_reps = number of replicates to run
  # OUTPUT
  # list with three elements. Each element has 'n_test_reps' sub-elements
  #    for each test.
  #  'a1' = indices where allele 1 is unique/shared amongst populations
  #  'a2' = indices where allele 2 is unique/shared amonst populations
  #  'all_shared' = indices where minor allele is shared by all populations;
  #                  ie: all populations have both allele 1 and 2
  ######################
  tot_out_list <- list()
  a1_rep_list <- list()
  a2_rep_list <- list()
  all_rep_list <- list()
  min_g_num <- min(unlist(lapply(group_ind_list, length)))
  #
  for(i in seq(n_test_reps)){
    # subsample indices so equal number in each group
    group_sub_inds <- lapply(group_ind_list, sample, size = min_g_num)
    # generate genotype object for each group
    tmp_genos <- lapply(group_sub_inds, function(x)
      genos[indNames(genos) %in% admix_res$V1[x], ])
    # find indices of SNPs missing allele 1 or allele 2 in each group
    a1_missing_inds <- lapply(tmp_genos, function(x) which(glMean(x) == 1))
    a2_missing_inds <- lapply(tmp_genos, function(x) which(glMean(x) == 0))
    # Find private alleles in each population
    a1_shared_alleles <- find_share_unique_inds(a1_missing_inds,
      n_snps = unlist(lapply(tmp_genos, nLoc))[1])
    a1_rep_list[[i]] <- a1_shared_alleles
    #
    a2_shared_alleles <- find_share_unique_inds(a2_missing_inds,
      n_snps = unlist(lapply(tmp_genos, nLoc))[1])
    a2_rep_list[[i]] <- a2_shared_alleles
    #
    # need to find SNPs that are not missing either a1 or a2 in any of the
    #   populations
    all_miss_inds <- unique(c(unlist(a1_missing_inds),
      unlist(a2_missing_inds)))
    n_snps <- unlist(lapply(tmp_genos, nLoc))[1]
    all_shared <- setdiff(seq(n_snps), all_miss_inds)
    all_rep_list[[i]] <- all_shared
  }
  tot_out_list[['a1']] <- a1_rep_list
  tot_out_list[['a2']] <- a2_rep_list
  tot_out_list[['all_shared']] <- all_rep_list
  return(tot_out_list)
}

###########

gen_share_unique_useNames_list <- function(group_name_list, genos, 
  n_test_reps){
  # Function to run replicate tests identifying unique alleles and alleles
  #   shared by combinations of groups and all groups
  #  Adapted from gen_share_unique_list, which uses group indices instead of
  #    names
  # INPUTS
  # group_name_list = list containing the sample names
  #                   that belong to different groups;
  # genos = genlight genotype object
  # n_test_reps = number of replicates to run
  # OUTPUT
  # list with three elements. Each element has 'n_test_reps' sub-elements
  #    for each test.
  #  'a1' = indices where allele 1 is unique/shared amongst populations
  #  'a2' = indices where allele 2 is unique/shared amonst populations
  #  'all_shared' = indices where minor allele is shared by all populations;
  #                  ie: all populations have both allele 1 and 2
  ######################
  tot_out_list <- list()
  a1_rep_list <- list()
  a2_rep_list <- list()
  all_rep_list <- list()
  min_g_num <- min(unlist(lapply(group_name_list, length)))
  #
  for(i in seq(n_test_reps)){
    # subsample indices so equal number in each group
    group_sub_names <- lapply(group_name_list, sample, size = min_g_num)
    # generate genotype object for each group
    tmp_genos <- lapply(group_sub_inds, function(x)
      genos[indNames(genos) %in% x, ])
    # find indices of SNPs missing allele 1 or allele 2 in each group
    a1_missing_inds <- lapply(tmp_genos, function(x) which(glMean(x) == 1))
    a2_missing_inds <- lapply(tmp_genos, function(x) which(glMean(x) == 0))
    # Find private alleles in each population
    a1_shared_alleles <- find_share_unique_inds(a1_missing_inds,
      n_snps = unlist(lapply(tmp_genos, nLoc))[1])
    a1_rep_list[[i]] <- a1_shared_alleles
    #
    a2_shared_alleles <- find_share_unique_inds(a2_missing_inds,
      n_snps = unlist(lapply(tmp_genos, nLoc))[1])
    a2_rep_list[[i]] <- a2_shared_alleles
    #
    # need to find SNPs that are not missing either a1 or a2 in any of the
    #   populations
    all_miss_inds <- unique(c(unlist(a1_missing_inds),
      unlist(a2_missing_inds)))
    n_snps <- unlist(lapply(tmp_genos, nLoc))[1]
    all_shared <- setdiff(seq(n_snps), all_miss_inds)
    all_rep_list[[i]] <- all_shared
  }
  tot_out_list[['a1']] <- a1_rep_list
  tot_out_list[['a2']] <- a2_rep_list
  tot_out_list[['all_shared']] <- all_rep_list
  return(tot_out_list)
}

###########

process_share_results <- function(share_result_list){
  # Calculate the mean values for unique/shared alleles across tests
  # INPUTS
  # share_result_list = unique/shared allele list generated with 
  #                       'gen_share_unique_list' function
  # OUTPUT
  # vector of mean values of unique alleles in individual and groups of
  #   populations as well as SNPs where minor allele is found in all
  #   populations. Names of vector correspond to the pops used.
  ############
  mean_a1_shared <- list()
  mean_a2_shared <- list()
  a1_rep_list <- share_result_list[[1]]
  a2_rep_list <- share_result_list[[2]]
  all_rep_list <- share_result_list[[3]]
  for(gn in names(a1_rep_list[[1]])){
    tmp_a1_val <- mean(unlist(lapply(lapply(a1_rep_list, function(x) x[[gn]]),
      length)))
    mean_a1_shared[[gn]] <- tmp_a1_val
    #
    tmp_a2_val <- mean(unlist(lapply(lapply(a2_rep_list, function(x) x[[gn]]),
      length)))
    mean_a2_shared[[gn]] <- tmp_a2_val
    #
  }
  tot_unique_vec <- unlist(mean_a1_shared) + unlist(mean_a2_shared)
  tot_shared <- mean(unlist(lapply(all_rep_list, length)))
  tot_unique_vec[['all']] <- tot_shared
  return(tot_unique_vec)
}

############

get_shared_alleles_multicutoff <- function(multicut_group_ind_list, genos,
  n_test_reps){
  # Function to calculate the number of private/shared alleles at multiple
  #  "purity" cutoffs
  # INPUTS
  # multicut_group_ind_list = list, each element being a list of indices for
  #    samples belonging to different groups, using different "purity" cutoffs
  #    for each element
  # genos = genlight genotype object
  # n_test_reps = number of replicates to run
  # OUTPUTS
  # list, with each element showing the mean value of private/shared alleles
  #  in each population given the "purity" cutoff
  ###############
  tmp_res_list <- list()
  for(gco in names(multicut_group_ind_list)){
    cutoff_val <- as.numeric(unlist(strsplit(gco, split = '_'))[2])
    tmp_results <- gen_share_unique_list(
      group_ind_list = multicut_group_ind_list[[gco]], genos = genos, 
      n_test_reps = n_test_reps)
    tmp_means <- process_share_results(tmp_results)
    tmp_res_list[[gco]] <- tmp_means
  }
  tmp_df <- data.frame(group = names(tmp_res_list[[1]]), stringsAsFactors = F)
  for(gco in names(tmp_res_list)){
    tmp_col_name = paste('n_alleles_', gco, sep = '')
    tmp_df[, tmp_col_name] <- tmp_res_list[[gco]]
  }
  return(tmp_df)
}






