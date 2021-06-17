# Functions for analysis of population-level introgression analysis based on
#   comparisons of training populations; originally inspired by the idea of 
#   chromosome painting

add_cumulative_pos <- function(ref_freq_tab, chrom_space = 1e5){
  #######
  # calculate the cumulative position of SNPs across all chromosome in 
  #  'ref_freq_tab' so can plot linearly
  ### INPUTS ###
  # ref_freq_tab = table with CHR, POS, and REF_freq columns for different
  #                   populations
  # chrom_space = the space to put between the last postion of a chromosome
  #                 and  position 1 of the next chromosome
  ### OUTPUT ###
  # same date table as 'ref_freq_tab' but with an additional column, POS_CUM
  #########
  max_chrom_pos <- c()
  for(chrm in unique(ref_freq_tab$CHR)){
    max_pos <- max(ref_freq_tab$POS[ref_freq_tab$CHR == chrm])
    max_chrom_pos[[chrm]] <- max_pos
  }
  max_chrom_pos <- unlist(max_chrom_pos)
  # figure out the cumulative amount to add to each chromosome
  cum_chrom_amount <- rep(0, times = length(max_chrom_pos))
  for(i in c(2:length(cum_chrom_amount))){
    cum_chrom_amount[i] <- sum(max_chrom_pos[1:(i-1)])+(chrom_space*(i-1))
  }
  names(cum_chrom_amount) <- names(max_chrom_pos)
  #
  ref_freq_tab[, POS_CUM := as.numeric(NA)]
  for(chrm in unique(ref_freq_tab$CHR)){
    tmp_inds <- which(ref_freq_tab$CHR == chrm)
    tmp_new_pos <- ref_freq_tab$POS[tmp_inds] + cum_chrom_amount[chrm]
    ref_freq_tab[tmp_inds, POS_CUM := tmp_new_pos]
  }
  return(ref_freq_tab)
}

gen_Fval_nonoverlap_window_tab <- function(ref_freq_tab, window_size){
  ###########
  # within NON-overlapping windows, calculate the cumulative 
  #   F_sub_v_pop1 values for all the SNPs in the window
  #  Note: F_sub_v_pop1 values are (subgroup REF freq - pop1 REF freq)
  #   corrected so that positive = in the direction of pop2 REF freq; and 
  #   divided by the range between pop1 REF freq and pop2 REF freq 
  ### INPUTS ###
  # ref_freq_tab = table that contains 'F_sub_v_pop1' column with the values
  #                  to sum up in the table
  # window_size = the number of SNPs in each window
  ### OUTPUT ###
  # data.table with
  #  CHROM = chromosome with window, 
  #  POS = position of beginning of each window,
  #  POS_CUM = cumulative position of beginning of each window, for plotting
  #  F_SUB_V_POP1 = sum of F_subgrp_v_pop1 in each window
  ################
  tot_list <- list()
  for(chrom in unique(ref_freq_tab$CHR)){
    chr_inds <- grep(chrom, ref_freq_tab$CHR)
    window_start_inds <- chr_inds[seq(1, length(chr_inds), by = window_size)]
    window_start_inds <- window_start_inds[-length(window_start_inds)]
    wind_nonov_Fval_vec <- c()
    for(i in window_start_inds){
      test_inds <- c(i:(i+(window_size - 1)))
      tmp_Fval_sum <- sum(ref_freq_tab$F_subgrp_v_pop1[test_inds], na.rm = T)
      wind_nonov_Fval_vec <- c(wind_nonov_Fval_vec, tmp_Fval_sum)    
    }
    nonoverlap_plot_tab <- data.table(
      CHROM = ref_freq_tab$CHR[window_start_inds],
      POS = ref_freq_tab$POS[window_start_inds], 
      POS_CUM = ref_freq_tab$POS_CUM[window_start_inds],
      F_SUB_V_POP1 = wind_nonov_Fval_vec)
    tot_list[[chrom]] <- nonoverlap_plot_tab
  }
  tot_tab <- rbindlist(tot_list)
  return(tot_tab)
}

gen_Fval_overlap_window_tab <- function(ref_freq_tab, window_size){
  ###########
  # within OVERLAPPING windows, calculate the cumulative 
  #   F_sub_v_pop1 values for all the SNPs in the window
  #  Note: F_sub_v_pop1 values are (subgroup REF freq - pop1 REF freq)
  #   corrected so that positive = in the direction of pop2 REF freq; and 
  #   divided by the range between pop1 REF freq and pop2 REF freq 
  ### INPUTS ###
  # ref_freq_tab = table that contains 'F_sub_v_pop1' column with the values
  #                  to sum up in the table
  # window_size = the number of SNPs in each window
  ### OUTPUT ###
  # data.table with
  #  CHROM = chromosome with window, 
  #  POS = position of beginning of each window,
  #  POS_CUM = cumulative position of beginning of each window, for plotting
  #  F_SUB_V_POP1 = sum of F_subgrp_v_pop1 in each window
  ################
  tot_list <- list()
  for(chrom in unique(ref_freq_tab$CHR)){
    chr_inds <- grep(chrom, ref_freq_tab$CHR)
    window_inds <- chr_inds[1:(length(chr_inds)-(window_size - 1))]
    wind_ov_Fval_vec <- c()
    for(i in window_inds){
      test_inds <- c(i:(i+(window_size - 1)))
      tmp_Fval_sum <- sum(ref_freq_tab$F_subgrp_v_pop1[test_inds], na.rm = T)
      wind_ov_Fval_vec <- c(wind_ov_Fval_vec, tmp_Fval_sum)    
    }
    overlap_plot_tab <- data.table(
      CHROM = ref_freq_tab$CHR[window_inds],
      POS = ref_freq_tab$POS[window_inds], 
      POS_CUM = ref_freq_tab$POS_CUM[window_inds],
      F_SUB_V_POP1 = wind_ov_Fval_vec)
    tot_list[[chrom]] <- overlap_plot_tab
  }
  tot_tab <- rbindlist(tot_list)
  return(tot_tab)
}







