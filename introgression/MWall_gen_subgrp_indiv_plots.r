# Generate plots that are color colded by MW8X subgroups

### LOAD ENVIORNMENT ###
# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)
library(gridExtra)

### INPUT DATA ###
samp_info_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'sg_8X_result_tabs/', 'natv2filt_res_tab_v4.0.txt', sep = '')
samp_info <- fread(samp_info_file)

# introgression analysis results
mw_v_g_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_GULF_into_MW_10SNP_window_results.rds', sep = '')
mw_v_g_list <- readRDS(mw_v_g_file)

mw_v_a_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_ATL_into_MW_10SNP_window_results.rds', sep = '')
mw_v_a_list <- readRDS(mw_v_a_file)

# training sample info
train_samp_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/', 'introgression_samps_and_groups.txt', sep = '')
train_samps <- fread(train_samp_file)

### SET OUTPUTS ###
mw_v_g_adj_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_GULF_into_MW_10SNP_window_ADJ_N_POP2_subgrp_colors.pdf', sep = '')

mw_v_g_score_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_GULF_into_MW_10SNP_window_POP2_SCORE_subgrp_colors.pdf', sep = '')

mw_v_a_adj_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_ATL_into_MW_10SNP_window_ADJ_N_POP2_subgrp_colors.pdf', sep = '')

mw_v_a_score_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MWall_ATL_into_MW_10SNP_window_POP2_SCORE_subgrp_colors.pdf', sep = '')

### Set plotting variables ###
# set y-axis limits for plotting
max_tot_adj_n_pop2 <- 50
#max_tot_adj_n_pop2 <- max(unlist(lapply(tot_pop2_list, 
#  function(x) max(x[[4]]$ADJ_N_POP2))))

max_tot_score <- 12.5
#max_tot_score <- max(unlist(lapply(tot_pop2_list,
#  function(x) max(x[[4]]$POP2_SCORE))))

###############################

# set where to draw Chr lines
cum_chrom_amount <- c()
for(j in unique(mw_v_g_list[[1]][[4]]$CHROM)){
  tmp_amount <- max(mw_v_g_list[[1]][[4]]$POS_CUM[
    which(mw_v_g_list[[1]][[4]]$CHROM == j)])
  cum_chrom_amount <- c(cum_chrom_amount, tmp_amount)
}

##
mw_inds <-  grep('MW_', samp_info$sub_grp)

mw_train_names <- train_samps[train_samps$GROUP == 'Midwest', LIB]

mw_subgrp_names <- sort(unique(samp_info[mw_inds, sub_grp]))

mw_colors <- rainbow(n=length(mw_subgrp_names))
names(mw_colors) <- mw_subgrp_names

mw_g_adj_gglist <- list()
mw_g_score_gglist <- list()
mw_a_adj_gglist <- list()
mw_a_score_gglist <- list()
for(sbg in mw_subgrp_names){
  test_names <- samp_info[which(samp_info$sub_grp == sbg), samp_name]
  # set palette
  test_palette <- rep(c(mw_colors[sbg], 'gray30'), times = 9)
  names(test_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
  test_chr_pal <- scale_colour_manual(name = 'CHROM',  values = test_palette)
  # 
  subgrp_m_g_adj_gglist <- list()
  for(i in seq(length(test_names))){
    tmp_list_ind <- which(names(mw_v_g_list) == test_names[i])
    test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
    train_status <- c()
    if(test_names[i] %in% mw_train_names){train_status <- 'Training Sample; '}
    subgrp_m_g_adj_gglist[[i]] <- ggplot(test_tab,
        aes(x = POS_CUM, y = ADJ_N_POP2)) +
      geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
      test_chr_pal +
      ylim(0,max_tot_adj_n_pop2) +
      geom_vline(xintercept = cum_chrom_amount[c(1:17)],
        color = 'gray50', linetype = 'dashed') +
      xlab('') +
      ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
        '; ', sbg, '; ', train_status, 'ADJ_N_POP2 score; ', 
        'GULF into MW background; 10 SNP windows', sep = ''))
  }
  mw_g_adj_gglist <- c(mw_g_adj_gglist, subgrp_m_g_adj_gglist)
  #
  subgrp_m_g_score_gglist <- list()
  for(i in seq(length(test_names))){
    tmp_list_ind <- which(names(mw_v_g_list) == test_names[i])
    test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
    train_status <- c()
    if(test_names[i] %in% mw_train_names){train_status <- 'Training Sample; '}
    subgrp_m_g_score_gglist[[i]] <- ggplot(test_tab,
        aes(x = POS_CUM, y = POP2_SCORE)) +
      geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
      test_chr_pal +
      ylim(0,max_tot_score) +
      geom_vline(xintercept = cum_chrom_amount[c(1:17)],
        color = 'gray50', linetype = 'dashed') +
      xlab('') +
      ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
        '; ', sbg, '; ', train_status, 'POP2_SCORE; ', 
        'GULF into MW background; 10 SNP windows', sep = ''))
  }
  mw_g_score_gglist <- c(mw_g_score_gglist, subgrp_m_g_score_gglist)
  #
  subgrp_m_a_adj_gglist <- list()
  for(i in seq(length(test_names))){
    tmp_list_ind <- which(names(mw_v_a_list) == test_names[i])
    test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
    train_status <- c()
    if(test_names[i] %in% mw_train_names){train_status <- 'Training Sample; '}
    subgrp_m_a_adj_gglist[[i]] <- ggplot(test_tab,
        aes(x = POS_CUM, y = ADJ_N_POP2)) +
      geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
      test_chr_pal +
      ylim(0,max_tot_adj_n_pop2) +
      geom_vline(xintercept = cum_chrom_amount[c(1:17)],
        color = 'gray50', linetype = 'dashed') +
      xlab('') +
      ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
        '; ', sbg, '; ', train_status, 'ADJ_N_POP2 score; ',
        'ATL into MW background; 10 SNP windows', sep = ''))
  }
  mw_a_adj_gglist <- c(mw_a_adj_gglist, subgrp_m_a_adj_gglist)
  #
  subgrp_m_a_score_gglist <- list()
  for(i in seq(length(test_names))){
    tmp_list_ind <- which(names(mw_v_a_list) == test_names[i])
    test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
    train_status <- c()
    if(test_names[i] %in% mw_train_names){train_status <- 'Training Sample; '}
    subgrp_m_a_score_gglist[[i]] <- ggplot(test_tab,
        aes(x = POS_CUM, y = POP2_SCORE)) +
      geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
      test_chr_pal +
      ylim(0,max_tot_score) +
      geom_vline(xintercept = cum_chrom_amount[c(1:17)],
        color = 'gray50', linetype = 'dashed') +
      xlab('') +
      ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
        '; ', sbg, '; ', train_status, 'POP2_SCORE; ',
        'ATL into MW background; 10 SNP windows', sep = ''))
  }
  mw_a_score_gglist <- c(mw_a_score_gglist, subgrp_m_a_score_gglist)
}

#mw_v_g_adj_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
#  'introgression_v3/MW_analysis/',
#  'MWall_GULF_into_MW_10SNP_window_POP2_SCORE_subgrp_colors.pdf', sep = '')

pdf(mw_v_g_adj_out,  width = 18, height = 2*length(mw_g_adj_gglist))
grid.arrange(grobs = mw_g_adj_gglist, ncol = 1)
dev.off()

pdf(mw_v_g_score_out,  width = 18, height = 2*length(mw_g_score_gglist))
grid.arrange(grobs = mw_g_score_gglist, ncol = 1)
dev.off()

pdf(mw_v_a_adj_out,  width = 18, height = 2*length(mw_a_adj_gglist))
grid.arrange(grobs = mw_a_adj_gglist, ncol = 1)
dev.off()

pdf(mw_v_a_score_out,  width = 18, height = 2*length(mw_a_score_gglist))
grid.arrange(grobs = mw_a_score_gglist, ncol = 1)
dev.off()
   
quit(save = 'no')


