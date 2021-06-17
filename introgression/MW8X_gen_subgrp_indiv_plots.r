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
  'MW_8X_GULF_into_MW_10SNP_window_results.rds', sep = '')
mw_v_g_list <- readRDS(mw_v_g_file)

mw_v_a_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MW_8X_ATL_into_MW_10SNP_window_results.rds', sep = '')
mw_v_a_list <- readRDS(mw_v_a_file)

#####

mw8X_si_inds <- intersect(which(samp_info$ploidy == '8X'), 
  grep('MW_', samp_info$sub_grp))

mw01_names <- samp_info[intersect(mw8X_si_inds, 
   which(samp_info$sub_grp == 'MW_01')), samp_name]
mw02_names <- samp_info[intersect(mw8X_si_inds, 
   which(samp_info$sub_grp == 'MW_02')), samp_name]
mw03_names <- samp_info[intersect(mw8X_si_inds,
   which(samp_info$sub_grp == 'MW_03')), samp_name]
mw04_names <- samp_info[intersect(mw8X_si_inds,
   which(samp_info$sub_grp == 'MW_04')), samp_name]
mw05_names <- samp_info[intersect(mw8X_si_inds,
   which(samp_info$sub_grp == 'MW_05')), samp_name]
mw09_names <- samp_info[intersect(mw8X_si_inds,
   which(samp_info$sub_grp == 'MW_09')), samp_name]

###
blue_chrom_palette <- rep(c('blue3', 'gray30'), times = 9)
names(blue_chrom_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
blue_chr_pal <- scale_colour_manual(name = 'CHROM',
  values = blue_chrom_palette)

purple_chrom_palette <- rep(c('purple3', 'gray30'), times = 9)
names(purple_chrom_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
purple_chr_pal <- scale_colour_manual(name = 'CHROM',
  values = purple_chrom_palette)

green_chrom_palette <- rep(c('green4', 'gray30'), times = 9)
names(green_chrom_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
green_chr_pal <- scale_colour_manual(name = 'CHROM',
  values = green_chrom_palette)

yellow_chrom_palette <- rep(c('yellow4', 'gray30'), times = 9)
names(yellow_chrom_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
yellow_chr_pal <- scale_colour_manual(name = 'CHROM',
  values = yellow_chrom_palette)

brown_chrom_palette <- rep(c('brown3', 'gray30'), times = 9)
names(brown_chrom_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
brown_chr_pal <- scale_colour_manual(name = 'CHROM',
  values = yellow_chrom_palette)

palegreen_chrom_palette <- rep(c('palegreen4', 'gray30'), times = 9)
names(palegreen_chrom_palette) <- unique(mw_v_g_list[[1]][[4]]$CHROM)
palegreen_chr_pal <- scale_colour_manual(name = 'CHROM',
  values = palegreen_chrom_palette)

###
# set y-axis limits for plotting
max_tot_adj_n_pop2 <- 50
#max_tot_adj_n_pop2 <- max(unlist(lapply(tot_pop2_list, 
#  function(x) max(x[[4]]$ADJ_N_POP2))))

max_tot_score <- 12.5
#max_tot_score <- max(unlist(lapply(tot_pop2_list,
#  function(x) max(x[[4]]$POP2_SCORE))))

###
cum_chrom_amount <- c()
for(j in unique(mw_v_g_list[[1]][[4]]$CHROM)){
  tmp_amount <- max(mw_v_g_list[[1]][[4]]$POS_CUM[
    which(mw_v_g_list[[1]][[4]]$CHROM == j)])
  cum_chrom_amount <- c(cum_chrom_amount, tmp_amount)
}

######
# MW_01

# MW vs GULF
mw01_mw_v_g_inds <- which(names(mw_v_g_list) %in% mw01_names)

mw01_m_v_g_adj_gglist <- list()
for(i in seq(length(mw01_mw_v_g_inds))){
  tmp_list_ind <- mw01_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw01_m_v_g_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    purple_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], 
    '; MW_01; ', 'ADJ_N_POP2 score; ', 
    'GULF into MW background; 10 SNP windows', sep = ''))
}

mw01_m_v_g_score_gglist <- list()
for(i in seq(length(mw01_mw_v_g_inds))){
  tmp_list_ind <- mw01_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw01_m_v_g_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    purple_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
    '; MW_01; ', 'POP2 score; ',
    'GULF into MW background; 10 SNP windows', sep = ''))
}

# MW vs ATL
mw01_mw_v_a_inds <- which(names(mw_v_a_list) %in% mw01_names)

mw01_m_v_a_adj_gglist <- list()
for(i in seq(length(mw01_mw_v_a_inds))){
  tmp_list_ind <- mw01_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw01_m_v_a_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    purple_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
    '; MW_01; ', 'ADJ_N_POP2 score; ',
    'ATL into MW background; 10 SNP windows', sep = ''))
}

mw01_m_v_a_score_gglist <- list()
for(i in seq(length(mw01_mw_v_a_inds))){
  tmp_list_ind <- mw01_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw01_m_v_a_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    purple_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind],
    '; MW_01; ', 'POP2 score; ',
    'ATL into MW background; 10 SNP windows', sep = ''))
}


######
# MW_02

# MW vs GULF
mw02_mw_v_g_inds <- which(names(mw_v_g_list) %in% mw02_names)

mw02_m_v_g_adj_gglist <- list()
for(i in seq(length(mw02_mw_v_g_inds))){
  tmp_list_ind <- mw02_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw02_m_v_g_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    brown_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_02; ', 
      'ADJ_N_POP2 score; ', 'GULF into MW background; 10 SNP windows', 
      sep = ''))
}

mw02_m_v_g_score_gglist <- list()
for(i in seq(length(mw02_mw_v_g_inds))){
  tmp_list_ind <- mw02_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw02_m_v_g_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    brown_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_02; ',
      'POP2 score; ', 'GULF into MW background; 10 SNP windows',
      sep = ''))
}

# MW vs ATL
mw02_mw_v_a_inds <- which(names(mw_v_a_list) %in% mw02_names)

mw02_m_v_a_adj_gglist <- list()
for(i in seq(length(mw02_mw_v_a_inds))){
  tmp_list_ind <- mw02_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw02_m_v_a_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    brown_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_02; ',
      'ADJ_N_POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

mw02_m_v_a_score_gglist <- list()
for(i in seq(length(mw02_mw_v_a_inds))){
  tmp_list_ind <- mw02_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw02_m_v_a_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    brown_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_02; ',
      'POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

#######
# MW_03

# MW v GULF
mw03_mw_v_g_inds <- which(names(mw_v_g_list) %in% mw03_names)

mw03_m_v_g_adj_gglist <- list()
for(i in seq(length(mw03_mw_v_g_inds))){
  tmp_list_ind <- mw03_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw03_m_v_g_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    green_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_03; ', 
      'ADJ_N_POP2 score; ', 'GULF into MW background; 10 SNP windows', 
      sep = ''))
}

mw03_m_v_g_score_gglist <- list()
for(i in seq(length(mw03_mw_v_g_inds))){
  tmp_list_ind <- mw03_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw03_m_v_g_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    green_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_03; ',
      'POP2 score; ', 'GULF into MW background; 10 SNP windows',
      sep = ''))
}

# MW v ATL
mw03_mw_v_a_inds <- which(names(mw_v_a_list) %in% mw03_names)

mw03_m_v_a_adj_gglist <- list()
for(i in seq(length(mw03_mw_v_a_inds))){
  tmp_list_ind <- mw03_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw03_m_v_a_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    green_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_03; ',
      'ADJ_N_POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

mw03_m_v_a_score_gglist <- list()
for(i in seq(length(mw03_mw_v_a_inds))){
  tmp_list_ind <- mw03_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw03_m_v_a_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    green_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_03; ',
      'POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

#############
# MW_04
# MW v GULF
mw04_mw_v_g_inds <- which(names(mw_v_g_list) %in% mw04_names)

mw04_m_v_g_adj_gglist <- list()
for(i in seq(length(mw04_mw_v_g_inds))){
  tmp_list_ind <- mw04_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw04_m_v_g_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    yellow_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_04; ', 
      'ADJ_N_POP2 score; ', 'GULF into MW background; 10 SNP windows', 
      sep = ''))
}

mw04_m_v_g_score_gglist <- list()
for(i in seq(length(mw04_mw_v_g_inds))){
  tmp_list_ind <- mw04_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw04_m_v_g_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    yellow_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_04; ',
      'POP2 score; ', 'GULF into MW background; 10 SNP windows',
      sep = ''))
}

# MW v ATL
mw04_mw_v_a_inds <- which(names(mw_v_a_list) %in% mw04_names)

mw04_m_v_a_adj_gglist <- list()
for(i in seq(length(mw04_mw_v_a_inds))){
  tmp_list_ind <- mw04_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw04_m_v_a_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    yellow_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_04; ',
      'ADJ_N_POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

mw04_m_v_a_score_gglist <- list()
for(i in seq(length(mw04_mw_v_a_inds))){
  tmp_list_ind <- mw04_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw04_m_v_a_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    yellow_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_04; ',
      'POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

###########
# MW_05

# MW v GULF
mw05_mw_v_g_inds <- which(names(mw_v_g_list) %in% mw05_names)

mw05_m_v_g_adj_gglist <- list()
for(i in seq(length(mw05_mw_v_g_inds))){
  tmp_list_ind <- mw05_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw05_m_v_g_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    palegreen_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_05; ', 
      'ADJ_N_POP2 score; ', 'GULF into MW background; 10 SNP windows', 
      sep = ''))
}

mw05_m_v_g_score_gglist <- list()
for(i in seq(length(mw05_mw_v_g_inds))){
  tmp_list_ind <- mw05_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw05_m_v_g_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    palegreen_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_05; ',
      'POP2 score; ', 'GULF into MW background; 10 SNP windows',
      sep = ''))
}

# MW v ATL
mw05_mw_v_a_inds <- which(names(mw_v_a_list) %in% mw05_names)

mw05_m_v_a_adj_gglist <- list()
for(i in seq(length(mw05_mw_v_a_inds))){
  tmp_list_ind <- mw05_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw05_m_v_a_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    palegreen_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_05; ',
      'ADJ_N_POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

mw05_m_v_a_score_gglist <- list()
for(i in seq(length(mw05_mw_v_a_inds))){
  tmp_list_ind <- mw05_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw05_m_v_a_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    palegreen_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_05; ',
      'POP2 score; ', 'ATL into MW background; 10 SNP windows',
      sep = ''))
}

#############
# MW_09

# MW v GULF
mw09_mw_v_g_inds <- which(names(mw_v_g_list) %in% mw09_names)

mw09_m_v_g_adj_gglist <- list()
for(i in seq(length(mw09_mw_v_g_inds))){
  tmp_list_ind <- mw09_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw09_m_v_g_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    blue_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_09; ', 
      'ADJ_N_POP2 score; ', 'GULF into MW background; 10 SNP windows', 
       sep = ''))
}

mw09_m_v_g_score_gglist <- list()
for(i in seq(length(mw09_mw_v_g_inds))){
  tmp_list_ind <- mw09_mw_v_g_inds[i]
  test_tab <- mw_v_g_list[[tmp_list_ind]][[4]]
  mw09_m_v_g_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    blue_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_09; ',
      'POP2 score; ', 'GULF into MW background; 10 SNP windows',
       sep = ''))
}

# MW v ATL
mw09_mw_v_a_inds <- which(names(mw_v_a_list) %in% mw09_names)

mw09_m_v_a_adj_gglist <- list()
for(i in seq(length(mw09_mw_v_a_inds))){
  tmp_list_ind <- mw09_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw09_m_v_a_adj_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    blue_chr_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_09; ',
      'ADJ_N_POP2 score; ', 'ATL into MW background; 10 SNP windows',
       sep = ''))
}

mw09_m_v_a_score_gglist <- list()
for(i in seq(length(mw09_mw_v_a_inds))){
  tmp_list_ind <- mw09_mw_v_a_inds[i]
  test_tab <- mw_v_a_list[[tmp_list_ind]][[4]]
  mw09_m_v_a_score_gglist[[i]] <- ggplot(test_tab,
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    blue_chr_pal +
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(1:17)],
      color = 'gray50', linetype = 'dashed') +
    xlab('') +
    ggtitle(paste(names(mw_v_g_list)[tmp_list_ind], '; MW_09; ',
      'POP2 score; ', 'ATL into MW background; 10 SNP windows',
       sep = ''))
}


#############
# MW v GULF ADJ_N_POP2 plot
tot_mw_v_g_adj_gglist <- c(mw01_m_v_g_adj_gglist, mw02_m_v_g_adj_gglist, 
  mw03_m_v_g_adj_gglist, mw04_m_v_g_adj_gglist, mw05_m_v_g_adj_gglist, 
  mw09_m_v_g_adj_gglist)

mw_v_g_adj_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MW8X_GULF_into_MW_10SNP_window_ADJ_N_POP2_subgrp_colors.pdf', sep = '')

pdf(mw_v_g_adj_out,  width = 18, height = 2*length(tot_mw_v_g_adj_gglist))
grid.arrange(grobs = tot_mw_v_g_adj_gglist, ncol = 1)
dev.off()

############
# MW v GULF POP2_SCORE plot
tot_mw_v_g_score_gglist <- c(mw01_m_v_g_score_gglist, mw02_m_v_g_score_gglist,
  mw03_m_v_g_score_gglist, mw04_m_v_g_score_gglist, mw05_m_v_g_score_gglist,
  mw09_m_v_g_score_gglist)

mw_v_g_score_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MW8X_GULF_into_MW_10SNP_window_POP2_SCORE_subgrp_colors.pdf', sep = '')

pdf(mw_v_g_score_out,  width = 18, height = 2*length(tot_mw_v_g_score_gglist))
grid.arrange(grobs = tot_mw_v_g_score_gglist, ncol = 1)
dev.off()

###########
# MW v ATL ADJ_N_POP2 plot
tot_mw_v_a_adj_gglist <- c(mw01_m_v_a_adj_gglist, mw02_m_v_a_adj_gglist,
  mw03_m_v_a_adj_gglist, mw04_m_v_a_adj_gglist, mw05_m_v_a_adj_gglist,
  mw09_m_v_a_adj_gglist)

mw_v_a_adj_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MW8X_ATL_into_MW_10SNP_window_ADJ_N_POP2_subgrp_colors.pdf', sep = '')

pdf(mw_v_a_adj_out,  width = 18, height = 2*length(tot_mw_v_a_adj_gglist))
grid.arrange(grobs = tot_mw_v_a_adj_gglist, ncol = 1)
dev.off()

############
# MW v ATL POP2_SCORE plot
tot_mw_v_a_score_gglist <- c(mw01_m_v_a_score_gglist, mw02_m_v_a_score_gglist,
  mw03_m_v_a_score_gglist, mw04_m_v_a_score_gglist, mw05_m_v_a_score_gglist,
  mw09_m_v_a_score_gglist)

mw_v_a_score_out <- paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/MW_analysis/',
  'MW8X_ATL_into_MW_10SNP_window_POP2_SCORE_subgrp_colors.pdf', sep = '')

pdf(mw_v_a_score_out,  width = 18, height = 2*length(tot_mw_v_a_score_gglist))
grid.arrange(grobs = tot_mw_v_a_score_gglist, ncol = 1)
dev.off()



