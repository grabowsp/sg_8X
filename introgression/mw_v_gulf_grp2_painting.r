# Assign genotypes and "paint" chromosomes using windows

### LOAD MODULES ###
# bash
# source activate R_analysis

# args = commandArgs(trailingOnly=T)

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)
library(gridExtra)

chrom_paint_func_file <- paste('/home/grabowsky/tools/workflows/sg_8X/', 
  'introgression/chrom_painting_files/', 'chrom_painting_functions.r', 
  sep = '')
source(chrom_paint_func_file)

### INPUT DATA ###
# vcf_file <- args[1]
vcf_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/', 'grp2_MWvGULF.vcf', sep = '')
vcf_in <- read.table(vcf_file, header = F, stringsAsFactors = F)

vcf_head_tmp <- system(paste('grep CHR ', vcf_file, sep = ''), intern = T)
vcf_head_2 <- sub('#CHROM', 'CHROM', vcf_head_tmp)
vcf_head_3 <- unlist(strsplit(vcf_head_2, split = '\t'))

colnames(vcf_in) <- vcf_head_3
vcf_in <- data.table(vcf_in)

# allele_state_file <- args[2]
allele_state_file <- paste('/home/f2p1/work/grabowsk/data/switchgrass/', 
  'introgression_v3/', 'MW_v_GULF_allele_states.txt', sep = '')
allele_states <- fread(allele_state_file)

### SET OUPUTS ###

window_result_file_out <-  paste('/home/f2p1/work/grabowsk/data/switchgrass/',
  'introgression_v3/', 'grp2_MW_v_GULF_10SNP_window_results.rds', sep = '')

### SET VARIABLES ###
# pop1_name <- args[3]
# pop1 = main background group
pop1_name <- 'MW'

# pop2_name <- args[4]
# pop2 = cadidate introgression group
pop2_name <- 'GULF'

# test_window_size <- args[5]
# test_window_size = the number of SNPs in each window
test_window_size <- 10

###
# Genotype Categories and Color options
## These colors were chosen based on pop1 = MW and pop2 = GULF and should
##   be adjusted accordingly
# cat_a = 1:1 = 3:3 = dark blue (HOM pop1_ONLY)
# cat_b = 2:2 = 4:4 = medium blue (HOM pop1_MAINLY)
# cat_c = 1:Z = 3:Y = dark blue, less saturation (HET pop1_ONLY:NO_INFO)
# cat_d = 2:Z = 4:Y = medium blue, less saturation (HET pop1_MAINLY:NO_INFO)

# cat_e = A:A = C:C = dark red (HOM pop2_ONLY)
# cat_f = B:B = D:D = medium red (HOM pop2_MAINLY)
# cat_g = A:Z = C:Y = dark red, less saturaion (HET pop2_ONLY:NO_INFO)
# cat_h = B:Z = D:Y = medium red, less saturation (HET pop2_MAINLY:NO_INFO)

# cat_i = 1:C = 3:A =  magenta (HET pop1_ONLY:pop2_ONLY)
# cat_j = 1:D = 3:B =  blue-er magenta (more blue) (HET pop1_ONLY:pop2_MAINLY)
# cat_k = 2:C = 4:A =  red-er magenta (more red) (HET pop1_MAINLY:pop2_ONLY)
# cat_l = 2:D = 4:B =  light magenta (HET pop1_MAINLY:pop2_MAINLY)

# cat_m = Z:Z = Y:Y = grey (HOM NO_INFO)
# cat_n = missing data

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
#######

### Introgression weights
pop1_only_weight <- -1
pop1_mainly_weight <- -0.5
pop2_only_weight <- 1
pop2_mainly_weight <- 0.5

intro_wt_vec <- c()
intro_wt_vec$cat_e <- 2*pop2_only_weight
intro_wt_vec$cat_f <- 2*pop2_mainly_weight
intro_wt_vec$cat_g <- pop2_only_weight
intro_wt_vec$cat_h <- pop2_mainly_weight
intro_wt_vec$cat_i <- pop2_only_weight
intro_wt_vec$cat_j <- pop2_mainly_weight
intro_wt_vec$cat_k <- pop2_only_weight
intro_wt_vec$cat_l <- pop2_mainly_weight
intro_wt_vec[c('cat_a', 'cat_b', 'cat_c', 'cat_d', 'cat_m', 'cat_n')] <- 0

#####################
# adjust vcf so SNPs are ordered and there is a cumulative position for
#   plotting
# vector of sample names
samp_vec <- colnames(vcf_in)[c(10:ncol(vcf_in))]
# find max position for each chromosome
max_chrom_pos <- c()
for(chrm in unique(vcf_in$CHROM)){
  max_pos <- max(vcf_in$POS[vcf_in$CHROM == chrm])
  max_chrom_pos[[chrm]] <- max_pos
}
max_chrom_pos <- unlist(max_chrom_pos)
# figure out the cumulative amount to add to each chromosome
cum_chrom_amount <- rep(0, times = length(max_chrom_pos))
for(i in c(2:length(cum_chrom_amount))){
#  cum_chrom_amount[i] <- cum_chrom_amount[(i-1)] + max_chrom_pos[(i-1)]
  cum_chrom_amount[i] <- sum(max_chrom_pos[1:(i-1)])+(100000*(i-1))
}
names(cum_chrom_amount) <- names(max_chrom_pos)

vcf_in[, POS_CUM := as.numeric(NA)]
for(chrm in unique(vcf_in$CHROM)){
  tmp_inds <- which(vcf_in$CHROM == chrm)
  tmp_new_pos <- vcf_in$POS[tmp_inds] + cum_chrom_amount[chrm]
  vcf_in[tmp_inds, POS_CUM := tmp_new_pos]
}

vcf_in <- vcf_in[order(vcf_in$POS_CUM)]

####################
samp_1 <- 'J567.A'
tmp_state_vec <- gen_geno_state_vec(vcf = vcf_in, samp_name = samp_1,
  allele_states = allele_states)
geno_cat_vec <- gen_cat_state_vec(geno_state_vec = tmp_state_vec, 
  pop1_name = 'MW', pop2_name = 'GULF')

geno_score_vec <- rep(NA, times = length(geno_cat_vec))
for(cn in unique(geno_cat_vec)){
  geno_score_vec[which(geno_cat_vec == cn)] <- intro_wt_vec[[cn]]
}

geno_color_vec <- rep(NA, times = nrow(vcf_in))
for(cn in unique(geno_cat_vec)){
  geno_color_vec[which(geno_cat_vec == cn)] <- intro_col_vec[[cn]]
}

#####################

# tot_pop2_list <- readRDS(window_result_file_out)

tot_pop2_list <- list()
for(sv in samp_vec){
  tot_pop2_list[[sv]] <- list()
  # genotype state based on allele states
  tmp_state_vec <- gen_geno_state_vec(vcf = vcf_in, samp_name = sv,
    allele_states = allele_states)
  tot_pop2_list[[sv]][['geno_state_vec']] <- tmp_state_vec
  #
  # genotype categories
  geno_cat_vec <- gen_cat_state_vec(geno_state_vec = tmp_state_vec,
    pop1_name = 'MW', pop2_name = 'GULF')
  tot_pop2_list[[sv]][['geno_cat_vec']] <- geno_cat_vec
  #
  # calculate genotype scores
  geno_score_vec <- rep(NA, times = length(geno_cat_vec))
  for(cn in unique(geno_cat_vec)){
    geno_score_vec[which(geno_cat_vec == cn)] <- intro_wt_vec[[cn]]
  }
  # Generate non-overlapping window table 
  tmp_nonoverlap_ind <- list()
  for(chrm in unique(vcf_in$CHROM)){
    tmp_tab <- gen_pop2_presence_nonoverlap_window_tab(
      samp_state_vec = tmp_state_vec, vcf = vcf_in, chrom = chrm, 
      window_size = test_window_size)
    tmp_score_tab <- gen_pop2_score_nonoverlap_window_tab(
      samp_score_vec = geno_score_vec, vcf = vcf_in, chrom = chrm,
    window_size = test_window_size)  
    tmp_tab[, POP2_SCORE := tmp_score_tab$POP2_SCORE]
    tmp_nonoverlap_ind[[chrm]] <- tmp_tab
  }
  tmp_nonoverlap_1 <- rbindlist(tmp_nonoverlap_ind)
  tot_pop2_list[[sv]][['non_overlap_wind']] <- tmp_nonoverlap_1
  #
  # Generate overlapping window table
  tmp_overlap_ind <- list()
  for(chrm in unique(vcf_in$CHROM)){
    tmp_tab <- gen_pop2_presence_overlap_window_tab(
      samp_state_vec = tmp_state_vec, vcf = vcf_in, chrom = chrm, 
      window_size = test_window_size)
    tmp_score_tab <- gen_pop2_score_overlap_window_tab(
      samp_score_vec = geno_score_vec, vcf = vcf_in, chrom = chrm,
    window_size = test_window_size)  
    tmp_tab[, POP2_SCORE := tmp_score_tab$POP2_SCORE]
    tmp_overlap_ind[[chrm]] <- tmp_tab
  }
  tmp_overlap_1 <- rbindlist(tmp_overlap_ind)
  tot_pop2_list[[sv]][['overlap_wind']] <- tmp_overlap_1
}

saveRDS(tot_pop2_list, file = window_result_file_out)


source(chrom_paint_func_file)

# 3) wrap to do this across many samples


#max_tot_adj_n_pop2 <- max(unlist(lapply(tot_pop2_list, 
#  function(x) max(x[[4]]$ADJ_N_POP2))))
max_tot_adj_n_pop2 <- 50

#max_tot_score <- max(unlist(lapply(tot_pop2_list,
#  function(x) max(x[[4]]$POP2_SCORE))))
max_tot_score <- 12.5

test_tab <- tot_pop2_list[[1]][[4]]

gg_test_1 <- ggplot(test_tab, aes(x = POS_CUM, y = N_ANY_POP2)) +
  geom_line(aes(color = factor(CHROM), group = CHROM)) +
  ylim(0, 10) +
  geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))], 
    color = 'gray50', linetype = 'dashed')

gg_test_2 <- ggplot(test_tab, aes(x = POS_CUM, y = ADJ_N_POP2)) +
  geom_line(aes(color = 'black', group = CHROM)) +
  ylim(0, max_tot_adj_n_pop2) +
  geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))], 
    color = 'gray50', linetype = 'dashed')

gg_test_3 <- ggplot(test_tab, aes(x = POS_CUM, y = POP2_SCORE)) +
  geom_line(aes(color = 'black', group = CHROM)) +
  ylim(0, max_tot_score) +
  geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
    color = 'gray50', linetype = 'dashed')



#test_plot_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/J567.A_Chr01K_introgression_test_plot.pdf'

test_plot_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/J567.A_Chr01K_introgression_test_plot_1.pdf'

pdf(test_plot_file, width = 18, height = 4*3)
gg_test_1 / gg_test_2 / gg_test_3
dev.off()

test_tab_2 <- tot_pop2_list[[2]][[4]]

gg_test_4 <- ggplot(test_tab_2, aes(x = POS_CUM, y = N_ANY_POP2)) +
  geom_line(color = 'blue') +
  ylim(0,10) +
  geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
    color = 'gray50', linetype = 'dashed')

gg_test_5 <- ggplot(test_tab_2, aes(x = POS_CUM, y = ADJ_N_POP2)) +
  geom_line(color = 'blue') +
  ylim(0,max_tot_adj_n_pop2) +
  geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))], 
    color = 'gray50', linetype = 'dashed')

gg_test_6 <- ggplot(test_tab_2, aes(x = POS_CUM, y = POP2_SCORE)) +
  geom_line(color = 'blue') +
  ylim(0, max_tot_score) +
  geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
    color = 'gray50', linetype = 'dashed')

test_plot_file_1 <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/multi_samp_test.pdf'

pdf(test_plot_file_1, width = 18, height = 4*3)
gg_test_1 / gg_test_2 / gg_test_3 / gg_test_4 / gg_test_5 / gg_test_6
dev.off()

test_plot_file_2 <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/multi_samp_test_2.pdf'

gg_list <- list(gg_test_1, gg_test_2, gg_test_3, gg_test_4, gg_test_5, 
  gg_test_6)

pdf(test_plot_file_2, width = 18, height = 4*3)
grid.arrange(grobs = gg_list, ncol = 1)
dev.off()

### Generate ADJ_N_POP2 

black_chr_palette <- rep('black', 
  times = length(unique(tot_pop2_list[[1]][[4]]$CHROM)))
names(black_chr_palette) <- unique(tot_pop2_list[[1]][[4]]$CHROM)
bl_pal <- scale_colour_manual(name = 'CHROM', values = black_chr_palette)

blue_chr_palette <- rep('blue', 
  times = length(unique(tot_pop2_list[[1]][[4]]$CHROM)))
names(blue_chr_palette) <- unique(tot_pop2_list[[1]][[4]]$CHROM)
blue_pal <- scale_colour_manual(name = 'CHROM', values = blue_chr_palette)

tot_adj_gglist <- list()
for(i in seq(length(tot_pop2_list))){
  test_tab <- tot_pop2_list[[i]][[4]]
  tot_adj_gglist[[i]] <- ggplot(test_tab, 
      aes(x = POS_CUM, y = ADJ_N_POP2)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    bl_pal +
    ylim(0,max_tot_adj_n_pop2) +
    geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
      color = 'gray50', linetype = 'dashed')
}

# geom_line(aes(color = factor(CHROM), group = CHROM))

grp2_adj_plot_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/grp2_mw_v_gulf_ADJ_N_GULF_v1.pdf'

pdf(grp2_adj_plot_file, width = 18, height = 2*length(tot_adj_gglist))
grid.arrange(grobs = tot_adj_gglist, ncol = 1)
dev.off()

tot_score_gglist <- list()
for(i in seq(length(tot_pop2_list))){
  test_tab <- tot_pop2_list[[i]][[4]]
  tot_score_gglist[[i]] <- ggplot(test_tab,   
      aes(x = POS_CUM, y = POP2_SCORE)) +
    geom_line(aes(color = CHROM, group = CHROM), show.legend = F) +
    blue_pal + 
    ylim(0,max_tot_score) +
    geom_vline(xintercept = cum_chrom_amount[c(2:length(cum_chrom_amount))],
      color = 'gray50', linetype = 'dashed')
}

grp2_score_plot_file <- '/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/grp2_mw_v_gulf_GULF_SCORE_v1.pdf'

pdf(grp2_score_plot_file, width = 18, height = 2*length(tot_adj_gglist))
grid.arrange(grobs = tot_score_gglist, ncol = 1)
dev.off()


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

gg_col_test <- ggplot(tmp_tab) +
  theme_void() + 
#  geom_point() + 
  xlim(-10000, (max(tmp_tab$POS_CUM)+10000)) + 
  ylim(0, 1) +
  geom_vline(xintercept = tmp_tab$POS_CUM, color = tmp_tab$PLOT_COLOR)


pdf(test_color_plot_1, width = 18, height = 4)
gg_col_test
dev.off()

# Want to add the geom_line(aes(group = CHROM)) addition to make sure
#   the windows from different chromosomes are separated



tmp_tab <- tot_pop2_list[[1]][[4]]
tmp_tab[, PLOT_COLOR := geno_color_vec]


