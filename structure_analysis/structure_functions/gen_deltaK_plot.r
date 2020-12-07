# Script to generate deltaK and associated plots based on STRUCTURE results

# INPUTS FOR SCRIPT
# [1] in_file 		= the full path to the file with the Ln_Prob_data 
#                           information
#                          This file needs to be generated using the output 
#			    from all the STRUCTURE runs
# [2] nKs 		= the top K tested; should also be the total number 
#			    of K's tested
# [3] nreps 		= the number of STRUCTURE replicates for each K
# [4] samp_set_name 	= the name of the sample set to be used in the title
#			    of the plots 
#############

# bash
# source activate R_analysis

args <- commandArgs(trailingOnly = TRUE)

rundir_args <- commandArgs(trailingOnly = T)


### LOAD LIBRARIES ###
library(ggplot2)
library(patchwork)

### LOAD INPUTS ###

in_file <- args[1]
#in_file <- '/home/t4c1/WORK/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_tet/expandgeo_alltet_LnProbData.txt'

ln_prob <- read.table(in_file, header = F, sep = '\t', stringsAsFactors = F)

### SET OUTPUT ###
out_fig_suf <- '_deltaK.pdf'
out_fig_file <- gsub('.txt', out_fig_suf, in_file, fixed = T)

out_tab_suf <- '_deltaK.txt'
out_tab_file <- gsub('.txt', out_tab_suf, in_file, fixed = T)

### SET VARIABLES ###
nKs <- as.numeric(args[2])
#nKs <- 10

nreps <- as.numeric(args[3])
#nreps <- 3

samp_set_name <- args[4]
#samp_set_name <- 'expandgeo_alltet'
samp_set_name <- sub('_$', '', samp_set_name)
##################

ln_prob_vals <- unlist(lapply(strsplit(ln_prob[,1], split = ' '),
  function(x) as.numeric(rev(x)[1])))

ln_prob_df <- data.frame(K = rep(seq(nKs), times = 1, each = nreps),
  rep = rep(seq(nreps), times = nKs), ln_prob = ln_prob_vals,
  stringsAsFactors = F)

ln_prob_df$ord_1 <- NA
ln_prob_df$ord_2 <- NA
for(i in seq(nreps)){
  tmp_inds <- which(ln_prob_df$rep == i)
  tmp_vals <- ln_prob_df$ln_prob[tmp_inds]
  tmp_change_1 <- tmp_vals[2:nKs] - tmp_vals[1:(nKs-1)]
  change_1 <- c(NA, tmp_change_1)
  ln_prob_df$ord_1[tmp_inds] <- change_1
  tmp_change_2 <- abs(tmp_change_1[2:(nKs-1)] - tmp_change_1[1:(nKs-2)])
  change_2 <- c(NA, tmp_change_2, NA)
  ln_prob_df$ord_2[tmp_inds] <- change_2
}

mean_ln_prob <- tapply(ln_prob_df$ln_prob, ln_prob_df$K, mean)
sd_ln_prob <- tapply(ln_prob_df$ln_prob, ln_prob_df$K, sd)

mean_change1_K <- tapply(ln_prob_df$ord_1, ln_prob_df$K, mean)
sd_change1_K <- tapply(ln_prob_df$ord_1, ln_prob_df$K, sd)

mean_change2_K <- tapply(ln_prob_df$ord_2, ln_prob_df$K, mean)
sd_change2_K <- tapply(ln_prob_df$ord_2, ln_prob_df$K, sd)

delta_K <- mean_change2_K / sd_ln_prob

delta_df <- data.frame(K = seq(nKs), mean_ln_prob, sd_ln_prob, mean_change1_K,
  sd_change1_K, mean_change2_K, sd_change2_K, delta_K, stringsAsFactors = F)

####

gg_lnprob <- ggplot(delta_df, aes(x = K, y = mean_ln_prob)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_ln_prob - sd_ln_prob,
    ymax = mean_ln_prob + sd_ln_prob)) +
  xlab('K') +
  ylab('ln Prob Data') +
  ggtitle(paste('Ln Prob Data for ', samp_set_name, sep = ''))

gg_change1 <- ggplot(delta_df, aes(x = K, y = mean_change1_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change1_K - sd_change1_K,
    ymax = mean_change1_K + sd_change1_K)) +
  xlab('K') +
  ylab("L'(K)") +
  ggtitle(paste("L'(K) for ", samp_set_name, sep = ''))

gg_change2 <- ggplot(delta_df, aes(x = K, y = mean_change2_K)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_change2_K - sd_change2_K,
    ymax = mean_change2_K + sd_change2_K)) +
  xlab('K') +
  ylab("|L''(K)|") +
  ggtitle(paste("|L''(K)| for ", samp_set_name, sep = ''))

gg_deltaK <- ggplot(delta_df, aes(x = K, y = delta_K)) +
  geom_point() +
  geom_line() +
  xlab('K') +
  ylab("delta K") +
  ggtitle(paste("delta K for ", samp_set_name, sep = ''))

pdf(out_fig_file, width = 12, height = 10)
(gg_lnprob + gg_change1) / (gg_change2 + gg_deltaK)
dev.off()

write.table(ln_prob_df, file = out_tab_file, quote = F, sep= '\t',
  row.names = F, col.names = T)

quit(save = 'no')



