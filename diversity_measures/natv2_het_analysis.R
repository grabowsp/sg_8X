# Look at Heterozygosity patterns in and natv2

### NEED TO ADJUST TO NATV2 ###

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)

### INPUT DATA ###
samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

natv2_samp_file <- '/Users/grabowsk/Analysis/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_samp_file, header = F)$V1

het_res_file <- '/Users/grabowsk/data/sg_8X_analysis/allsamps.heterozygosity_v1.txt'
allsamp_het_tab <- fread(het_res_file)

allsamp_k3_file <- '/Users/grabowsk/data/sg_8X_analysis/ADMIX_res/GW_200k_allsamps.3.results.txt'
allsamp_k3_res <- fread(allsamp_k3_file)
# $V2 = MW, $V3 = Atlantic, $V4 = Gulf

natv2_k3_file <- '/Users/grabowsk/data/sg_8X_analysis/ADMIX_res/GW_50k_natv2.3.results.txt'
natv2_k3_res <- fread(natv2_k3_file)
# $V2 = Atlantic, $V3 = MW, $V4 = Gulf

### SET VARIABLES ###
pure_cut <- 0.9
adm_cut <- 0.05

### SET OUTPUT ###
out_pdf_1 <- '/Users/grabowsk/Documents/Switchgrass_8X/het/natv2_RawAndAdj_het_violin_v1.pdf'
out_pdf_2 <- '/Users/grabowsk/Documents/Switchgrass_8X/het/natv2_Adj_het_violin_v1.pdf'
####

plot_tab <- allsamp_het_tab[allsamp_het_tab$samp_name %in% natv2_samps, ]
plot_tab[ , genepool := as.character(NA)]

atl_col <- 1
gulf_col <- 3
mw_col <- 2

admix_k3_res <- natv2_k3_res

ancestry_order <- apply(admix_k3_res[, c(2:4)], 1, order, decreasing = T)

mw_names <- admix_k3_res[
  which(admix_k3_res[,c(1+mw_col), with = F] > pure_cut), V1]
plot_tab[plot_tab$samp_name %in% mw_names, genepool := 'Midwest']

atl_names <- admix_k3_res[
  which(admix_k3_res[, c(1+atl_col), with = F] > pure_cut), V1]
plot_tab[plot_tab$samp_name %in% atl_names, genepool := 'Atlantic']

gulf_names <- admix_k3_res[
  which(admix_k3_res[, c(1+gulf_col), with = F] > pure_cut), V1]
plot_tab[plot_tab$samp_name %in% gulf_names, genepool := 'Gulf']

# Samples with ancestry from all 3 gene pools
multi_adm <- admix_k3_res$V1[which(admix_k3_res$V2 > adm_cut & 
                                       admix_k3_res$V3 > adm_cut & 
                                       admix_k3_res$V3 > adm_cut)]
plot_tab[plot_tab$samp_name %in% multi_adm, genepool := 'MW+ATL+GULF']
# admixed samples
mw_atl_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == mw_col & 
                            ancestry_order[2, ] == atl_col)], 
  c(mw_names, multi_adm))
mw_gulf_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == mw_col & 
                            ancestry_order[2, ] == gulf_col)], 
  c(mw_names,multi_adm))

atl_mw_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == atl_col & 
                            ancestry_order[2, ] == mw_col)], 
  c(atl_names, multi_adm))
atl_gulf_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == atl_col & 
                            ancestry_order[2, ] == gulf_col)], 
  c(atl_names, multi_adm))

gulf_mw_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == gulf_col & 
                            ancestry_order[2, ] == mw_col)], 
  c(gulf_names, multi_adm))
gulf_atl_names <- setdiff(
  admix_k3_res$V1[which(ancestry_order[1, ] == gulf_col & 
                            ancestry_order[2, ] == atl_col)], 
  c(gulf_names, multi_adm))

plot_tab[plot_tab$samp_name %in% mw_atl_names, genepool := 'MW+ATL']
plot_tab[plot_tab$samp_name %in% mw_gulf_names, genepool := 'MW+GULF']
plot_tab[plot_tab$samp_name %in% atl_mw_names, genepool := 'ATL+MW']
plot_tab[plot_tab$samp_name %in% atl_gulf_names, genepool := 'ATL+GULF']
plot_tab[plot_tab$samp_name %in% gulf_mw_names, genepool := 'GULF+MW']
plot_tab[plot_tab$samp_name %in% gulf_atl_names, genepool := 'GULF+ATL']

het_v_1 <- ggplot(plot_tab, aes(x = genepool, y = het_raw)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_tab_sub1 <- plot_tab[, -c('het_corrected')]
colnames(plot_tab_sub1)[2] <- 'HET'
plot_tab_sub1$genepool <- paste(plot_tab_sub1$genepool, '_RAW', sep = '')

plot_tab_sub2 <- plot_tab[, -c('het_raw')]
colnames(plot_tab_sub2)[2] <- 'HET'
plot_tab_sub2$genepool <- paste(plot_tab_sub2$genepool, '_ADJ', sep = '')

plot_tab_2 <- rbind(plot_tab_sub1, plot_tab_sub2)

het_v_2 <- ggplot(plot_tab_2, aes(x = genepool, y = HET)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1) +
  ggtitle('Nat.v2 Raw and Adjusted HET') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = out_pdf_1, height = 5, width = 12)
het_v_2
dev.off()

het_v_3 <- ggplot(plot_tab, aes(x = genepool, y = het_corrected)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle('Nat.v2 Adjusted HET')

pdf(file = out_pdf_2, height = 5, width = 8)
het_v_3
dev.off()
