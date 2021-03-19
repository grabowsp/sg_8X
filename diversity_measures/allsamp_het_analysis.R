# Look at Heterozygosity patterns in allsamps and natv2

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

### SET VARIABLES ###
pure_cut <- 0.9
adm_cut <- 0.05

### SET OUTPUT ###
out_pdf <- '/Users/grabowsk/Documents/Switchgrass_8X/het/allsamp_het_violin_v1.pdf'

####

plot_tab <- allsamp_het_tab
plot_tab$samp_name <- gsub('^X', '', plot_tab$samp_name)
plot_tab[ , genepool := as.character(NA)]

allsamp_k3_res$V1 <- gsub('-', '.', allsamp_k3_res$V1)

ancestry_order <- apply(allsamp_k3_res[, c(2:4)], 1, order, decreasing = T)

mw_names <- allsamp_k3_res[V2 > pure_cut, V1]
plot_tab[plot_tab$samp_name %in% mw_names, genepool := 'Midwest']

atl_names <- allsamp_k3_res[V3 > pure_cut, V1]
plot_tab[plot_tab$samp_name %in% atl_names, genepool := 'Atlantic']

gulf_names <- allsamp_k3_res[V4 > pure_cut, V1]
plot_tab[plot_tab$samp_name %in% gulf_names, genepool := 'Gulf']

# Samples with ancestry from all 3 gene pools
multi_adm <- allsamp_k3_res$V1[which(allsamp_k3_res$V2 > adm_cut & 
                                       allsamp_k3_res$V3 > adm_cut & 
                                       allsamp_k3_res$V3 > adm_cut)]
plot_tab[plot_tab$samp_name %in% multi_adm, genepool := 'MW+ATL+GULF']
# admixed samples
mw_atl_names <- setdiff(
  allsamp_k3_res$V1[which(ancestry_order[1, ] == 1 & 
                            ancestry_order[2, ] == 2)], 
  c(mw_names, multi_adm))
mw_gulf_names <- setdiff(
  allsamp_k3_res$V1[which(ancestry_order[1, ] == 1 & 
                            ancestry_order[2, ] == 3)], 
  c(mw_names,multi_adm))

atl_mw_names <- setdiff(
  allsamp_k3_res$V1[which(ancestry_order[1, ] == 2 & 
                            ancestry_order[2, ] == 1)], 
  c(atl_names, multi_adm))
atl_gulf_names <- setdiff(
  allsamp_k3_res$V1[which(ancestry_order[1, ] == 2 & 
                            ancestry_order[2, ] == 3)], 
  c(atl_names, multi_adm))

gulf_mw_names <- setdiff(
  allsamp_k3_res$V1[which(ancestry_order[1, ] == 3 & 
                            ancestry_order[2, ] == 1)], 
  c(gulf_names, multi_adm))
gulf_atl_names <- setdiff(
  allsamp_k3_res$V1[which(ancestry_order[1, ] == 3 & 
                            ancestry_order[2, ] == 2)], 
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
  ggtitle('Allsamps Het') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pdf(file = out_pdf, height = 5, width = 10)
het_v_2
dev.off()
