### LOAD PACKAGES ###
library(phangorn)
library(data.table)
### INPUT DATA ###

#nex.net.file <- '/Users/grabowsk/Desktop/GW.natv2.ancestral.sub.min4.NNet_Ignore.nexus.nex'
nex.net.file <- '/Users/grabowsk/Desktop/GW.natv2.ancestral.sub.min4.NNet_Ignore.RootedEqAngle.nexus.nex'
nnet_info <- read.nexus.networx(nex.net.file)

samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

natv2_samp_file <- '/Users/grabowsk/Analysis/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt'
natv2_samps <- fread(natv2_samp_file, header = F)$V1

allsamp_k3_file <- '/Users/grabowsk/data/sg_8X_analysis/ADMIX_res/GW_200k_allsamps.3.results.txt'
allsamp_k3_res <- fread(allsamp_k3_file)
# $V2 = MW, $V3 = Atlantic, $V4 = Gulf

natv2_k3_file <- '/Users/grabowsk/data/sg_8X_analysis/ADMIX_res/GW_50k_natv2.3.results.txt'
natv2_k3_res <- fread(natv2_k3_file)
# $V2 = Atlantic, $V3 = MW, $V4 = Gulf

### SET VARIABLES ###
pure_cut <- 0.9
adm_cut <- 0.05
  
#########

# divide samples into groups

plot_tab <- data.table(samp_name = nnet_info[['tip.label']], 
                       genepool = as.character(NA) )

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

plot_tab[which(is.na(plot_tab$genepool)), genepool := 'NONE']

# assign colors to groups
gp_col_vec <- c()
gp_col_vec['Atlantic'] <- 'yellow2'
gp_col_vec['Gulf'] <- 'red2'
gp_col_vec['Midwest'] <- 'blue2'
gp_col_vec['ATL+GULF'] <-'goldenrod3'
gp_col_vec['ATL+MW'] <- 'green2'
gp_col_vec['GULF+ATL'] <- 'orangered1'
gp_col_vec['GULF+MW'] <- 'orchid2'
gp_col_vec['MW+ATL'] <- 'forestgreen'
gp_col_vec['MW+GULF'] <- 'purple3'
gp_col_vec['MW+ATL+GULF'] <- 'sienna4'
gp_col_vec['NONE'] <- 'black'

plot_tab[ , gp_col := as.character(NA)]
for(gpname in names(gp_col_vec)){
  plot_tab[which(plot_tab$genepool == gpname), gp_col := gp_col_vec[gpname]]
}

plot_tab[, ploidy := as.character(NA)]
tet_names <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '4X'],
                       plot_tab$samp_name)
hex_names <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '6X'],
                       plot_tab$samp_name)
oct_names <- intersect(samp_meta$VCF_NAME[samp_meta$NQUIRE_PLOIDY == '8X'],
                       plot_tab$samp_name)

plot_tab[plot_tab$samp_name %in% tet_names, ploidy := '4X']
plot_tab[plot_tab$samp_name %in% hex_names, ploidy := '6X']
plot_tab[plot_tab$samp_name %in% oct_names, ploidy := '8X']
plot_tab[which(is.na(plot_tab$ploidy)), ploidy := '4X']

ploidy_col_vec <- c()
ploidy_col_vec['4X'] <- 'black'
ploidy_col_vec['8X'] <- 'red'
ploidy_col_vec['6X'] <- 'cyan3'

plot_tab[ , ploidy_col := as.character(NA)]
for(pcv in names(ploidy_col_vec)){
  plot_tab[which(plot_tab$ploidy == pcv), ploidy_col := ploidy_col_vec[pcv]]
}

# plot network
test_pdf_file <- '/Users/grabowsk/Desktop/test_NNet.pdf'
size_mult <- 15
pdf(file = test_pdf_file, width = 10*size_mult, height = 10*size_mult)
plot(nnet_info, '2D', show.tip.label = T, use.edge.length = T, edge.width = 1,
     tip.color = plot_tab$gp_col)
dev.off()

gp_pdf_file <- '/Users/grabowsk/Desktop/natv2_NNet_rooted_genepools.pdf'
size_mult <- 5
pdf(file = gp_pdf_file, width = 10*size_mult, height = 10*size_mult)
plot(nnet_info, '2D', show.tip.label = T, use.edge.length = T, edge.width = 1,
     tip.color = plot_tab$gp_col)
dev.off()

test_pdf_2 <- '/Users/grabowsk/Desktop/test_NNet_ploidy.pdf'
ploidy_pdf_file <- '/Users/grabowsk/Desktop/natv2_NNet_rooted_ploidy.pdf'
pdf(file = ploidy_pdf_file, width = 10*size_mult, height = 10*size_mult)
plot(nnet_info, '2D', show.tip.label = T, use.edge.length = T, edge.width = 1,
     tip.color = plot_tab$ploidy_col)
dev.off()





