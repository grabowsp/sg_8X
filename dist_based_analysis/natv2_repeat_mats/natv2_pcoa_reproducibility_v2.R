# Goal:
# For each PCo X:
#  Plot PCo for disomic genotypes
#     PCo1vPCoX (or PCo1 v PCo2)
#     Color by gene pool
#  Mean value boxplots separated by ploidy:
#    abs(tet - dip)
#    abs(poly - tet)
#    abs(poly - ploidy-specific)
#  Mean value boxplotes separated by genepool:
#    abs(poly - ploidy-specific)
#  Identify mean-change outliers
#  Variance/SD boxplots separated by ploidy
#    total variance for diploid genotypes
#    total change-in-SD (see if this plot is better than var)
#    tet - dip
#    poly - dip
#    poly - ploidy-specific
#  Variance boxplots separated by genepool:
#    poly - ploidy-specific
#  Identify variance-change outliers
# PCoA figure with change from diploid to polyploid genotype results
#    shown as line with different start/stop symbols

# delta(sqrt(abs(var1 - var2)))

### LOAD MODULES ###
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(data.table)
library(ggplot2)
library(patchwork)

### INPUT DATA ###

pcoa_list_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2.100k.pcoa_mats.list.rds'
pcoa_list <- readRDS(pcoa_list_file)

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

### SET OUTPUTS ###
out_dir <- '/Users/grabowsk/Documents/Switchgrass_8X/pcoa_reproducibility/'

### SET VARIABLES ###
pure_cut <- 0.9
adm_cut <- 0.05

tet_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '4X')]
oct_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '8X' |
                                        samp_meta$NQUIRE_PLOIDY == '6X')]
oct_names_exact <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '8X')]
hex_names <- samp_meta$VCF_NAME[which(samp_meta$NQUIRE_PLOIDY == '6X')]
### ANALYSIS FUNCTIONS ###
standardize_pco <- function(pco_mat, pco_num, control_name){
  # Standardize a PCo vector from 0 to 1; includes setting a specific sample
  #   as always having a positive value, so that PCo's from replicate runs
  #   (should) have the same sign
  # INPUTS
  # pco_mat = results from cmdscale or some other PCoA function
  # pco_num = the PCo number to analyze
  # control_name = the "control" sample that will always be positive
  # OUTPUT
  # vector of standardized PCo$pco_num values
  ###########
  control_ind <- which(rownames(pco_mat) == control_name)
  tmp_vec <- pco_mat[, pco_num]
  if(tmp_vec[control_ind] < 0){
    tmp_vec <- tmp_vec * -1
  }
  tmp_vec <- tmp_vec - min(tmp_vec)
  tmp_vec <- tmp_vec / max(tmp_vec)
  return(tmp_vec)
}

gen_pco_comp_mat <- function(pco_res_list, pco_num, control_name){
  tmp_pco_list <- list()
  for(dpl in seq(length(pco_res_list))){
    tmp_pco_list[[dpl]] <- standardize_pco(pco_mat = pco_res_list[[dpl]],
                                           pco_num = pco_num, control_name = control_name)
  }
  tmp_pco_mat <- matrix(unlist(tmp_pco_list), nrow = length(tmp_pco_list),
                        byrow = T)
  colnames(tmp_pco_mat) <- names(tmp_pco_list[[1]])
  return(tmp_pco_mat)
}

################
# Asign samples to genepools for plotting
plot_tab <- data.table(samp_name = rownames(pcoa_list[[1]][[1]]), 
                       genepool = as.character(NA), ploidy = as.character(NA))
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

plot_tab[plot_tab$samp_name %in% tet_names, ploidy := '4X']
plot_tab[plot_tab$samp_name %in% hex_names, ploidy := '6X']
plot_tab[plot_tab$samp_name %in% oct_names_exact, ploidy := '8X']
plot_tab[which(is.na(plot_tab$ploidy)), ploidy := '4X']

ploidy_col_vec <- c()
ploidy_col_vec['4X'] <- 'black'
ploidy_col_vec['8X'] <- 'red'
ploidy_col_vec['6X'] <- 'cyan3'

plot_tab[ , ploidy_col := as.character(NA)]
for(pcv in names(ploidy_col_vec)){
  plot_tab[which(plot_tab$ploidy == pcv), ploidy_col := ploidy_col_vec[pcv]]
}

### Add PCoA results to plot_tab
dip_pcoa_list <- pcoa_list[['dip_pcoa_list']]
tet_pcoa_list <- pcoa_list[['tet_pcoa_list']]
poly_pcoa_list <- pcoa_list[['poly_pcoa_list']]

# PCo1
dip_pco_1_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list, 
                                        pco_num = 1, control_name = 'J430.A')
plot_tab[ , dip_pco_1_mean := apply(dip_pco_1_stand_mat, 2, mean)]
plot_tab[ , dip_pco_1_var := apply(dip_pco_1_stand_mat, 2, var)]

tet_pco_1_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list, 
                                        pco_num = 1, control_name = 'J430.A')
plot_tab[ , tet_pco_1_mean := apply(tet_pco_1_stand_mat, 2, mean)]
plot_tab[ , tet_pco_1_var := apply(tet_pco_1_stand_mat, 2, var)]

poly_pco_1_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list, 
                                        pco_num = 1, control_name = 'J430.A')
plot_tab[ , poly_pco_1_mean := apply(poly_pco_1_stand_mat, 2, mean)]
plot_tab[ , poly_pco_1_var := apply(poly_pco_1_stand_mat, 2, var)]

# PCo2
dip_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
                                        pco_num = 2, control_name = 'J660.B')
plot_tab[ , dip_pco_2_mean := apply(dip_pco_2_stand_mat, 2, mean)]
plot_tab[ , dip_pco_2_var := apply(dip_pco_2_stand_mat, 2, var)]

tet_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
                                        pco_num = 2, control_name = 'J660.B')
plot_tab[ , tet_pco_2_mean := apply(tet_pco_2_stand_mat, 2, mean)]
plot_tab[ , tet_pco_2_var := apply(tet_pco_2_stand_mat, 2, var)]

poly_pco_2_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
                                        pco_num = 2, control_name = 'J660.B')
plot_tab[ , poly_pco_2_mean := apply(poly_pco_2_stand_mat, 2, mean)]
plot_tab[ , poly_pco_2_var := apply(poly_pco_2_stand_mat, 2, var)]
# PCo3
dip_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
                                        pco_num = 3, control_name = 'J456.B')
plot_tab[ , dip_pco_3_mean := apply(dip_pco_3_stand_mat, 2, mean)]
plot_tab[ , dip_pco_3_var := apply(dip_pco_3_stand_mat, 2, var)]

tet_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
                                        pco_num = 3, control_name = 'J456.B')
plot_tab[ , tet_pco_3_mean := apply(tet_pco_3_stand_mat, 2, mean)]
plot_tab[ , tet_pco_3_var := apply(tet_pco_3_stand_mat, 2, var)]

poly_pco_3_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
                                         pco_num = 3, control_name = 'J456.B')
plot_tab[ , poly_pco_3_mean := apply(poly_pco_3_stand_mat, 2, mean)]
plot_tab[ , poly_pco_3_var := apply(poly_pco_3_stand_mat, 2, var)]

# PCo4
dip_pco_4_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
                                        pco_num = 4, control_name = 'J645.C')
plot_tab[ , dip_pco_4_mean := apply(dip_pco_4_stand_mat, 2, mean)]
plot_tab[ , dip_pco_4_var := apply(dip_pco_4_stand_mat, 2, var)]

tet_pco_4_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
                                        pco_num = 4, control_name = 'J645.C')
plot_tab[ , tet_pco_4_mean := apply(tet_pco_4_stand_mat, 2, mean)]
plot_tab[ , tet_pco_4_var := apply(tet_pco_4_stand_mat, 2, var)]

poly_pco_4_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
                                        pco_num = 4, control_name = 'J645.C')
plot_tab[ , poly_pco_4_mean := apply(poly_pco_4_stand_mat, 2, mean)]
plot_tab[ , poly_pco_4_var := apply(poly_pco_4_stand_mat, 2, var)]

# PCo5
dip_pco_5_stand_mat <- gen_pco_comp_mat(pco_res_list = dip_pcoa_list,
                                        pco_num = 5, control_name = 'J635.C')
plot_tab[ , dip_pco_5_mean := apply(dip_pco_5_stand_mat, 2, mean)]
plot_tab[ , dip_pco_5_var := apply(dip_pco_5_stand_mat, 2, var)]

tet_pco_5_stand_mat <- gen_pco_comp_mat(pco_res_list = tet_pcoa_list,
                                        pco_num = 5, control_name = 'J635.C')
plot_tab[ , tet_pco_5_mean := apply(tet_pco_5_stand_mat, 2, mean)]
plot_tab[ , tet_pco_5_var := apply(tet_pco_5_stand_mat, 2, var)]

poly_pco_5_stand_mat <- gen_pco_comp_mat(pco_res_list = poly_pcoa_list,
                                        pco_num = 5, control_name = 'J635.C')
plot_tab[ , poly_pco_5_mean := apply(poly_pco_5_stand_mat, 2, mean)]
plot_tab[ , poly_pco_5_var := apply(poly_pco_5_stand_mat, 2, var)]

### Generate PCoA dotplot objects
## PCo1 and PCo2
gg_pco_1_2_dip_gp_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                           y = dip_pco_2_mean, 
                                           color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo2 using disomic genotypes')

gg_pco_1_2_dip_pl_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                              y = dip_pco_2_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo2 using disomic genotypes')
# PCo1 separates MW from Lowland (Gulf + Atlantic)
# PCo2 separates Gulf from non-gulf

# tetrasomic
gg_pco_1_2_tet_gp_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_2_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo2 using tetrasomic genotypes')

gg_pco_1_2_tet_pl_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_2_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo2 using tetrasomic genotypes')

# polyploidy
gg_pco_1_2_poly_gp_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_2_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo2 using polyploid genotypes')

gg_pco_1_2_poly_pl_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_2_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo2 using polyploid genotypes')

pco_1v2_dot_out <- paste(out_dir, 'pco_1v2_threetypes_dotplots.pdf', sep = '')
pdf(pco_1v2_dot_out, height = 10, width = 10)
(gg_pco_1_2_dip_gp_dot + gg_pco_1_2_dip_pl_dot)/
  (gg_pco_1_2_tet_gp_dot + gg_pco_1_2_tet_pl_dot)/
  (gg_pco_1_2_poly_gp_dot + gg_pco_1_2_poly_pl_dot)
dev.off()
# very little noticable change looking at PCo1 and PCo2

## PCo3
# disomic
gg_pco_1_3_dip_gp_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                           y = dip_pco_3_mean, 
                                           color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo3 using disomic genotypes')

gg_pco_1_3_dip_pl_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                              y = dip_pco_3_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo3 using disomic genotypes')
# PCo3 separates Gulf into MX (and TX) at top and FL at bottom

# tetrasomic
gg_pco_1_3_tet_gp_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_3_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo3 using tetrasomic genotypes')

gg_pco_1_3_tet_pl_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_3_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo3 using tetrasomic genotypes')

# polyploid
gg_pco_1_3_poly_gp_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_3_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo3 using polyploid genotypes')

gg_pco_1_3_poly_pl_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_3_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo3 using polyploid genotypes')

pco_1v3_dot_out <- paste(out_dir, 'pco_1v3_threetypes_dotplots.pdf', sep = '')
pdf(pco_1v3_dot_out, height = 10, width = 10)
(gg_pco_1_3_dip_gp_dot + gg_pco_1_3_dip_pl_dot)/
  (gg_pco_1_3_tet_gp_dot + gg_pco_1_3_tet_pl_dot)/
  (gg_pco_1_3_poly_gp_dot + gg_pco_1_3_poly_pl_dot)
dev.off()
# very little noticable change looking at PCo3

## PCo4
# disomic
gg_pco_1_4_dip_gp_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                           y = dip_pco_4_mean, 
                                           color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo4 using disomic genotypes')

gg_pco_1_4_dip_pl_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                              y = dip_pco_4_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo4 using disomic genotypes')
# PCo4 seprates MW by ploidy AND separates N and S Atlantic AND
#   separates E Gulf from W Gulf

# tetrasomic
gg_pco_1_4_tet_gp_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_4_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo4 using tetrasomic genotypes')

gg_pco_1_4_tet_pl_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                                  y = tet_pco_4_mean, 
                                                  color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo4 using tetrasomic genotypes')

# polyploid
gg_pco_1_4_poly_gp_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_4_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo4 using polyploid genotypes')

gg_pco_1_4_poly_pl_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                                  y = poly_pco_4_mean, 
                                                  color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo4 using polyploidy genotypes')

pco_1v4_dot_out <- paste(out_dir, 'pco_1v4_threetypes_dotplots.pdf', sep = '')
pdf(pco_1v4_dot_out, height = 10, width = 10)
(gg_pco_1_4_dip_gp_dot + gg_pco_1_4_dip_pl_dot)/
  (gg_pco_1_4_tet_gp_dot + gg_pco_1_4_tet_pl_dot)/
  (gg_pco_1_4_poly_gp_dot + gg_pco_1_4_poly_pl_dot)
dev.off()
# tetrasomic genotypes remove variation in PCo4 in Midwest genepool;

## PCo5
# disomic
gg_pco_1_5_dip_gp_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                           y = dip_pco_5_mean, 
                                           color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo5 using disomic genotypes')

gg_pco_1_5_dip_pl_dot <- ggplot(plot_tab, aes(x = dip_pco_1_mean, 
                                              y = dip_pco_5_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo5 using disomic genotypes')
# PCo4 seprates MW by ploidy AND separates N and S Atlantic AND
#   separates E Gulf from W Gulf

# tetrasomic
gg_pco_1_5_tet_gp_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_5_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo5 using tetrasomic genotypes')

gg_pco_1_5_tet_pl_dot <- ggplot(plot_tab, aes(x = tet_pco_1_mean, 
                                              y = tet_pco_5_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo5 using tetrasomic genotypes')

# polyploid
gg_pco_1_5_poly_gp_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_5_mean*-1, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo1 vs PCo5 using polyploid genotypes')

gg_pco_1_5_poly_pl_dot <- ggplot(plot_tab, aes(x = poly_pco_1_mean, 
                                              y = poly_pco_5_mean*-1, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo1 vs PCo5 using polyploid genotypes')

pco_1v5_dot_out <- paste(out_dir, 'pco_1v5_threetypes_dotplots.pdf', sep = '')
pdf(pco_1v5_dot_out, height = 10, width = 10)
(gg_pco_1_5_dip_gp_dot + gg_pco_1_5_dip_pl_dot)/
  (gg_pco_1_5_tet_gp_dot + gg_pco_1_5_tet_pl_dot)/
  (gg_pco_1_5_poly_gp_dot + gg_pco_1_5_poly_pl_dot)
dev.off()
# tetrasomic genotypes remove variation in PCo5 in lowland genepools;

## PCo4 v PCo5
# disomic
gg_pco_4_5_dip_gp_dot <- ggplot(plot_tab, aes(x = dip_pco_4_mean, 
                     y = dip_pco_5_mean, 
                     color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo4 vs PCo5 using disomic genotypes')

gg_pco_4_5_dip_pl_dot <- ggplot(plot_tab, aes(x = dip_pco_4_mean, 
                                          y = dip_pco_5_mean, 
                                          color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo4 vs PCo5 using disomic genotypes')

# tetrasomic
gg_pco_4_5_tet_gp_dot <- ggplot(plot_tab, aes(x = tet_pco_4_mean, 
                                              y = tet_pco_5_mean, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo4 vs PCo5 using tetrasomic genotypes')

gg_pco_4_5_tet_pl_dot <- ggplot(plot_tab, aes(x = tet_pco_4_mean, 
                                              y = tet_pco_5_mean, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo4 vs PCo5 using tetrasomic genotypes')

# polyploid
gg_pco_4_5_poly_gp_dot <- ggplot(plot_tab, aes(x = poly_pco_4_mean, 
                                              y = poly_pco_5_mean*-1, 
                                              color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo4 vs PCo5 using polyploid genotypes')

gg_pco_4_5_poly_pl_dot <- ggplot(plot_tab, aes(x = poly_pco_4_mean, 
                                              y = poly_pco_5_mean*-1, 
                                              color = as.factor(ploidy))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$ploidy)), 
                     values = ploidy_col_vec, name = 'Ploidy') +
  labs(title = 'PCo4 vs PCo5 using polyploid genotypes')

pco_4v5_dot_out <- paste(out_dir, 'pco_4v5_threetypes_dotplots.pdf', sep = '')
pdf(pco_4v5_dot_out, height = 10, width = 10)
(gg_pco_4_5_dip_gp_dot + gg_pco_4_5_dip_pl_dot)/
  (gg_pco_4_5_tet_gp_dot + gg_pco_4_5_tet_pl_dot)/
  (gg_pco_4_5_poly_gp_dot + gg_pco_4_5_poly_pl_dot)
dev.off()
# tetrasomic genotypes remove variation in PCo5 in lowland genepools;




test_mat <- as.data.frame(tet_pcoa_list[[1]], stringsAsFactors = F)

ggplot(test_mat, aes(x = V4, y = V5)) + geom_point()




ggplot(plot_tab, aes(x = poly_pco_4_mean, 
                     y = poly_pco_5_mean, 
                     color = as.factor(genepool))) +
  geom_point() +
  scale_color_manual(labels = levels(as.factor(plot_tab$genepool)), 
                     values = gp_col_vec, name = 'Genepool') +
  labs(title = 'PCo4 vs PCo5 using disomic genotypes')
