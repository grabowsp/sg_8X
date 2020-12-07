# Script for generating PCA figures for 'geobig_northinland' results

#module load python/3.7-anaconda-2019.10
#source activate adegenet_2_env

### LOAD PACKAGES ###
library(adegenet)
library(data.table)
library(parallel)
library(ggplot2)
library(patchwork)

### LOAD INPUTS ###
pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_northinland_pca/GW.50kSNPs.tetrasomic.CDS.geobig_northinland.genlight.PCAresults.rds'
pca_res <- readRDS(pca_res_file)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

### SET OUTPUTS ###
out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_northinland_pca/'
out_short <- 'GW.50kSNPs.geobig_northinland.PCA.pdf'
out_full <- paste(out_dir, out_short, sep = '')

### SET VARIABLES ###
ploidy_color_vec <- c()
ploidy_color_vec['4X'] <- 'red2'
ploidy_color_vec['8X'] <- 'blue2'
ploidy_color_vec['6X'] <- 'gray40'

ploidy_palette <- scale_colour_manual(name = 'Ploidy',
  values = ploidy_color_vec)

############################

pca_eig <- pca_res$eig
per_var_vec <- pca_eig / sum(pca_eig) * 100

pca_mat <- pca_res$scores
pca_df <- data.frame(lib = rownames(pca_mat), pca_mat, stringsAsFactors = F)

meta_ord <- c()
for(j in seq(nrow(pca_df))){
  tmp_ind <- which(samp_meta$VCF_NAME == pca_df$lib[j])
  meta_ord <- c(meta_ord, tmp_ind)
}

pca_df$ploidy <- samp_meta$NQUIRE_PLOIDY[meta_ord]

pca_plot_list <- list()

pcs_to_plot <- c(2:5)

for(i in seq(length(pcs_to_plot))){
  pcX <- 1
  pcY <- pcs_to_plot[i]
  #
  pca_plot_list[[i]] <- ggplot(pca_df, 
    aes_string(x = paste('PC', pcX, sep = ''), 
    y = paste('PC', pcY, sep = ''))) +
    geom_point(aes(color = ploidy)) + ploidy_palette +
    xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
    ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
    ggtitle(paste('PCA for 265 geobig NorthInland samps\nusing 50k SNPs, PC', 
      pcX, ' vs PC', pcY, sep = ''))
}


pdf(out_full, height = 5*2, width = 6*2)
wrap_plots(pca_plot_list, nrow = 2)
dev.off()

quit(save = 'no')

# examine patterns...
samp_meta$STATE[which(samp_meta$VCF_NAME %in% 
  pca_df$lib[which(pca_df$PC1 > -5 & pca_df$PC1 < 2)])]

pca_df$PC3[which(pca_df$PC1 > -5 & pca_df$PC1 < 2)]
pca_df$lib[which(pca_df$PC1 > -5 & pca_df$PC1 < 2)]
samp_meta$STATE[samp_meta$VCF_NAME == 'J400.A']

pca_df$lib[which(pca_df$PC4 < -10)]

samp_meta$STATE[which(samp_meta$VCF_NAME %in% 
  pca_df$lib[which(pca_df$PC1 < -10)])]

table(samp_meta$STATE[which(samp_meta$VCF_NAME %in% 
  pca_df$lib[which(pca_df$PC1 > -10 & pca_df$PC1 < -5)])])

table(samp_meta$STATE[which(samp_meta$VCF_NAME %in%
  pca_df$lib[which(pca_df$PC1 > 2 & pca_df$PC2 < -10)])])

quit(save = 'no')

