# Chr01K geobig PCA results prelim analysis

# module load python3/3.7-anaconda-2019.10
# source activate r_adegenet_env

### LOAD LIBRARIES ###
library(adegenet)
library(parallel)
library(ggplot2)
library(patchwork)
#library(data.table)

### LOAD DATA ###

pca_res_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_pca/Chr01K.tet.CDS.geobig.prelim.PCA.rds'
pca_res <- readRDS(pca_res_in)

meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v1.0.csv'
samp_meta <- read.table(meta_file, header = T, stringsAsFactors = F, sep = ',')

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

ploidy_color_vec <- c()
ploidy_color_vec['4X'] <- 'red2'
ploidy_color_vec['8X'] <- 'blue2'
ploidy_color_vec['6X'] <- 'gray40'

ploidy_palette <- scale_colour_manual(name = 'Ploidy', 
  values = ploidy_color_vec)


pca_plot_list <- list()

pcs_to_plot <- c(2:5)

for(i in seq(length(pcs_to_plot))){
  pcX <- 1
  pcY <- pcs_to_plot[i]

pca_plot_list[[i]] <- ggplot(pca_df, aes_string(x = paste('PC', pcX, sep = ''), 
  y = paste('PC', pcY, sep = ''))) +
  geom_point(aes(color = ploidy)) + ploidy_palette +
  xlab(paste('PC', pcX, ' (', round(per_var_vec[pcX], 2), '%)', sep = '')) +
  ylab(paste('PC', pcY, ' (', round(per_var_vec[pcY], 2), '%)', sep = '')) +
  ggtitle(paste('Prelim PCA for geobig samps and Chr01K\nPC', pcX, ' vs PC',
    pcY, sep = ''))
}

plot_out <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_pca/Chr01K_geobig_prelim_pca_plots.pdf'

pdf(plot_out, height = 5 * 2, width = 6 * 2)
wrap_plots(pca_plot_list, nrow = 2)
dev.off()

quit(save = 'no')

