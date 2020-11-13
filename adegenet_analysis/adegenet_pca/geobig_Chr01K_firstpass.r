# First pass on PCA using Chr01K SNPs for `geobig` sampleset

# module load python3/3.7-anaconda-2019.10
# source activate r_adegenet_env

library(adegenet)
library(parallel)

gl_file_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/Chr01K.tetrasomic.CDS.geobig.genlight.rds'

tot_gl <- readRDS(gl_file_in)

sub_set_inds <- sort(sample(seq(nLoc(tot_gl)), size = 5e4))

sub_gl <- tot_gl[, sub_set_inds]

n_eig <- nInd(sub_gl) -1

sub_pca <- glPca(sub_gl, nf = n_eig, loadings = F, alleleAsUnit = F, useC = F)

out_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/geobig_pca/'

out_file <- paste(out_dir, 'Chr01K.tet.CDS.geobig.prelim.PCA.rds', sep = '')

saveRDS(sub_pca, out_file)

quit(save = 'no')

