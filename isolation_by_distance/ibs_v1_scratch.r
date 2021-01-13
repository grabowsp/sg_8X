# on HA, to start

library(geosphere)
library(data.table)
library(ecodist)


meta_file <- '/home/f2p1/work/grabowsk/data/switchgrass/sg_8X_metadata_v1.0.csv'
samp_meta <- fread(meta_file)

struc_NI_k3_file <- '/home/f2p1/work/grabowsk/data/switchgrass/GW_50k_geobigNorthInland.3.results.txt'
k3_res <- fread(struc_NI_k3_file, header = F)

ld_sym_dist_res_file <- '/home/f2p1/work/grabowsk/data/switchgrass/plink_files/all_samps/all_samps_Chr01K_ld0.3_symmetric.mdist'
ld_sym_dist_res <- fread(ld_sym_dist_res_file, header = F)
ld_sym_dist_ids <- fread(paste(ld_sym_dist_res_file, 'id', sep = '.'),
  header = F)
ld_sym_dist_mat <- as.matrix(ld_sym_dist_res)
colnames(ld_sym_dist_mat) <- rownames(ld_sym_dist_mat) <- ld_sym_dist_ids$V1

g1_names <- k3_res$V1[k3_res$V2 > 0.9]
g3_names <- k3_res$V1[k3_res$V4 > 0.9]
g2_names <- k3_res$V1[k3_res$V3 > 0.9]
g2_names <- setdiff(g2_names, 'J020.A')


g1_meta_inds <- c()
for(i in g1_names){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  g1_meta_inds <- c(g1_meta_inds, tmp_ind)
}

g1_geo_dist <- distm(x = samp_meta[g1_meta_inds, c('LONGITUDE', 'LATITUDE')])
colnames(g1_geo_dist) <- rownames(g1_geo_dist) <- g1_names

g1_gen_inds <- c()
for(j in g1_names){
  tmp_ind <- which(colnames(ld_sym_dist_mat) == j)
  g1_gen_inds <- c(g1_gen_inds, tmp_ind)
  }

g1_gen_mat <- ld_sym_dist_mat[g1_gen_inds, g1_gen_inds]
g1_mantel <- mantel(as.dist(g1_gen_mat) ~ as.dist(g1_geo_dist), nperm = 10000)
g1_log_mantel <- mantel(as.dist(g1_gen_mat) ~ as.dist(log10(g1_geo_dist + 1)),
  nperm = 10000)

#####

g3_meta_inds <- c()
for(i in g3_names){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  g3_meta_inds <- c(g3_meta_inds, tmp_ind)
}

g3_geo_dist <- distm(x = samp_meta[g3_meta_inds, c('LONGITUDE', 'LATITUDE')])
colnames(g3_geo_dist) <- rownames(g3_geo_dist) <- g3_names

g3_gen_inds <- c()
for(j in g3_names){
  tmp_ind <- which(colnames(ld_sym_dist_mat) == j)
  g3_gen_inds <- c(g3_gen_inds, tmp_ind)
  }

g3_gen_mat <- ld_sym_dist_mat[g3_gen_inds, g3_gen_inds]
g3_mantel <- mantel(as.dist(g3_gen_mat) ~ as.dist(g3_geo_dist), nperm = 10000)
g3_log_mantel <- mantel(as.dist(g3_gen_mat) ~ as.dist(log10(g3_geo_dist + 1)),
  nperm = 10000)

#####

g2_meta_inds <- c()
for(i in g2_names){
  tmp_ind <- which(samp_meta$VCF_NAME == i)
  g2_meta_inds <- c(g2_meta_inds, tmp_ind)
}

g2_geo_dist <- distm(x = samp_meta[g2_meta_inds, c('LONGITUDE', 'LATITUDE')])
colnames(g2_geo_dist) <- rownames(g2_geo_dist) <- g2_names

g2_gen_inds <- c()
for(j in g2_names){
  tmp_ind <- which(colnames(ld_sym_dist_mat) == j)
  g2_gen_inds <- c(g2_gen_inds, tmp_ind)
  }

g2_gen_mat <- ld_sym_dist_mat[g2_gen_inds, g2_gen_inds]
g2_mantel <- mantel(as.dist(g2_gen_mat) ~ as.dist(g2_geo_dist), nperm = 10000)
g2_log_mantel <- mantel(as.dist(g2_gen_mat) ~ as.dist(log10(g2_geo_dist + 1)),
  nperm = 10000)

> g1_mantel
   mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
0.12708323 0.02620000 0.97390000 0.04500000 0.07086086 0.20961011 
> g1_log_mantel
   mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
 0.3586419  0.0001000  1.0000000  0.0001000  0.1938424  0.4598651 
> g2_mantel
    mantelr       pval1       pval2       pval3   llim.2.5%  ulim.97.5% 
-0.00328246  0.48430000  0.51580000  0.97460000 -0.03316773  0.04004052 
> g2_log_mantel
   mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
 0.2473571  0.0001000  1.0000000  0.0001000  0.2050020  0.2917821 
> g3_mantel
   mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
 0.4637020  0.0001000  1.0000000  0.0001000  0.4262223  0.4981558 
> g3_log_mantel
   mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
 0.5661548  0.0001000  1.0000000  0.0001000  0.5383278  0.5933893 



