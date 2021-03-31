# Generate results table for natv2filt population structure analysis

#module load python/3.7-anaconda-2019.10
#source activate adegenet_2_env

### LOAD PACKAGES ###
library(adegenet)
library(data.table)

### INPUT DATA ###
samp_meta_file <- '/global/homes/g/grabowsp/data/switchgrass/metadata_8X/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

pca_res_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/pca_analysis/natv2filt_pca/GW.100kSNPs.tetrasomic.CDS.natv2filt.genlight.PCAresults.rds'
pca_res <- readRDS(pca_res_file)

admix_k2_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/natv2filt_admix/GW_100k_natv2filt.2.results.txt'
admix_k2 <- fread(admix_k2_file)

admix_k3_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/natv2filt_admix/GW_100k_natv2filt.3.results.txt'
admix_k3 <- fread(admix_k3_file)

### SET OUTPUT ###

out_file <- '/global/homes/g/grabowsp/data/switchgrass/results_tables_8X/natv2filt_res_tab_v1.0.txt'


### FUNCTIONS ###
standardize_pca_res <- function(pca_res_vec){
  tmp_vec <- pca_res_vec - min(pca_res_vec)
  tmp_vec_2 <- tmp_vec / max(tmp_vec)
  return(tmp_vec_2)
}

########

meta_order <- c()
for(j in admix_k2$V1){
  tmp_ind <- which(samp_meta$VCF_NAME == j)
  meta_order <- c(meta_order, tmp_ind)
}

res_tab <- data.table(samp_name = admix_k2$V1, 
             ploidy = samp_meta$NQUIRE_PLOIDY[meta_order], 
             full_admix_k2_MW = admix_k2$V2,
             full_admix_k2_LOW = admix_k2$V3, full_admix_k3_ATL = admix_k3$V2, 
             full_admix_k3_MW = admix_k3$V3, full_admix_k3_GULF = admix_k3$V4)

pca_res_order <- c()
for(i in res_tab$samp_name){
  tmp_ind <- which(rownames(pca_res$scores) == i)
  pca_res_order <- c(pca_res_order, tmp_ind)
}

res_tab[ , full_pc01_raw := pca_res$scores[pca_res_order, 1]]
res_tab[ , full_pc01_stand := 
  standardize_pca_res(pca_res$scores[,1])[pca_res_order]]

res_tab[ , full_pc02_raw := pca_res$scores[pca_res_order, 2]]
res_tab[ , full_pc02_stand := 
  standardize_pca_res(pca_res$scores[,2])[pca_res_order]]

res_tab[ , full_pc03_raw := pca_res$scores[pca_res_order, 3]]
res_tab[ , full_pc03_stand :=
  standardize_pca_res(pca_res$scores[,3])[pca_res_order]]

res_tab[ , full_pc04_raw := pca_res$scores[pca_res_order, 4]]
res_tab[ , full_pc04_stand :=
  standardize_pca_res(pca_res$scores[,4])[pca_res_order]]

res_tab[ , full_pc05_raw := pca_res$scores[pca_res_order, 5]]
res_tab[ , full_pc05_stand :=
  standardize_pca_res(pca_res$scores[,5])[pca_res_order]]

fwrite(res_tab, file = out_file, sep = '\t')

# explore cutoffs

summary(res_tab$full_pc01_raw[res_tab$full_admix_k3_MW > 0.9])




