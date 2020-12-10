# scratch for working out measures of homozygosity

# module load python/3.7-anaconda-2019.10
# source activate adegenet_2_env

library(adegenet)
library(parallel)

div_funct_file <- '/home/grabowsky/tools/workflows/sg_8X/diversity_measures/div_measures_functions.r'
source(div_funct_file)

#geno_file_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.rds'

geno_file_in <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/Chr09N.tetrasomic.CDS.geobig.genlight.rds'

genos <- readRDS(geno_file_in)

tot_divers_list <- lapply(seploc(genos, block.size = 100000, parallel = T),
  calc_gl_het_stats)





tet_inds <- which(ploidy(genos) == 2)
oct_inds <- which(ploidy(genos) == 4)

geno_tab <- apply(as.matrix(genos), 1, table)
n_g_tot <- unlist(lapply(geno_tab, sum))

n_hom_tet <- unlist(lapply(geno_tab[tet_inds], function(x) sum(x['0'],x['2'])))
n_hom_oct <- unlist(lapply(geno_tab[oct_inds], function(x) sum(x['0'],x['4'])))

hom_per_tet <- n_hom_tet / n_g_tot[tet_inds]
hom_per_oct <- n_hom_oct / n_g_tot[oct_inds]

n_1_oct <- unlist(lapply(geno_tab[oct_inds], function(x) sum(x['1'],x['3'])))
n_2_oct <- unlist(lapply(geno_tab[oct_inds], function(x) x['2']))

corrected_het_per_oct <- ((n_2_oct * 2/3)+(n_1_oct * 0.5)) / n_g_tot[oct_inds]


calc_gl_het_stats <- function(gl_genos){
  stat_list <- list()
  tet_inds <- which(ploidy(gl_genos) == 2)
  oct_inds <- which(ploidy(gl_genos) == 4)
  #
  geno_tab <- apply(as.matrix(gl_genos), 1, table)
  n_g_tot <- unlist(lapply(geno_tab, sum))
  n_hom_tet <- unlist(lapply(geno_tab[tet_inds], function(x) 
    sum(x['0'],x['2'])))
  n_hom_oct <- unlist(lapply(geno_tab[oct_inds], function(x) 
    sum(x['0'],x['4'])))
  #
  hom_per_tet <- n_hom_tet / n_g_tot[tet_inds]
  hom_per_oct <- n_hom_oct / n_g_tot[oct_inds]
  #
  n_1_oct <- unlist(lapply(geno_tab[oct_inds], function(x) sum(x['1'],x['3'])))
  n_2_oct <- unlist(lapply(geno_tab[oct_inds], function(x) x['2']))
  corrected_het_per_oct <- ((n_2_oct*2/3)+(n_1_oct*0.5)) / n_g_tot[oct_inds]
  #
  tot_snps <- nLoc(gl_genos)
  #
  stat_list[['tet_hom']] <- hom_per_tet
  stat_list[['oct_hom']] <- hom_per_oct
  stat_list[['corrected_oct_het']] <- corrected_het_per_oct
  stat_list[['n_snps']] <- tot_snps
  return(stat_list)
}

test <- calc_gl_het_stats(genos[, c(1:10000)])

tot_test <- lapply(seploc(genos, block.size = 100000, parallel = T), 
  calc_gl_het_stats)

n_snp_vec <- unlist(lapply(tot_test, function(x) x[['n_snps']]))

norm_fac_vec <- n_snp_vec / max(n_snp_vec)

tet_hom_df <- data.frame(lapply(tot_test, function(x) x[['tet_hom']]))

tet_hom_test <- data.frame(lapply(tot_test, function(x) 
  x[['tet_hom']] * x[['n_snps']]))
tot_tet_hom_test <- apply(tet_hom_test, 1, sum) / nLoc(genos)

tot_tet_hom <- apply(data.frame(lapply(tot_test, function(x) 
  x[['tet_hom']] * x[['n_snps']])), 1, sum) / nLoc(genos)

tot_oct_hom <- apply(data.frame(lapply(tot_test, function(x) 
  x[['oct_hom']] * x[['n_snps']])), 1, sum) / nLoc(genos)

tot_correct_oct_het <- apply(data.frame(lapply(tot_test, function(x)
  x[['corrected_oct_het']] * x[['n_snps']])), 1, sum) / nLoc(genos)



for(i in seq(ncol(tet_hom_df))){
  tet_hom_df[, i] <- tet_hom_df[, i] * norm_fac_vec[i]
}

