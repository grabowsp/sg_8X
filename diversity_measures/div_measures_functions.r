# Functions used for calculating diversity measures

calc_gl_het_stats <- function(gl_genos){
  # Calculate heterozygosity and homozygosity stats using genlight object
  #   this helps to take ploidy into account
  #  # Note: calculations take NA's into account and get adjusted appropriately
  # INPUTS
  # gl_genos = genlight genotype object
  # OUTPUT
  # list with 4 elements: 
  #  [['tet_hom']] is percentage of homozygous loci in 4X samples
  #  [['oct_hom']] is percentage of homozygous loci in 8X samples
  #  [['corrected_oct_het']] is the number of heterozygos SNPs standardized
  #  by the pairwise comparisons that would show heterozygosity, ex: 3:1 = 3/6,
  #  2:2 = 4/6
  # [['n_snps']] is the number of SNPs in gl_genos
  ##############
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


