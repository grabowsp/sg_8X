# This analysis was fun on my laptop
install.packages('LDheatmap')

library(LDheatmap)
library(data.table)

all_chr8N_ld_file <- '/Users/grabowsk/Desktop/AvM_Chr08N_all.geno.ld'
all_chr8N_ld <- fread(all_chr8N_ld_file)

natv2_chr8N_ld_file <- '/Users/grabowsk/Desktop/AvM_Chr08N_natv2.geno.ld'
natv2_chr8N_ld <- fread(natv2_chr8N_ld_file)

mwall_chr8N_ld_file <- '/Users/grabowsk/Desktop/AvM_Chr08N_MWall.geno.ld'
mwall_chr8N_ld <- fread(mwall_chr8N_ld_file)

mwintro_chr8N_ld_file <- '/Users/grabowsk/Desktop/AvM_Chr08N_MWintrogressed.geno.ld'
mwintro_chr8N_ld <- fread(mwintro_chr8N_ld_file)

atl_chr8N_ld_file <- '/Users/grabowsk/Desktop/AvM_Chr08N_ATL.geno.ld'
atl_chr8N_ld <- fread(atl_chr8N_ld_file)

gulf_chr8N_ld_file <- '/Users/grabowsk/Desktop/AvM_Chr08N_GULF.geno.ld'
gulf_chr8N_ld <- fread(gulf_chr8N_ld_file)
#

fval_file <- '/Users/grabowsk/Desktop/MW8X_all_Fval_GULF_into_MW_allsnps.txt' 
  fvals <- fread(fval_file)

atl_fval_file <- '/Users/grabowsk/Desktop/MW8X_all_Fval_ATL_into_MW_allsnps.txt'
atl_fvals <- fread(atl_fval_file)

rda_res_file <- '/Users/grabowsk/Desktop/Results_Gulf_Midwest_8x_4x_rda.csv'
rda_res <- fread(rda_res_file)
# need to transfer from Downloads to HA, then load

atl_rda_res_file <- '/Users/grabowsk/Desktop/Results_Atlantic_Midwest_8x_4x_rda.csv'
atl_rda_res <- fread(atl_rda_res_file)


# Make table of Fval results

sd_cut <- 3

##########

## Determine cutoff and rda significant SNPs
# Gulf
outlier_cut <- mean(fvals$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(fvals$F_subgrp_v_pop1, na.rm = T))
outlier_inds <- which(fvals$F_subgrp_v_pop1 > outlier_cut)

fvals[, SNP_NAME := paste(fvals$CHR, fvals$POS, sep = '_')]
rda_overlap_inds <- which(fvals$SNP_NAME %in% rda_res$snp)

fvals[, PLOT_CLASS := fvals$CHR]
fvals[rda_overlap_inds, PLOT_CLASS := 'RDA_SIG']
fvals[intersect(outlier_inds, rda_overlap_inds), PLOT_CLASS := 'OVERLAP']

fvals[, EXP_FVAL := exp(fvals$F_subgrp_v_pop1)]

# Atlantic
atl_outlier_cut <- mean(atl_fvals$F_subgrp_v_pop1, na.rm = T) + (
  sd_cut*sd(atl_fvals$F_subgrp_v_pop1, na.rm = T))
atl_outlier_inds <- which(atl_fvals$F_subgrp_v_pop1 > atl_outlier_cut)

atl_fvals[, SNP_NAME := paste(atl_fvals$CHR, atl_fvals$POS, sep = '_')]
atl_rda_overlap_inds <- which(atl_fvals$SNP_NAME %in% atl_rda_res$snp)

atl_fvals[, PLOT_CLASS := atl_fvals$CHR]
atl_fvals[atl_rda_overlap_inds, PLOT_CLASS := 'RDA_SIG']
atl_fvals[intersect(atl_outlier_inds, atl_rda_overlap_inds),
          PLOT_CLASS := 'OVERLAP']

atl_fvals[, EXP_FVAL := exp(atl_fvals$F_subgrp_v_pop1)]

# get Overlap hit positions

chr08N_hit_pos <- atl_fvals[intersect(which(atl_fvals$CHR == 'Chr08N'), 
                      which(atl_fvals$PLOT_CLASS == 'OVERLAP')), POS]

# all sample LD analysis

all_chr8N_ld_hits <- all_chr8N_ld[which(all_chr8N_ld$POS1 %in% chr08N_hit_pos 
                                        & all_chr8N_ld$POS2 %in% chr08N_hit_pos), ]

all_chr8N_ld_mat <- matrix(data = as.numeric(NA), 
        nrow = length(chr08N_hit_pos), ncol = length(chr08N_hit_pos))
rownames(all_chr8N_ld_mat) <- colnames(all_chr8N_ld_mat) <- chr08N_hit_pos

for(i in seq(length(chr08N_hit_pos)-1)){
  for(j in c((i+1):length(chr08N_hit_pos))){
    all_chr8N_ld_mat[i, j] <- all_chr8N_ld_mat[j, i] <- unlist(
      all_chr8N_ld_hits[which(all_chr8N_ld_hits$POS1 == chr08N_hit_pos[i]
                              & all_chr8N_ld_hits$POS2 == chr08N_hit_pos[j]), 'R^2'])
  }
}

all_chr8N_ld_mat[which(is.na(all_chr8N_ld_mat))] <- 1

allsamps_LDheatmap <- LDheatmap(all_chr8N_ld_mat, 
                                genetic.distances = chr08N_hit_pos)

allsamps_out_file <- '/Users/grabowsk/Desktop/allsamps_Chr08N_overlap_LD_heatmap.pdf'

pdf(allsamps_out_file, width = 5, height = 5)
LDheatmap(all_chr8N_ld_mat, genetic.distances = chr08N_hit_pos)
dev.off()

all_chr8N_pos <- sort(unique(union(all_chr8N_ld$POS1, all_chr8N_ld$POS2)))

sub_chr8N_pos <- sort(unique(c(sample(all_chr8N_pos, size = 200), 
                               chr08N_hit_pos)))

allsnp_chr8N_ld_mat <- matrix(data = as.numeric(NA), 
                           nrow = length(sub_chr8N_pos), 
                           ncol = length(sub_chr8N_pos))
rownames(allsnp_chr8N_ld_mat) <- colnames(allsnp_chr8N_ld_mat) <- sub_chr8N_pos

for(i in seq(length(sub_chr8N_pos)-1)){
    pos1_vec <- which(all_chr8N_ld$POS1 == sub_chr8N_pos[i])
    print(i)
        for(j in c((i+1):length(sub_chr8N_pos))){
      pos2_vec <- which(all_chr8N_ld$POS2 == sub_chr8N_pos[j])
      keep_ind <- intersect(pos1_vec, pos2_vec)
    allsnp_chr8N_ld_mat[i, j] <- allsnp_chr8N_ld_mat[j, i] <- unlist(
      all_chr8N_ld[keep_ind, 'R^2'])
  }
}

allsnp_chr8N_ld_mat[which(is.na(allsnp_chr8N_ld_mat))] <- 1

LDheatmap(allsnp_chr8N_ld_mat, genetic.distances = all_chr8N_pos, add.map = F,
          SNP.name = chr08N_hit_pos)

##########
# natv2 LD analysis

all_chr8N_pos <- sort(unique(union(natv2_chr8N_ld$POS1, natv2_chr8N_ld$POS2)))

sum(chr08N_hit_pos %in% all_chr8N_pos == F) == 0
# TRUE

sub_chr8N_pos <- sort(unique(c(sample(all_chr8N_pos, size = 200), 
                               chr08N_hit_pos)))

natv2_chr8N_ld_mat <- matrix(data = as.numeric(NA), 
                              nrow = length(sub_chr8N_pos), 
                              ncol = length(sub_chr8N_pos))
rownames(natv2_chr8N_ld_mat) <- colnames(natv2_chr8N_ld_mat) <- sub_chr8N_pos

for(i in seq(length(sub_chr8N_pos)-1)){
  pos1_vec <- which(natv2_chr8N_ld$POS1 == sub_chr8N_pos[i])
  print(i)
  for(j in c((i+1):length(sub_chr8N_pos))){
    pos2_vec <- which(natv2_chr8N_ld$POS2 == sub_chr8N_pos[j])
    keep_ind <- intersect(pos1_vec, pos2_vec)
    natv2_chr8N_ld_mat[i, j] <- natv2_chr8N_ld_mat[j, i] <- unlist(
      natv2_chr8N_ld[keep_ind, 'R^2'])
  }
}

natv2_chr8N_ld_mat[which(is.na(natv2_chr8N_ld_mat))] <- 1

LDheatmap(natv2_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F,
          SNP.name = chr08N_hit_pos)

natv2_LD_out <- paste('/Users/grabowsk/Desktop/', 
                      'natv2_Chr08N_LD_heatmap.pdf', sep = '')
pdf(natv2_LD_out, height = 5, width = 5)
LDheatmap(natv2_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F,
          SNP.name = chr08N_hit_pos)
dev.off()

### MW LD analysis
mwall_chr8N_pos <- sort(unique(union(mwall_chr8N_ld$POS1, mw_chr8N_ld$POS2)))

sum(chr08N_hit_pos %in% all_chr8N_pos == F) == 0
# TRUE

sub_chr8N_pos <- sort(unique(c(sample(all_chr8N_pos, size = 200), 
                               chr08N_hit_pos)))

mwall_chr8N_ld_mat <- matrix(data = as.numeric(NA), 
                             nrow = length(sub_chr8N_pos), 
                             ncol = length(sub_chr8N_pos))
rownames(mwall_chr8N_ld_mat) <- colnames(mwall_chr8N_ld_mat) <- sub_chr8N_pos

for(i in seq(length(sub_chr8N_pos)-1)){
  pos1_vec <- which(mwall_chr8N_ld$POS1 == sub_chr8N_pos[i])
  print(i)
  for(j in c((i+1):length(sub_chr8N_pos))){
    pos2_vec <- which(mwall_chr8N_ld$POS2 == sub_chr8N_pos[j])
    keep_ind <- intersect(pos1_vec, pos2_vec)
    mwall_chr8N_ld_mat[i, j] <- mwall_chr8N_ld_mat[j, i] <- unlist(
      mwall_chr8N_ld[keep_ind, 'R^2'])
  }
}

mwall_chr8N_ld_mat[which(is.na(mw_chr8N_ld_mat))] <- 1

LDheatmap(mw_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F,
  SNP.name = chr08N_hit_pos)

mwall_LD_out <- paste('/Users/grabowsk/Desktop/', 
                      'mwall_Chr08N_LD_heatmap.pdf', sep = '')
pdf(mwall_LD_out, height = 5, width = 5)
LDheatmap(mw_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F,
          SNP.name = chr08N_hit_pos)
dev.off()

mwall_hit_inds <- which(rownames(mwall_chr8N_ld_mat) %in% chr08N_hit_pos)

LDheatmap(mwall_chr8N_ld_mat[-mwall_hit_inds, -mwall_hit_inds],
          genetic.distances = sub_chr8N_pos[-mwall_hit_inds], add.map = F)

mwall_noRDA_LD_out <- paste('/Users/grabowsk/Desktop/', 
                            'mwall_Chr08N_noRDA_LD_heatmap.pdf', sep = '')

pdf(mwall_noRDA_LD_out, height = 5, width = 5)
LDheatmap(mwall_chr8N_ld_mat[-mwall_hit_inds, -mwall_hit_inds],
          genetic.distances = sub_chr8N_pos[-mwall_hit_inds], add.map = F)
dev.off()

### Atlantic LD analysis
all_chr8N_pos <- sort(unique(union(atl_chr8N_ld$POS1, atl_chr8N_ld$POS2)))

sum(chr08N_hit_pos %in% all_chr8N_pos == F) == 0
# FALSE

sum(chr08N_hit_pos %in% all_chr8N_pos == F)

sub_chr8N_pos <- sort(sample(all_chr8N_pos, size = 200))
#23 of 27 are missing - I think because they are missing

atl_chr8N_ld_mat <- matrix(data = as.numeric(NA), 
                           nrow = length(sub_chr8N_pos), 
                           ncol = length(sub_chr8N_pos))
rownames(atl_chr8N_ld_mat) <- colnames(atl_chr8N_ld_mat) <- sub_chr8N_pos

for(i in seq(length(sub_chr8N_pos)-1)){
  pos1_vec <- which(atl_chr8N_ld$POS1 == sub_chr8N_pos[i])
  print(i)
  for(j in c((i+1):length(sub_chr8N_pos))){
    pos2_vec <- which(atl_chr8N_ld$POS2 == sub_chr8N_pos[j])
    keep_ind <- intersect(pos1_vec, pos2_vec)
    atl_chr8N_ld_mat[i, j] <- atl_chr8N_ld_mat[j, i] <- unlist(
      atl_chr8N_ld[keep_ind, 'R^2'])
  }
}

atl_chr8N_ld_mat[which(is.na(atl_chr8N_ld_mat))] <- 1

LDheatmap(atl_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F)

atl_LD_out <- paste('/Users/grabowsk/Desktop/', 
                    'atl_Chr08N_LD_heatmap.pdf', sep = '')
pdf(atl_LD_out, height = 5, width = 5)
LDheatmap(atl_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F)
dev.off()

### GULF LD analysis
all_chr8N_pos <- sort(unique(union(gulf_chr8N_ld$POS1, gulf_chr8N_ld$POS2)))

sum(chr08N_hit_pos %in% all_chr8N_pos == F) == 0
# FALSE

sum(chr08N_hit_pos %in% all_chr8N_pos == F)
#6 of 27 are missing - I think because they are fixed

gulf_hit_pos <- chr08N_hit_pos[which(chr08N_hit_pos %in% all_chr8N_pos)]

sub_chr8N_pos <- sort(unique(c(sample(all_chr8N_pos, size = 200), 
                              gulf_hit_pos)))

gulf_chr8N_ld_mat <- matrix(data = as.numeric(NA), 
                           nrow = length(sub_chr8N_pos), 
                           ncol = length(sub_chr8N_pos))
rownames(gulf_chr8N_ld_mat) <- colnames(gulf_chr8N_ld_mat) <- sub_chr8N_pos

for(i in seq(length(sub_chr8N_pos)-1)){
  pos1_vec <- which(gulf_chr8N_ld$POS1 == sub_chr8N_pos[i])
  print(i)
  for(j in c((i+1):length(sub_chr8N_pos))){
    pos2_vec <- which(gulf_chr8N_ld$POS2 == sub_chr8N_pos[j])
    keep_ind <- intersect(pos1_vec, pos2_vec)
    gulf_chr8N_ld_mat[i, j] <- gulf_chr8N_ld_mat[j, i] <- unlist(
      gulf_chr8N_ld[keep_ind, 'R^2'])
  }
}

gulf_chr8N_ld_mat[which(is.na(gulf_chr8N_ld_mat))] <- 1

LDheatmap(gulf_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F,
          SNP.name = gulf_hit_pos)

gulf_LD_out <- paste('/Users/grabowsk/Desktop/', 
                    'gulf_Chr08N_LD_heatmap.pdf', sep = '')
pdf(gulf_LD_out, height = 5, width = 5)
LDheatmap(gulf_chr8N_ld_mat, genetic.distances = sub_chr8N_pos, add.map = F,
          SNP.name = gulf_hit_pos)
dev.off()


