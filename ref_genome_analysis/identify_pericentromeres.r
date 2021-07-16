# Use gene density to identify pericentromeric regions

## Goal
# use annotation file to calcuate the percentage of sequence in a window size
#  that is genic, then use those metrics to identify pericentromeric regions
# what window size?
# What cutoff?

# bash
# source activate R_analysis

### LOAD PACKAGES ###
library(data.table)

### INPUT DATA ###
gene_bed_file <- '/home/f1p1/tmp/switchgrass_8X/sg_v5_genic.bed'
gene_bed <- fread(gene_bed_file)


### SET VARIABLES ###
window_size <- 1e6

########

calculate_window_gene_density <- function(test_bed, window_size){
  window_start <- seq(min(test_bed[, V2]), 
    max(test_bed[, V3]), by = window_size)
  tmp_per_vec <- c()
  for(i in seq(length(window_start)-1)){
    tmp_wind_inds <- which(test_bed$V2 >= window_start[i] & 
      test_bed$V3 <= window_start[(i+1)])
    tmp_per <- length(unique(unlist(apply(test_bed[tmp_wind_inds], 1,
      function(x) c(x[2]:x[3])))))/window_size
    tmp_per_vec <- c(tmp_per_vec, tmp_per)
  }
  final_window <- c((rev(window_start)[1]+1), max(test_bed$V3))
  fin_wind_inds <- which(test_bed$V2 >= final_window[1] &
    test_bed$V3 <= final_window[2])
  fin_tmp_per <- length(unique(unlist(apply(test_bed[fin_wind_inds], 1,
    function(x) c(x[2]:x[3])))))/(final_window[2]-final_window[1])
  tmp_per_vec <- c(tmp_per_vec, fin_tmp_per)
  names(tmp_per_vec) <- window_start
  return(tmp_per_vec)
}

chr_names <- unique(gene_bed$V1)[grep('Chr', unique(gene_bed$V1))]

gene_density_list <- list()

for(chrom in chr_names){
  tmp_bed <- gene_bed[V1 == chrom,]
  gene_density_list[[chrom]] <- calculate_window_gene_density(
    test_bed = tmp_bed, window_size = window_size)
}

lapply(gene_density_list, function(x) sum(x < 0.05))

#### scratch space below this


chr01K_inds <- which(gene_bed$V1 == 'Chr01K')
chr01K_bed <- gene_bed[chr01K_inds,]

chr01K_density <- calculate_window_gene_density(test_bed = chr01K_bed, 
  window_size = window_size)


chr01K_density <- calculate_window_gene_density(test_bed = chr01K_bed, 
  window_size = 5e5)

window_start <- seq(min(gene_bed[chr01K_inds, V2]), 
  max(gene_bed[chr01K_inds, V3]), by = window_size)

window_1_inds <- which(chr01K_bed$V2 >= window_start[1] & 
  chr01K_bed$V3 <= window_start[2])

test <- length(unique(unlist(apply(chr01K_bed[window_1_inds], 1, 
  function(x) c(x[2]:x[3])))))/window_size

window_2_inds <- which(chr01K_bed$V2 >= window_start[2] & 
  chr01K_bed$V3 <= window_start[3])

test_2 <- length(unique(unlist(apply(chr01K_bed[window_2_inds], 1,
  function(x) c(x[2]:x[3])))))/window_size

tmp_bed <- chr01K_bed
tmp_per_vec <- c()
for(i in seq(length(window_start)-1)){
  tmp_wind_inds <- which(tmp_bed$V2 >= window_start[i] & 
    tmp_bed$V3 <= window_start[(i+1)])
  tmp_per <- length(unique(unlist(apply(tmp_bed[tmp_wind_inds], 1, 
    function(x) c(x[2]:x[3])))))/window_size
  tmp_per_vec <- c(tmp_per_vec, tmp_per)
}
final_window <- c((rev(window_start)[1]+1), max(tmp_bed$V3))
fin_wind_inds <- which(tmp_bed$V2 >= final_window[1] & 
    tmp_bed$V3 <= final_window[2])
fin_tmp_per <- length(unique(unlist(apply(tmp_bed[fin_wind_inds], 1,
    function(x) c(x[2]:x[3])))))/(final_window[2]-final_window[1])
tmp_per_vec <- c(tmp_per_vec, fin_tmp_per)

names(tmp_per_vec) <- window_start

