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
library(ggplot2)
library(patchwork)

### INPUT DATA ###
gene_bed_file <- '/home/f1p1/tmp/switchgrass_8X/sg_v5_genic.bed'
gene_bed <- fread(gene_bed_file)

### SET OUTPUT ###

peri_out <- '/home/f1p1/tmp/switchgrass_8X/sg_v5_pericentromere_pos_.075density.txt'

density_out <- '/home/f1p1/tmp/switchgrass_8X/sg_v5_gene_density_5e5_window.rds'

### SET VARIABLES ###
# window_size <- 1e6
# window sizes to test 5e5, 1e6, 5e6, 1e7

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

# 5e5 window size
gene_density_5e5 <- list()
for(chrom in chr_names){
  tmp_bed <- gene_bed[V1 == chrom,]
  gene_density_5e5[[chrom]] <- calculate_window_gene_density(
    test_bed = tmp_bed, window_size = 5e5)
}

# Decision: 
# 5e5 bp window
# 0.075 density cutoff
# end with 2+ spaces between windows below cutoff
# add "+1" to the end index of the pericentromere because that is essentially
#  the end chromosomal position of the last chosen window

peri_list <- list()

which(gene_density_5e5[[1]] < 0.075)
peri_list[[1]] <- as.numeric(names(gene_density_5e5[[1]])[c(46,(52+1))])

which(gene_density_5e5[[2]] < 0.075)
peri_list[[2]] <- as.numeric(names(gene_density_5e5[[2]])[c(43,(57+1))])

which(gene_density_5e5[[3]] < 0.075)
peri_list[[3]] <- as.numeric(names(gene_density_5e5[[3]])[c(50,(65+1))])

which(gene_density_5e5[[4]] < 0.075)
peri_list[[4]] <- as.numeric(names(gene_density_5e5[[4]])[c(56,(67+1))])

which(gene_density_5e5[[5]] < 0.075)
peri_list[[5]] <- as.numeric(names(gene_density_5e5[[5]])[c(79,(82+1))])

which(gene_density_5e5[[6]] < 0.075)
peri_list[[6]] <- as.numeric(names(gene_density_5e5[[6]])[c(79,(88+1))])

which(gene_density_5e5[[7]] < 0.075)
peri_list[[7]] <- as.numeric(names(gene_density_5e5[[7]])[c(51,(53+1))])

which(gene_density_5e5[[8]] < 0.075)
peri_list[[8]] <- as.numeric(names(gene_density_5e5[[8]])[c(49,(51+1))])

which(gene_density_5e5[[9]] < 0.075)
peri_list[[9]] <- as.numeric(names(gene_density_5e5[[9]])[c(51,(56+1))])

which(gene_density_5e5[[10]] < 0.075)
peri_list[[10]] <- as.numeric(names(gene_density_5e5[[10]])[c(57,(66+1))])

which(gene_density_5e5[[11]] < 0.075)
peri_list[[11]] <- as.numeric(names(gene_density_5e5[[11]])[c(42,(50+1))])

which(gene_density_5e5[[12]] < 0.075)
peri_list[[12]] <- as.numeric(names(gene_density_5e5[[12]])[c(47,(53+1))])

which(gene_density_5e5[[13]] < 0.075)
peri_list[[13]] <- as.numeric(names(gene_density_5e5[[13]])[c(28,(34+1))])

which(gene_density_5e5[[14]] < 0.075)
peri_list[[14]] <- as.numeric(names(gene_density_5e5[[14]])[c(23,(29+1))])

which(gene_density_5e5[[15]] < 0.075)
peri_list[[15]] <- as.numeric(names(gene_density_5e5[[15]])[c(70,(73+1))])

which(gene_density_5e5[[16]] < 0.075)
peri_list[[16]] <- as.numeric(names(gene_density_5e5[[16]])[c(47,(58+1))])

which(gene_density_5e5[[17]] < 0.075)
peri_list[[17]] <- as.numeric(names(gene_density_5e5[[17]])[c(64,(75+1))])

which(gene_density_5e5[[18]] < 0.075)
peri_list[[18]] <- as.numeric(names(gene_density_5e5[[18]])[c(69,(76+1))])

peri_tab <- data.table(CHROM = chr_names, 
  PERI_START = unlist(lapply(peri_list, function(x) x[1])),
  PERI_END = unlist(lapply(peri_list, function(x) x[2]))
)

fwrite(peri_tab, peri_out, sep = '\t')

saveRDS(gene_density_5e5, density_out)

### Code used for choosing window size and density cutoff


# 1e6 window size
gene_density_1e6 <- list()
for(chrom in chr_names){
  tmp_bed <- gene_bed[V1 == chrom,]
  gene_density_1e6[[chrom]] <- calculate_window_gene_density(
    test_bed = tmp_bed, window_size = 1e6)
}

# 5e6 window size
gene_density_5e6 <- list()
for(chrom in chr_names){
  tmp_bed <- gene_bed[V1 == chrom,]
  gene_density_5e6[[chrom]] <- calculate_window_gene_density(
    test_bed = tmp_bed, window_size = 5e6)
} 

# 1e7 window size
gene_density_1e7 <- list()
for(chrom in chr_names){
  tmp_bed <- gene_bed[V1 == chrom,]
  gene_density_1e7[[chrom]] <- calculate_window_gene_density(
    test_bed = tmp_bed, window_size = 1e7)
} 

gd_5e5 <- data.table(gene_density = unlist(gene_density_5e5))
gg_5e5 <- ggplot(gd_5e5, aes(x = gene_density)) +
  geom_histogram() +
  ggtitle('Gene density using 5e5 window')

gd_1e6 <- data.table(gene_density = unlist(gene_density_1e6))
gg_1e6 <- ggplot(gd_1e6, aes(x = gene_density)) +
  geom_histogram() +
  ggtitle('Gene density using 1e6 window')

gd_5e6 <- data.table(gene_density = unlist(gene_density_5e6))
gg_5e6 <- ggplot(gd_5e6, aes(x = gene_density)) +
  geom_histogram() +
  ggtitle('Gene density using 5e6 window')

gd_1e7 <- data.table(gene_density = unlist(gene_density_1e7))
gg_1e7 <- ggplot(gd_1e7, aes(x = gene_density)) +
  geom_histogram() +
  ggtitle('Gene density using 1e7 window')


out_pdf <- '/home/f1p1/tmp/switchgrass_8X/sg_gene_density_windows_v1.pdf'

pdf(out_pdf, width = 8, height = 5)
(gg_5e5 + gg_1e6) / (gg_5e6 + gg_1e7)
dev.off()





# 2.5e5 window size
gene_density_2.5e5 <- list()
for(chrom in chr_names){
  tmp_bed <- gene_bed[V1 == chrom,]
  gene_density_2.5e5[[chrom]] <- calculate_window_gene_density(
    test_bed = tmp_bed, window_size = 2.5e5)
}

# I think 2.5e5 is too fine of a window

sum(gene_density_2.5e5[[1]] < 0.1)/length(gene_density_2.5e5[[1]])
# 15.9% of chromosome looks like pericentromer if use 10% as cutoff
15/length(gene_density_2.5e5[[1]])
# only 6.6% of the chromosome forms continuous pericentromeric region with
#  this cutoff
# 22517760 to 26017760 = 3.5 Mb

sum(gene_density_2.5e5[[1]] < 0.05)/length(gene_density_2.5e5[[1]])
# 8.4% of chromosome looks like pericentromere
5/length(gene_density_2.5e5[[1]])
# 2.2% of chromosome forms continuous pericentromeric regions with this cutoff
# 23017760 to 24017760
# 

sum(gene_density_5e5[[1]] < 0.1)/length(gene_density_5e5[[1]])
# 13.1% looks like pericentromere
7/length(gene_density_5e5[[1]])
# 6.1% forms continuous pericentromeric regions

sum(gene_density_5e5[[1]] < 0.05)/length(gene_density_5e5[[1]])
# 7.9% looks like pericentormere
7/length(gene_density_5e5[[1]])
# 6.1% forms continuous pericentromeric regions

sum(gene_density_5e5[[1]] < 0.075)/length(gene_density_5e5[[1]])
# 10.5% looks like pericentormere
7/length(gene_density_5e5[[1]])
# 6.1% forms continuous pericentromeric regions

sum(gene_density_1e6[[1]] < 0.075)/length(gene_density_1e6[[1]])
# 7.0% looks like pericentormere
3/length(gene_density_1e6[[1]])
# 5.3%

sum(gene_density_5e5[[3]] < 0.1)/length(gene_density_5e5[[3]])
# 14.5% looks like pericentromere
10/length(gene_density_5e5[[3]])
# 7.2% OR 11.6%  forms continuous pericentromeric regions

sum(gene_density_5e5[[3]] < 0.05)/length(gene_density_5e5[[3]])
# 10.9% looks like pericentromere
10/length(gene_density_5e5[[3]])
# 7.2% OR 11.6%  forms continuous pericentromeric regions

sum(gene_density_5e5[[16]] < 0.05)/length(gene_density_5e5[[16]])
# 11.7% looks like pericentromere
10/length(gene_density_5e5[[3]])

which(gene_density_2.5e5[[1]] < 0.1)

which(gene_density_5e5[[1]] < 0.05)

which(gene_density_1e6[[1]] < 0.05)

which(gene_density_5e6[[1]] < 0.05)
# too big of a window




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

