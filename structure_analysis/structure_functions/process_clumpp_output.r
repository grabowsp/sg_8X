# Script for processing CLUMPP output

# bash
# source activate R_analysis

args <- commandArgs(trailingOnly = TRUE)

### LOAD INPUTS ###
#pre_file <- args[1]
#pre_info <- read.table(pre_file, header = F, stringsAsFactors = F)

res_file <- args[1]
res_info_tmp <- read.table(res_file, header = F, stringsAsFactors = F)

struc_res_1_file <- args[2]
struc_res_1 <- read.table(struc_res_1_file, header = F, stringsAsFactors = F)

### SET OUTPUT ###

out_file <- gsub('.clumpp_', '_', gsub('out', 'clumpp.processed', res_file), 
  fixed = T)

######

nreps <- 3
nsamps <- nrow(struc_res_1)

out_info <- data.frame(LIB = struc_res_1[,1],
  res_info_tmp[seq(nsamps), c(6:ncol(res_info_tmp))], stringsAsFactors = F)

write.table(out_info, out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)

quit(save = 'no')





