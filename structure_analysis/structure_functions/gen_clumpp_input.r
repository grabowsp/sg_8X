# Script for generating CLUMPP inpute

# bash
# source activate R_analysis

args <- commandArgs(trailingOnly = TRUE)

### LOAD INPUTS ###
in_file <- args[1]
in_res <- read.table(in_file, header = F, stringsAsFactors = F)

### SET OUTPUT ###

out_file <- gsub('_q_tmp', '.clumpp.indfile', in_file)

######

nreps <- 3
nsamps <- nrow(in_res)/3

samp_int_vec <- rep(seq(nsamps), times = 3)

third_vec <- rep('(x)', times = nrow(in_res))

pop_vec <- rep(NA, times = nrow(in_res))
fifth_vec <- rep(':', times = nrow(in_res))

out_res_tmp <- data.frame( C1=samp_int_vec, C2 = samp_int_vec,
  C3 = third_vec, C4 = samp_int_vec, C5 = fifth_vec,
  stringsAsFactors = F)

out_res <- cbind(out_res_tmp, in_res[, c(2:ncol(in_res))])

write.table(out_res, out_file, quote=F, sep = '\t', row.names = F,
  col.names = F)

quit(save = 'no')





