

head_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tot_VCF_header.txt'

head_txt <- read.table(head_file, header = F, comment.char = '$', 
  stringsAsFactors = F, sep = '\t')

contig_vec <- head_txt[grep(',length=', head_txt[,1], fixed = T), 1]

contig_names <- gsub(',length', '', 
  unlist(lapply(strsplit(contig_vec, split = '='), function(x) x[3])),
  fixed = T)

contig_df <- data.frame(contig = contig_names, 
  c_num = seq(length(contig_names)), stringsAsFactors = F)

out_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/plink_chr_name_map.txt'

write.table(contig_df, file = out_file, quote = F, sep = '\t', row.names = F,
  col.names = F)

quit(save = 'no')

