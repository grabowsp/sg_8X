# Functions for VCF analysis

read_vcf <- function(vcf_file){
  # Function for inputing a VCF into R for analysis
  #   mainly for inputting VCFs generated for switchgrass analysis; not sure
  #   about the universality of these functions
  #   - does NOT work for .gz files - those need to be uncompressed first
  # INPUTS #
  # vcf_file = full path to the vcf file; must already be uncompressed
  # OUTPUT #
  # data.table with VCF info, each column still in the original formatting
  ###########
  vcf_in <- read.table(vcf_file, header = F, stringsAsFactors = F)
  vcf_head_tmp <- system(paste('grep CHR ', vcf_file, sep = ''),
    intern = T)
  vcf_head_2 <- sub('#CHROM', 'CHROM', vcf_head_tmp)
  vcf_head_3 <- unlist(strsplit(vcf_head_2, split = '\t'))
  colnames(vcf_in) <- vcf_head_3
  vcf_in <- data.table(vcf_in)
  return(vcf_in)
}

vcf_SNP_genotypes <- function(vcf, snp_name, samp_vec = c()){
  # Look at the genotypes at a specific SNP
  #  - assumes "FORMAT" is the last column before genotypes start
  # INPUTS #
  # vcf = vcf loaded into R
  # snp_name = name of SNP; usually in format of 'CHRNAME_CHRPOS'
  # samp_vec = vector of sample names. If blank, then returns genotypes for
  #  all samples
  # OUTPUT #
  # Vector with genotype, in character format, for each sample
  ######
  # if looking at all samples
  if(length(samp_vec) == 0){
    format_col <- which(colnames(vcf) == 'FORMAT')
    samp_vec <- colnames(vcf)[(format_col + 1):ncol(vcf)]
  }
  if(length(samp_vec) > 0){
    miss_names <- setdiff(samp_vec, colnames(vcf))
    if(length(miss_names) > 0){
      stop(paste(paste(miss_names, collapse = ','), 'missing from VCF'))
    }
  }
  vcf_ind <- which(vcf$ID == snp_name)
  tmp_geno_string <- vcf[vcf_ind, samp_vec, with = F]
  tmp_geno_list <- sapply(tmp_geno_string, strsplit, split = ':')
  tmp_geno_vec <- unlist(lapply(tmp_geno_list, function(x) x[1]))
  return(tmp_geno_vec)
}

vcf_many_SNP_genotypes <- function(vcf, snp_name_vec, samp_vec = c()){
  # Get genotypes for samples at multiple SNPs
  # INPUTS #
  # vcf = vcf loaded into R
  # snp_name_vec = vector of snp_names to be included, ; usually in format 
  #  of 'CHRNAME_CHRPOS'
  # samp_vec = vector of sample names. If blank, then returns genotypes for
  #  all samples
  # OUTPUT #
  # data.table of genotypes, first column = ID (snp_name), other columns
  #  are each sample in samp_vec
  ##########
  # check samp_vec
  if(length(samp_vec) == 0){
    format_col <- which(colnames(vcf) == 'FORMAT')
    samp_vec <- colnames(vcf)[(format_col + 1):ncol(vcf)]
  }
  if(length(samp_vec) > 0){
    miss_names <- setdiff(samp_vec, colnames(vcf))
    if(length(miss_names) > 0){
      stop(paste(paste(miss_names, collapse = ','), 'missing from VCF'))
    }
  }
  miss_snps <- setdiff(snp_name_vec, vcf$ID)
  if(length(miss_snps) > 0){
    stop(paste(paste(miss_snps, collapse = ','), 'missing from VCF'))
  }
  vcf_inds <- which(vcf$ID %in% snp_name_vec)
  geno_df <- apply(vcf[vcf_inds, samp_vec, with = F], 2, function(x)
    unlist(lapply(strsplit(x, split = ':'), function(x) x[1]))) 
  geno_tab <- data.table(ID = vcf$ID[vcf_inds], geno_df)
  return(geno_tab)
}

get_portion_genotype_1SNP <- function(geno_vec, group_info, geno_classes = c(),
    return_percent = T){
  # Calculate the portion of each "group" with each genotype class
  #####
  # check group_info
  miss_g_samp <- setdiff(names(group_info), names(geno_vec))
  if(length(miss_g_samp) > 0){
    stop(paste(paste(miss_g_samp, sep = ','), 'missing from group info'))
  }
  # check geno_classes
  if(length(geno_classes) == 0){
    geno_classes <- c('0/0', '0/1', '1/1', './.')
    names(geno_classes) <- c('homRef', 'het', 'homAlt', 'missing')
  }
  if(length(names(geno_classes)) == 0){
    names(geno_classes) <- paste('geno_class', seq(length(geno_classes)), 
      sep = '_')
  }
  group_vec <- group_info[names(geno_vec)]
  # next - tally up each genotype class by the groups 
  uni_info <- sort(unique(group_vec))
  geno_class_list <- list()
  for(j in seq(length(geno_classes))){
    tmp_count <- sapply(uni_info, function(x) length(intersect(
      which(geno_vec == geno_classes[j]), which(group_vec == x))))
    geno_class_list[[j]] <- tmp_count
  }
  geno_class_count_tab <- data.table(matrix(unlist(geno_class_list), 
    ncol = length(geno_class_list), byrow = F))
  group_count <- apply(geno_class_count_tab, 1, sum)
  geno_class_per_tab <- geno_class_count_tab / group_count
  colnames(geno_class_per_tab) <- names(geno_classes)
  colnames(geno_class_count_tab) <- names(geno_classes)
  if(return_percent){
    per_tab <- data.table(group = uni_info, n_samps = group_count, 
      geno_class_per_tab)
    return(per_tab)
  } else{
    count_tab <- data.table(group = uni_info, n_samps = group_count,
      geno_class_count_tab)
    return(count_tab)
  }
}

get_portion_genotypes <- function(geno_tab, group_info, geno_classes = c()){
  # Calculate the percentage of each group with each genotype class across
  #   many SNPs
  #######
  # check group_info
  miss_g_samp <- setdiff(names(group_info), names(geno_vec))
  if(length(miss_g_samp) > 0){
    stop(paste(paste(miss_g_samp, sep = ','), 'missing from group info'))
  }
  # check geno_classes
  if(length(geno_classes) == 0){
    geno_classes <- c('0/0', '0/1', '1/1', './.')
    names(geno_classes) <- c('homRef', 'het', 'homAlt', 'missing')
  }
  if(length(names(geno_classes)) == 0){
    names(geno_classes) <- paste('geno_class', seq(length(geno_classes)),
      sep = '_')
  }
  samp_geno_count_mat <- apply(geno_tab[, c(2:ncol(geno_tab)), with = F], 2, 
    function(x) sapply(geno_classes, function(y) sum(x == y)))
  group_vec <- group_info[colnames(samp_geno_count_mat)]
  uni_info <- sort(unique(group_vec))
  info_count_list <- list()
  for(ui in uni_info){
    if(sum(group_vec == ui) == 1){
      info_count_list[[ui]] <- samp_geno_count_mat[, which(group_vec == ui)]
    } else{
      info_count_list[[ui]] <- apply(
        samp_geno_count_mat[, which(group_vec == ui)], 1, sum)
    } 
  }
  info_count_mat <- matrix(unlist(info_count_list), 
    ncol = length(info_count_list), byrow = F)
  colnames(info_count_mat) <- names(info_count_list)
  rownames(info_count_mat) <- names(info_count_list[[1]])
  info_count_mat_1 <- t(apply(info_count_mat/nrow(geno_tab), 2, function(x) 
    x / sum(x)))
  group_count <- sapply(rownames(info_count_mat_1), function(x) 
    sum(group_vec == x))
  tot_count_tab <- data.table(group = rownames(info_count_mat_1), 
    n_samps = group_count, info_count_mat_1)
  return(tot_count_tab)
}
