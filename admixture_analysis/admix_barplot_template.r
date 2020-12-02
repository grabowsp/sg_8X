# Generate bargraphs for ADMIXTURE results

# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###
library(tidyverse)
library(data.table)

gen_func_file <- '/global/homes/g/grabowsp/tools/sg_8X/general_r_tools/general_functions.r'
source(gen_func_file)

struc_func_file <- '/global/homes/g/grabowsp/tools/sg_8X/structure_analysis/structure_functions/struc_functions.r'
# * change to appropriate path
source(struc_func_file)

### LOAD INPUTS ###
data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/admix_analysis/geobig_admix/'
data_dir <- add_slash(data_dir)

admix_res_short <- 'GW_50k_geobig.3.Q'
admix_res_file <- paste(data_dir, admix_res_short, sep = '')
admix_res <- fread(admix_res_file, header = F)

fam_file_short <- 'GW_50k_geobig.fam'
fam_file <- paste(data_dir, fam_file_short, sep = '')
fam_info <- fread(fam_file, header = F)
samp_names <- unlist(fam_info[, 1])

### SET OUTPUTS ###

out_file <- gsub('.Q', '.ADMIXTURE.memb.pdf', admix_res_file, fixed = T)

### SET VARIABLES ###

base_width <- 12
base_height <- 4

# the amount that the number of samples are divided by to scale the width of
#   the plot
width_factor <- 125

tot_width <- base_width * nrow(admix_res)/width_factor

##############
admix_df <- data.frame(samp_names, admix_res, stringsAsFactors = F)

### evaluate what the groups are and set group plotting variables

# re-order samples based on 'pop_order'
test_samp_orders <- order_struc_samps(clumpp_result_df = 
  admix_df[, c(2:ncol(admix_df))], pop_order = c(1:ncol(admix_res)), 
  zero_cut = 0.1)

#####
res_ord <- admix_df[test_samp_orders, ]

samp_vec <- rep(res_ord[,1], times = ncol(res_ord)-1)
group_vec <- rep(seq(ncol(res_ord)-1), each = nrow(res_ord))

res_vec <- c()
for(j in c(2:ncol(res_ord))){
  res_vec <- c(res_vec, res_ord[,j])
}

plot_df <- data.frame(LIB = samp_vec, GROUP = group_vec, MEMB = res_vec,
  stringsAsFactors = F)

plot_df$LIB <- factor(plot_df$LIB, levels = res_ord[,1])

plot_df$GROUP <- factor(plot_df$GROUP, levels = c(1:ncol(admix_res)))

gg_bar <- ggplot(plot_df, aes(y = MEMB, x = LIB, fill = GROUP)) +
  geom_bar(position = 'fill', stat = 'identity') +
#  scale_fill_manual('Group', values = group_col_vec, 
#    labels = group_lab_vec) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = 'Group Membership')

pdf(out_file, width = tot_width, height = base_height)
gg_bar
dev.off()

quit(save = 'no')

