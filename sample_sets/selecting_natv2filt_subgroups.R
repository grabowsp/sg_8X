### LOAD MODULES ###
library(data.table)
library(ggplot2)

### INPUT DATA ###
samp_meta_file <- '/Users/grabowsk/Analysis/sg_8X/metadata_management/meta_files/sg_8X_metadata_v3.1.csv'
samp_meta <- fread(samp_meta_file)

res_tab_file <- '/Users/grabowsk/data/sg_8X_analysis/natv2filt_res_tab_v1.0.txt'
res_tab <- fread(res_tab_file)

res_tab$pch <- 0
res_tab[which(res_tab$ploidy == '8X'), pch := 1]
res_tab[which(res_tab$ploidy == '6X'), pch := 2]

res_tab[ , tmp_col := 'black']
res_tab[which(res_tab$full_admix_k3_MW < 0.7 & 
                res_tab$full_admix_k3_GULF < 0.7 &
                res_tab$full_admix_k3_ATL < 0.65),
        tmp_col := 'red']

gg_tmp <- ggplot(res_tab, aes(x = full_pc01_raw, y = full_pc02_raw, 
                              color = res_tab$tmp_col)) +
  geom_point(shape = res_tab$pch)

gg_tmp

# For MW > 0.9 (K=2 or K=3)
# res_tab$full_admix_k2_MW > 0.9
# res_tab$full_admix_k3_MW > 0.9
# J581.A is a weird sample that is > 0.9 MW with admixture but clusters
#  clusters on its own in PCA - I'd say DONT include this one in the MW cluster

# K=2 LOW > 0.9
# res_tab$full_admix_k2_LOW > 0.9
# Most of the Gulf not included with this group...

# > 0.9 for Gulf or Atlantic
# which(res_tab$full_admix_k3_ATL > 0.9 | res_tab$full_admix_k3_GULF > 0.9)
# better than K=2 LOW, but still not perfect...

##################

# The "main clusters" may not actually be true clusters - will have to see
#   what they look like with PCAs of subgroups

# For clustering within MW:
# 1) Select the "main MW cluster" based on PC1 and PC2 locations
### which((res_tab$full_pc01_raw > 22 & res_tab$full_pc02_raw > 2.5) | (res_tab$full_pc01_raw > 26 & res_tab$full_pc02_raw > 1))
# 2) Select all potential MW: 
### which(res_tab$full_admix_k3_MW > 0.7)
# 3) Pull out J581.A as it's own group

# For clustering within Gulf
# 1) Select "main Gulf cluster: based on PC1 and PC2
### which(res_tab$full_pc01_raw < -3 & res_tab$full_pc02_raw < -20)
# 2) Select all potential Gulf:
### which(res_tab$full_admix_k3_GULF > 0.7)

# For clustering withing Atlantic
# 1) Select "main Atlantic cluster" based on PC1 and PC2
### which(res_tab$full_pc01_raw < -15 & res_tab$full_pc02_raw > 4)
# 2) Select all potential Atlantic
### which(res_tab$full_admix_k3_ATL > 0.65)
# 3) Pull out J534.A (77% ATL, 22% MW) and J530.B (69% ATL, 29% MW, 2% Gulf)

# For clustering remaining samples
# 1) Samples that don't fall in any "potential" gene pool
### which(res_tab$full_admix_k3_MW < 0.7 & res_tab$full_admix_k3_GULF < 0.7 &
#   res_tab$full_admix_k3_ATL < 0.65)
# 2) Add outliers from other genepool cutoffs: J581.A, J534.A, and J530.B




