# Steps for making CDS and genic BED files to be used for filtering VCFs

## Location of files
* CDS .bed
  * `/global/homes/g/grabowsp/data/switchgrass/sg_v5_CDS.bed`
* Genes .bed
  * `/global/homes/g/grabowsp/data/switchgrass/sg_v5_genic.bed`

## Generate CDS .bed file
* in R
```
# module load python/3.7-anaconda-2019.07
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

library(tidyverse)
library(data.table)
# library(bit64)

### INPUT DATA ###
sg_gff_file <- '/global/homes/g/grabowsp/data/switchgrass/ref_stuff/Pvirgatum_516_v5.1/Pvirgatum_516_v5.1.gene.gff3'

sg_gff <- fread(sg_gff_file)
colnames(sg_gff)[1] <- 'V1'

### SET OUTPUT ###

out_dir <- '/global/homes/g/grabowsp/data/switchgrass/'

CDS_out_short <- 'sg_v5_CDS.bed'
CDS_out <- paste(out_dir, CDS_out_short, sep = '')

genic_out_short <- 'sg_v5_genic.bed'
genic_out <- paste(out_dir, genic_out_short, sep = '')

###################

sg_CDS_bed <- sg_gff[V3 == 'CDS', .(V1, V4, V5, V9)]

CDS_names_short <- gsub('ID=', '', unlist(lapply(
  strsplit(sg_CDS_bed$V9, split = ';'), function(x) x[1])))

sg_CDS_bed[ , SHORT_NAME := CDS_names_short]

sg_CDS_bed_df <- data.frame(sg_CDS_bed[, .(V1, V4, V5, SHORT_NAME)], 
  stringsAsFactors = F)

write.table(sg_CDS_bed_df, file = CDS_out, quote = F, sep = '\t', 
  row.names = F, col.names = F)

sg_gene_bed <- sg_gff[V3 == 'gene', .(V1, V4, V5, V9)]

genic_names_short <- gsub('Name=', '', unlist(lapply(
  strsplit(sg_gene_bed$V9, split = ';'), function(x) x[2])))

sg_gene_bed[, SHORT_NAME := genic_names_short]

sg_gene_bed_df <- data.frame(sg_gene_bed[, .(V1, V4, V5, SHORT_NAME)],
  stringsAsFactors = F)

write.table(sg_gene_bed_df, file = genic_out, quote = F, sep = '\t',
  row.names = F, col.names = F)

```
