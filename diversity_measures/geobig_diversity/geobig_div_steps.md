# Steps for measuring diversity measures in `geo_big` samples

## Overview
* NERSC Resulst directory
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity`

## Calculate homozygosity in CDS across chromosomes
### Run script
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity

sbatch calc_geobig_div.sh
```
### Example of script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 48:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/adegenet_2_env

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity
DATA_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/geobig_tet_vcfs/
DATA_FILE_SUF=.tetrasomic.CDS.geobig.genlight.rds

OUT_PRE_SUB=geobig.CDS

BLOCKSIZE=1e5

cd $OUT_DIR

for CHR_NUM in {01..09};
  do
  for CHR_LET in K N;
    do
    DATA_FILE=$DATA_DIR'Chr'$CHR_NUM$CHR_LET$DATA_FILE_SUF;
    OUT_PRE=Chr$CHR_NUM$CHR_LET'.'$OUT_PRE_SUB;
    Rscript /global/homes/g/grabowsp/tools/sg_8X/diversity_measures/calc_hom_het_vals.r \
    $DATA_FILE $OUT_DIR $OUT_PRE $BLOCKSIZE;
    done;
  done;
```

## Aggregate homozygosity for all chromosomes
### Test
```
module load python/3.7-anaconda-2019.07
source activate adegenet_2_env

data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/diversity_analysis/geobig_diversity/'

file_suf <- '.geobig.CDS.HETvals.rds'

res_file_vec <- system(paste('ls ', data_dir, '*HETvals.rds', sep = ''),
  intern = T)

tet_list <- list()
oct_list <- list()
oct_cor_list <- list()
snp_list <- list()

for(rfv in res_file_vec){
  tmp_name <- sub(data_dir, '', sub('.geobig.CDS.HETvals.rds', '', rfv, 
    fixed = T))
  tmp_res <- readRDS(rfv)
  tet_list[[tmp_name]] <- tmp_res[['tot_tet_hom']]
  oct_list[[tmp_name]] <- tmp_res[['tot_oct_hom']]
  oct_cor_list[[tmp_name]] <- tmp_res[['tot_correct_oct_het']]
  snp_list[[tmp_name]] <- tmp_res[['N_SNPs']]
}

tot_snps <- sum(unlist(snp_list))

tet_list_alt <- list()
oct_list_alt <- list()
oct_cor_list_alt <- list()

for(nchrom in seq(length(snp_list))){
  tmp_name <- names(snp_list)[nchrom]
  tet_list_alt[[tmp_name]] <- tet_list[[tmp_name]]*snp_list[[tmp_name]]
  oct_list_alt[[tmp_name]] <- oct_list[[tmp_name]]*snp_list[[tmp_name]]
  oct_cor_list_alt[[tmp_name]] <- oct_cor_list[[tmp_name]]*snp_list[[tmp_name]]
}

tet_sum_mat <- tet_list_alt[[1]]
oct_sum_mat <- oct_list_alt[[1]]
oct_cor_sum_mat <- oct_cor_list_alt[[1]]

for(i in c(2:length(snp_list))){
  tet_sum_mat <- tet_sum_mat + tet_list_alt[[i]]
  oct_sum_mat <- oct_sum_mat + oct_list_alt[[i]]
  oct_cor_sum_mat <- oct_cor_sum_mat + oct_cor_list_alt[[i]]
}

tet_hom_norm <- tet_sum_mat / tot_snps
oct_hom_norm <- oct_sum_mat / tot_snps
oct_cor_het_norm <- oct_cor_sum_mat / tot_snps

raw_het <- c((1-tet_hom_norm), (1-oct_hom_norm))
cor_het <- c((1-tet_hom_norm), oct_cor_het_norm)
samp_name_vec <- c(names(tet_hom_norm), names(oct_hom_norm))

het_df <- data.frame(samp_name = samp_name_vec, het_raw = raw_het, 
  het_corrected = cor_het, stringsAsFactors = F)

out_file <- paste(data_dir, 'geobig.heterozygosity_v1.txt', sep = '')
write.table(het_df, file = out_file, quote = F, sep = '\t', 
  row.names = F, col.names = T)

```


