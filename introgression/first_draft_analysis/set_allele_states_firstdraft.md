# Steps for setting the allele states for first draft analysis

## Steps
1. Calculate allele freq and HW frequency in training sets
2. Use R script to set allele states
3. Explore results

## Output files
### Explanation
* File types:
 * `..allele_states.txt` has the allele state designation for each SNP
 * `..keep_pos.txt` includes the SNPs that passed all the filtering of the r script
  * in theory, no SNPs should have been removed because of previous filtering
 * `..ref_freq.txt` contains the allele frequency for each training set 
### Files
* HAIB Directory
 * `/home/f1p1/tmp/switchgrass_8X/introgression_files/`
* NERSC Directory
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/`
* Midwest vs Atlantic
 * Fst > 0.5
  * `MW_v_ATL_Fst0.5_allele_states.txt`
  * `MW_v_ATL_Fst0.5_keep_pos.txt`
  * `MW_v_ATL_Fst0.5_ref_freq.txt`
 * Fst = 1.0
  * `MW_v_ATL_Fst1.0_allele_states.txt`
  * `MW_v_ATL_Fst1.0_keep_pos.txt`
  * `MW_v_ATL_Fst1.0_ref_freq.txt`
* Midwest vs Gulf
 * Fst > 0.5
  * `MW_v_GULF_Fst0.5_allele_states.txt`
  * `MW_v_GULF_Fst0.5_keep_pos.txt`
  * `MW_v_GULF_Fst0.5_ref_freq.txt`
 * Fst = 1.0
  * `MW_v_GULF_Fst1.0_allele_states.txt`
  * `MW_v_GULF_Fst1.0_keep_pos.txt`
  * `MW_v_GULF_Fst1.0_ref_freq.txt`
* Atlantic vs Gulf
 * Fst > 0.5
  * `ATL_v_GULF_Fst0.5_allele_states.txt`
  * `ATL_v_GULF_Fst0.5_keep_pos.txt`
  * `ATL_v_GULF_Fst0.5_ref_freq.txt`
 * Fst = 1.0
  * `ATL_v_GULF_Fst1.0_allele_states.txt`
  * `ATL_v_GULF_Fst1.0_keep_pos.txt`
  * `ATL_v_GULF_Fst1.0_ref_freq.txt`

## Summary of Allele States
* From R analysis below
 * N_BOTH = number of SNPs where both alleles are informative
 * N_POP1 = number of SNPs with only a POP1-informative allele
 * N_POP2 = number of SNPs with only a POP2-informative allele
 * N_NO_INFO = number of SNPs where neither allele is informative
  * This is because of the cutoff I am using for calling informative alleles
  * I could adjust this, but the number of SNPs is very small
```
   POP1 POP2 FST N_BOTH N_POP1 N_POP2 N_NO_INFO
1: GULF  ATL 0.5  10160  16464  11896        11
2: GULF  ATL 1.0    100      0      0         0
3:   MW  ATL 0.5  30715  31194  27207        27
4:   MW  ATL 1.0   1674      0      0         0
5: GULF   MW 0.5  27502  29936  24013        27
6: GULF   MW 1.0    718      0      0         0
```

## Generate HW and allele freq files from vcfTools
```
# at NERSC
module load python/3.7-anaconda-2019.07
source activate gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files/

cd $OUT_DIR

COMP=AtlanticvsGulf

for FST in 0.5 1.0;
 do
 VCF_IN=$VCF_DIR$COMP'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.vcf.gz';
 for TRAIN_SET in GULF ATL;
  do
  OUT_PRE=$TRAIN_SET'_train_'$COMP'_firstdraft_Fst'$FST;
  SAMP_FILE=$TRAIN_SET'_firstdraft_train_filt40.txt';
  #
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $SAMP_FILE --freq;
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $SAMP_FILE --hardy;
  done;
 done;

#

COMP=AtlanticvsMidwest
for FST in 0.5 1.0;
 do
 VCF_IN=$VCF_DIR$COMP'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.vcf.gz';
 for TRAIN_SET in MW ATL;
  do
  OUT_PRE=$TRAIN_SET'_train_'$COMP'_firstdraft_Fst'$FST;
  SAMP_FILE=$TRAIN_SET'_firstdraft_train_filt40.txt';
  #
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $SAMP_FILE --freq;
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $SAMP_FILE --hardy;
  done;
 done;

#

COMP=GulfvsMidwest
for FST in 0.5 1.0;
 do
 VCF_IN=$VCF_DIR$COMP'.hiFst.firstdraft.Fst'$FST'.disomic.CDS.vcf.gz';
 for TRAIN_SET in MW GULF;
  do
  OUT_PRE=$TRAIN_SET'_train_'$COMP'_firstdraft_Fst'$FST;
  SAMP_FILE=$TRAIN_SET'_firstdraft_train_filt40.txt';
  #
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $SAMP_FILE --freq;
  vcftools --gzvcf $VCF_IN --out $OUT_PRE --keep $SAMP_FILE --hardy;
  done;
 done;

```

## Set allele states
```
module load python/3.7-anaconda-2019.10
module swap PrgEnv-intel PrgEnv-gnu
source activate R_tidy

HW_FILE_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

MISS_CUT=0.5
N_TRAIN=40
MIN_KEEP_FREQ=0.5

cd $OUT_DIR

COMP=AtlanticvsMidwest
POP1=MW
POP2=ATL
#FST=0.5

for FST in 0.5 1.0;
 do
 POP1_HW_FILE=$HW_FILE_DIR$POP1'_train_'$COMP'_firstdraft_Fst'$FST'.hwe';
 POP2_HW_FILE=$HW_FILE_DIR$POP2'_train_'$COMP'_firstdraft_Fst'$FST'.hwe';
 #
 Rscript /global/homes/g/grabowsp/tools/sg_8X/introgression/set_introgression_allele_states.r \
  $POP1_HW_FILE $POP2_HW_FILE $POP1 $POP2 $OUT_DIR \
  $MISS_CUT $N_TRAIN $MIN_KEEP_FREQ;
 #
 mv $POP1'_v_'$POP2'_allele_states.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_allele_states.txt';
 #
 mv $POP1'_v_'$POP2'_keep_pos.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_keep_pos.txt';
 #
 mv $POP1'_v_'$POP2'_ref_freq.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_ref_freq.txt';
 done;

COMP=GulfvsMidwest
POP1=MW
POP2=GULF

for FST in 0.5 1.0;
 do
 POP1_HW_FILE=$HW_FILE_DIR$POP1'_train_'$COMP'_firstdraft_Fst'$FST'.hwe';
 POP2_HW_FILE=$HW_FILE_DIR$POP2'_train_'$COMP'_firstdraft_Fst'$FST'.hwe';
 #
 Rscript /global/homes/g/grabowsp/tools/sg_8X/introgression/set_introgression_allele_states.r \
  $POP1_HW_FILE $POP2_HW_FILE $POP1 $POP2 $OUT_DIR \
  $MISS_CUT $N_TRAIN $MIN_KEEP_FREQ;
 #
 mv $POP1'_v_'$POP2'_allele_states.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_allele_states.txt';
 #
 mv $POP1'_v_'$POP2'_keep_pos.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_keep_pos.txt';
 #
 mv $POP1'_v_'$POP2'_ref_freq.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_ref_freq.txt';
 done;

COMP=AtlanticvsGulf
POP1=ATL
POP2=GULF

for FST in 0.5 1.0;
 do
 POP1_HW_FILE=$HW_FILE_DIR$POP1'_train_'$COMP'_firstdraft_Fst'$FST'.hwe';
 POP2_HW_FILE=$HW_FILE_DIR$POP2'_train_'$COMP'_firstdraft_Fst'$FST'.hwe';
 #
 Rscript /global/homes/g/grabowsp/tools/sg_8X/introgression/set_introgression_allele_states.r \
  $POP1_HW_FILE $POP2_HW_FILE $POP1 $POP2 $OUT_DIR \
  $MISS_CUT $N_TRAIN $MIN_KEEP_FREQ;
 #
 mv $POP1'_v_'$POP2'_allele_states.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_allele_states.txt';
 #
 mv $POP1'_v_'$POP2'_keep_pos.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_keep_pos.txt';
 #
 mv $POP1'_v_'$POP2'_ref_freq.txt' \
  $POP1'_v_'$POP2'_Fst'$FST'_ref_freq.txt';
 done;
```

## Explore results
```
module load python/3.7-anaconda-2019.10
module swap PrgEnv-intel PrgEnv-gnu
source activate R_tidy

library(data.table)

data_dir <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/'

allele_state_files <- system(paste('ls ', data_dir, '*allele_states.txt', 
  sep = ''), intern = T)

tmp_data <- fread(allele_state_files[1])

tot_types <- setdiff(unique(c(tmp_data$REF_STATE, tmp_data$ALT_STATE)), 
  'NO_INFO')
pop_names <- unique(gsub('_MAINLY|_ONLY', '', tot_types))

no_info_inds <- intersect(which(tmp_data$REF_STATE == 'NO_INFO'),
  which(tmp_data$ALT_STATE == 'NO_INFO'))

both_informative_inds <- setdiff(seq(nrow(tmp_data)), 
  unique(c(which(tmp_data$REF_STATE == 'NO_INFO'),
    which(tmp_data$ALT_STATE == 'NO_INFO'))))

pop1_inds <- setdiff(sort(c(grep(pop_names[1], tmp_data$REF_STATE), 
  grep(pop_names[1], tmp_data$ALT_STATE))), both_informative_inds)

pop2_inds <- setdiff(sort(c(grep(pop_names[2], tmp_data$REF_STATE), 
  grep(pop_names[2], tmp_data$ALT_STATE))), both_informative_inds)
 

fst_val <- substr(rev(unlist(
  strsplit(gsub('_allele_states.txt', '', allele_state_files[1]), split = '_')
  ))[1], 4,6)

res_list <- list()

for(i in seq(allele_state_files)){
  tmp_data <- fread(allele_state_files[i])
  tot_types <- setdiff(unique(c(tmp_data$REF_STATE, tmp_data$ALT_STATE)), 
    'NO_INFO')
  pop_names <- unique(gsub('_MAINLY|_ONLY', '', tot_types))
  pop1 <- pop_names[1]
  pop2 <- pop_names[2]
  no_info_inds <- intersect(which(tmp_data$REF_STATE == 'NO_INFO'),
    which(tmp_data$ALT_STATE == 'NO_INFO'))
  both_informative_inds <- setdiff(seq(nrow(tmp_data)), 
    unique(c(which(tmp_data$REF_STATE == 'NO_INFO'),
    which(tmp_data$ALT_STATE == 'NO_INFO'))))
  pop1_inds <- setdiff(sort(c(grep(pop1, tmp_data$REF_STATE), 
    grep(pop1, tmp_data$ALT_STATE))), both_informative_inds)
  pop2_inds <- setdiff(sort(c(grep(pop2, tmp_data$REF_STATE), 
    grep(pop2, tmp_data$ALT_STATE))), both_informative_inds)
  fst_val <- substr(rev(unlist(
    strsplit(gsub('_allele_states.txt', '', allele_state_files[i]), split = '_')
    ))[1], 4,6)
  tmp_dt <- data.table(POP1 = pop1, POP2 = pop2, FST = fst_val, 
    N_BOTH = length(both_informative_inds), N_POP1 = length(pop1_inds),
    N_POP2 = length(pop2_inds), N_NO_INFO = length(no_info_inds))
  res_list[[i]] <- tmp_dt
}

res_tab <- rbindlist(res_list)

res_tab
#   POP1 POP2 FST N_BOTH N_POP1 N_POP2 N_NO_INFO
#1: GULF  ATL 0.5  10160  16464  11896        11
#2: GULF  ATL 1.0    100      0      0         0
#3:   MW  ATL 0.5  30715  31194  27207        27
#4:   MW  ATL 1.0   1674      0      0         0
#5: GULF   MW 0.5  27502  29936  24013        27
#6: GULF   MW 1.0    718      0      0         0

```

