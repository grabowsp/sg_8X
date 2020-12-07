#!/bin/bash

# STEPS: 
# 1. Concatenate the Ln Prob Data info from all structure results
# 2. Use Rscript to generate Evanno delta-K statistics, output figure and table

# INPUT 1: prefix of the structure results; example: expandgeo_pseudohap_
# INPUT 2: Maximum K; example: 10
# INPUT 3: Number of structure replicate runs: example: 3

OUT_PRE=$1
MAX_K=$2
N_R=$3

FILE_SHORT=$OUT_PRE'LnProbData.txt'

# Concatenate Ln Prob Data
for KT in {1..10};
  do
  for KR in {1..3};
    do
    grep 'Estimated Ln Prob of Data' $OUT_PRE$KT'.'$KR'_f' >> $FILE_SHORT;
    done;
  done;

source activate R_analysis

Rscript /home/grabowsky/tools/workflows/sg_ha_effort/polyploid_genos/structure_analysis/gen_deltaK_plot.r $FILE_SHORT $MAX_K $N_R $OUT_PRE 


