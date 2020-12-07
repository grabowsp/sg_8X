#!/bin/bash

# STEPS: 
# 1. Concatenate the structure output for the K value
# 2. Generate CLUMPP input
# 3. Run CLUMPP
# 4. Process CLUMPP output

# INPUT 1: prefix of the structure results; example: expandgeo_pseudohap_
# INPUT 2: K value; example: 2
# INPUT 3: Number of samples; example: 785
# INPUT 4: Number of structure replicate runs; example: 3

RES_PRE=$1
K_VAL=$2
N_C=$3
N_R=$4

# Generate CLUMPP input

cat $RES_PRE$K_VAL*'_q' > $RES_PRE$K_VAL'_combo_q_tmp'

source activate R_analysis

LNPROB_IN=$RES_PRE$K_VAL'_combo_q_tmp'

Rscript /home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/gen_clumpp_input.r $LNPROB_IN

# Run CLUMPP

source activate structure_env

PARAMFILE=/home/f2p1/work/grabowsk/data/switchgrass/polyploid_genos/\
struc/generic.clumpp.paramfile

ln -s $PARAMFILE ./paramfile

DATA_PRE=$RES_PRE$K_VAL'_combo.clumpp'
INDFILE=$DATA_PRE'.indfile'
MISCFILE=$DATA_PRE'.miscfile'

OUTFILE=$DATA_PRE'_k'$K_VAL'.out'

CLUMPP paramfile -i $INDFILE -o $OUTFILE -j $MISCFILE -k $K_VAL -c $N_C -r $N_R

# Process results

source activate R_analysis

STR_RES_1=$RES_PRE$K_VAL'.1_q'

Rscript /home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/process_clumpp_output.r $OUTFILE $STR_RES_1

#Rscript /home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/process_clumpp_output.r $INDFILE $OUTFILE $STR_RES_1




