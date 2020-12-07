# Steps for running STRUCTURE on the `geobig` sample set

## Overview
* Using tetrasomic genotypes for all samples
  * previous tests show high correlation between using all-tetrasomic, all-disomic, and all-pseudohaploid genotypes, but a slight ploidy-related bias when using part-NA genotypes 
  * results using pseudohaploid genotypes have less precision than the tetrasomic and disomic results

## Location of Generic param files
* HA Directory
  * `/home/f2p1/work/grabowsk/data/switchgrass/structure_8X`
### mainparam files
* Tetrasomic genotypes
  * `tet_generic.mainparams`
### extraparams
* `generic.extraparams`

## Generate genotypes
* generating on HA while Cori is down for maintenance
* location of 50k SNP 'genlight' genotype file on HA
  * `/home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.rds`
* location of 50k STRUCTURE genotype file
  * `/home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig/geobig_STRUC_genos.txt`
### Run script interactively
* `/home/grabowsky/tools/workflows/sg_8X/adegenet_analysis/adegenet_genotype_generation/generate_structure_input.r`
```
geno_in_file <- '/home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig/GW.50kSNPs.tetrasomic.CDS.geobig.genlight.rds'

out_short <- 'geobig_STRUC_genos.txt'

n_snps <- nLoc(gen_tot)

geno_type <- 'all_tet'
```

## Run 'geobig' K=1 to K=10
### Soft link the parameter files
```
cd /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig

ln -s /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/tet_generic.mainparams ./tet_generic.mainparams
ln -s /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/generic.extraparams ./generic.extraparams
```
### Generate submit scripts
```
cd /home/f2p1/work/grabowsk/data/switchgrass/polyploid_genos/struc/expand_geo/all_dip

bash

for SUB_FILE in expandgeo_structure_*sh;
  do
  FILE_SUF=`echo $SUB_FILE | awk -F'[_]' '{print $4}'`;
  sed 's/t4c1\/WORK\/grabowsk\/data\/switchgrass\/polyploid_genos\/struc\/expand_geo\/all_dip/f2p1\/work\/grabowsk\/data\/switchgrass\/structure_8X\/geobig/g; s/expandgeo_25k.strucgenos.txt_all_dip/geobig_STRUC_genos.txt/g; s/OUT_NAME=expandgeo_alldip_/OUT_NAME=geobig_alltet_/g; s/N_SAMPS=785/N_SAMPS=821/g; s/N_SNPS=25000/N_SNPS=50000/g; s/MAIN_PARAM=dip_generic.mainparams/MAIN_PARAM=tet_generic.mainparams/g; s/-N alldip_/-N geobig_tet_/g' \
  $SUB_FILE > /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig/geobig_alltet_structure_$FILE_SUF;
  done 
```
### Submit jobs
```
cd /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig

bash

for KT in {1..10};
  do
for KR in {1..3};
    do
    qsub geobig_alltet_structure_k$KT'.'$KR'.sh';
  done;
done
```

## Get Estimated Ln Prob of Data
* uses shell script to compile LnProb data and R script to generate figure
* Ln Prob data table
  * `/home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig/geobig_alltet_LnProbData_deltaK.txt`
* delta K figure
  * `/home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig/geobig_alltet_LnProbData_deltaK.pdf`
### Run interactively
```
bash

cd /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig

OUT_PRE=geobig_alltet_
MAX_K=10
N_R=3

/home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/run_deltaK_steps.sh $OUT_PRE $MAX_K $N_R
```

## Run CLUMPP
### K=2
```
bash

cd /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig

RES_PRE=geobig_alltet_
K_VAL=2
N_C=821
N_R=3

/home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

```
### K=3
```
bash

cd /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig

RES_PRE=geobig_alltet_
K_VAL=3
N_C=821
N_R=3

/home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

```
### K=4
```
bash

cd /home/f2p1/work/grabowsk/data/switchgrass/structure_8X/geobig

RES_PRE=geobig_alltet_
K_VAL=4
N_C=821
N_R=3

/home/grabowsky/tools/workflows/sg_8X/structure_analysis/structure_functions/run_clumpp_steps.sh $RES_PRE $K_VAL $N_C $N_R

```
