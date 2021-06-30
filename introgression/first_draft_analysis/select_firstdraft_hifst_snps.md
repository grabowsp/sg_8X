# Steps for Selecting high-Fst SNPs for updated analysis for first draft

## Steps
1. Select Training Sample Sets
2. Calcluate Fst between training sets
3. Choose SNPs based on Fst cutoff
4. Prune hiFst SNPs by LD
5. Filter hiFst SNPs by missing data in sample sets
6. Determine ancestry-informative alleles

## Select Training Sample Sets
### Criteria for training sets
1. Order samples based on PC results
  * ex: highest PC1 = "most MW"
2. Select top 60 samples that are also
  1. >99% genepool ancestry in ADMIXTURE results
  2. 4X  
3. Choose 40 samples selecting for maximum accessions and states
4. For bootstrapping:
  1. Select 40 samples
  2. Analyze the other 20 samples
  3. Repeat X number of times
### Sample set files
* R script used to choose samples and generate files:
  `~/sg_8X/introgression/first_draft_analysis/select_firstdraft_trainingsets.r`

## Calculate Fst between training sets
### Overview
* MAF
  * at least 5 samples = 5/(80*2) = 0.03125
* Fst values
  * use 'weir-fst-pop' flag in VCFtools
  * Going to start with Fst >= 0.5
  * Also make file of SNPs with Fst = 1
### Location of files
* Directory with Fst files:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft`
  * Raw Fst output file names end with `...Fst_CHROMOSOME_firstdraft.weir.fst`
  * non-pruned high Fst files end with `...Fst0.5.txt` or `...Fst1.0.txt` 
### Submit script to calculate Fst with VCFtools
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft

sbatch calc_fst_01_03_firstdraft.sh
sbatch calc_fst_04_06_firstdraft.sh
sbatch calc_fst_07_09_firstdraft.sh
```
#### Example Script
```
#!/bin/bash
#SBATCH -D .
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mem=32G
#SBATCH --qos=genepool_shared
#SBATCH -t 24:00:00
#SBATCH --mail-user pgrabowski@hudsonalpha.org
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END

module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/
ATLANTIC_FILE=ATL_firstdraft_train_filt40.txt
GULF_FILE=GULF_firstdraft_train_filt40.txt
MW_FILE=MW_firstdraft_train_filt40.txt

for CHR_NUM in {01..03};
  do
  for CHR_LET in K N;
    do
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    #
    MAF=0.03125;
    OUT_TXT=GULFvMW_Fst_Chr$CHR_NUM$CHR_LET'_firstdraft';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$GULF_FILE --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    OUT_TXT=ATLANTICvMW_Fst_Chr$CHR_NUM$CHR_LET'_firstdraft';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE --weir-fst-pop $SAMP_DIR$MW_FILE;
    #
    OUT_TXT=ATLANTICvGULF_Fst_Chr$CHR_NUM$CHR_LET'_firstdraft';
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT --maf $MAF \
      --weir-fst-pop $SAMP_DIR$ATLANTIC_FILE \
      --weir-fst-pop $SAMP_DIR$GULF_FILE;
    done;
  done;
```

## Select High-Fst SNPs
* Directory with Fst files:
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft`
 * non-pruned high Fst files end with `...Fst0.5.txt` or `...Fst1.0.txt` 
### Select SNPs based on Fst
* R script for filtering SNPs and makeing files:
  * `~/sg_8X/introgression/first_draft_analysis/select_Fst_SNPs_firstdraft.r`
#### SNP info and file locations
* Midwest-vs-Atlantic
  * 167,400 with Fst >= 0.5
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt`
  * 8,435 with Fst = 1
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt`
* Midwest-vs-Gulf
  * 149,467 with Fst >= 0.5
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt`
  * 3,848 with Fst = 1
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt`
  * Noticably fewer Fst=1 SNPs than for MWvsATL
    * Training set not as "selective" as desired
      * Many "pure" Gulf seem to be cultivars and were removed
    * One idea is to generate training set including cultivar samples
    * Gulf and MW might also just share more alleles due to hybridization
    * For paper, it might be better to focus on MW-vs-Atlantic
* Atlantic-vs-Gulf
  * 65,318 with Fst >= 0.5
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.txt`
  * 569 with Fst = 1
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.txt`
### LD-prune hiFst SNPs
#### Generate sample file
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft

SAMP_FILE=/global/homes/g/grabowsp/tools/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt

paste -d '\t' $SAMP_FILE $SAMP_FILE > natv2_samps.plink.txt 
```
#### Generate SNP list for plink
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft

# Atlantic vs Gulf

SNP_FILE=AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.txt
SNP_OUT=AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

SNP_FILE=AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.txt
SNP_OUT=AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

# Atlantic vs Midwest

SNP_FILE=AtlanticVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt
SNP_OUT=AtlanticVsMW_hiFst_SNPs_firstdraft_Fst0.5.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

SNP_FILE=AtlanticVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt
SNP_OUT=AtlanticVsMW_hiFst_SNPs_firstdraft_Fst1.0.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

# Gulf vs Midwest

SNP_FILE=GulfVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt
SNP_OUT=GulfVsMW_hiFst_SNPs_firstdraft_Fst0.5.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

SNP_FILE=GulfVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt
SNP_OUT=GulfVsMW_hiFst_SNPs_firstdraft_Fst1.0.plink.txt
cut -f 1 $SNP_FILE > tmp_chroms.txt
cut -f 2 $SNP_FILE > tmp_pos.txt
paste -d '_' tmp_chroms.txt tmp_pos.txt > $SNP_OUT

```
#### LD-prune using plink
* Atlantic vs Gulf
 * Fst 0.5 LD-pruned file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.LD.txt`
 * Fst 1.0 LD-pruned file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.LD.txt`
* Atlantic vs Midwest
 * Fst 0.5 LD-pruned file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMidwest_hiFst_SNPs_firstdraft_Fst0.5.LD.txt`
 * Fst 1.0 LD-pruned file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMidwest_hiFst_SNPs_firstdraft_Fst1.0.LD.txt`
* Gulf vs Midwest
 * Fst 0.5 LD-pruned file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMidwest_hiFst_SNPs_firstdraft_Fst0.5.LD.txt`
 * Fst 1.0 LD-pruned file:
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMidwest_hiFst_SNPs_firstdraft_Fst1.0.LD.txt`

```
module load python/3.7-anaconda-2019.10
source activate plink_1_env

START_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/all_samp_tpeds/

cd $START_DIR

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft

SAMP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/natv2_samps.plink.txt

## Atlantic vs Gulf
# Fst_0.5

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.plink.txt

OUT_TXT_1=AvG_Fst0.5_LD_1
OUT_TXT_2=AvG_Fst0.5_LD_2
OUT_TXT_3=AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

# Fst_1.0

cd $START_DIR

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.plink.txt

OUT_TXT_1=AvG_Fst1.0_LD_1
OUT_TXT_2=AvG_Fst1.0_LD_2
OUT_TXT_3=AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

## Atlantic vs Midwest
# Fst_0.5

cd $START_DIR

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMW_hiFst_SNPs_firstdraft_Fst0.5.plink.txt

OUT_TXT_1=AvM_Fst0.5_LD_1
OUT_TXT_2=AvM_Fst0.5_LD_2
OUT_TXT_3=AtlanticVsMidwest_hiFst_SNPs_firstdraft_Fst0.5.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

# Fst_1.0

cd $START_DIR

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMW_hiFst_SNPs_firstdraft_Fst1.0.plink.txt

OUT_TXT_1=AvM_Fst1.0_LD_1
OUT_TXT_2=AvM_Fst1.0_LD_2
OUT_TXT_3=AtlanticVsMidwest_hiFst_SNPs_firstdraft_Fst1.0.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

## Gulf vs Midwest
# Fst_0.5

cd $START_DIR

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMW_hiFst_SNPs_firstdraft_Fst0.5.plink.txt

OUT_TXT_1=GvM_Fst0.5_LD_1
OUT_TXT_2=GvM_Fst0.5_LD_2
OUT_TXT_3=GulfVsMidwest_hiFst_SNPs_firstdraft_Fst0.5.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

# Fst_1.0

cd $START_DIR

SNP_FILE=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMW_hiFst_SNPs_firstdraft_Fst1.0.plink.txt

OUT_TXT_1=GvM_Fst1.0_LD_1
OUT_TXT_2=GvM_Fst1.0_LD_2
OUT_TXT_3=GulfVsMidwest_hiFst_SNPs_firstdraft_Fst1.0.LD.txt

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $SNP_FILE \
--indep-pairwise 10 1 0.75 --out $OUT_TXT_1

plink --bfile GW_all_samps --keep $SAMP_FILE --extract $OUT_TXT_1'.prune.in' \
--indep-pairwise 1000 100 0.99 --out $OUT_TXT_2

cp $OUT_TXT_2'.prune.in' $OUT_DIR

cd $OUT_DIR

cut -d '_' -f 1 $OUT_TXT_2'.prune.in' > tmp_chroms_in.txt
cut -d '_' -f 2 $OUT_TXT_2'.prune.in' > tmp_pos_in.txt
paste tmp_chroms_in.txt tmp_pos_in.txt > $OUT_TXT_3

```

## Filter by missing data
# Test this interactively - should work, but will take some time
```
module load python/3.7-anaconda-2019.07
source activate /global/homes/g/grabowsp/.conda/envs/gen_bioinformatics

OUT_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

cd $OUT_DIR

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

SAMP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/
ATLANTIC_FILE=ATL_firstdraft_train_filt40.txt
GULF_FILE=GULF_firstdraft_train_filt40.txt
MW_FILE=MW_firstdraft_train_filt40.txt

SNP_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/

MISS_CUT=0.8
SNP_FILE_IN=$SNP_DIR'AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.LD.txt'
OUT_TXT_PRE=AvG_Chr
OUT_TXT_SUB=_Fst0.5_LD
SAMP_FILE_1=$SAMP_DIR$ATLANTIC_FILE
SAMP_FILE_2=$SAMP_DIR$GULF_FILE;

for CHR_NUM in {01..03};
  do
  for CHR_LET in K N;
    do
    VCF_FULL=$VCF_DIR$VCF_PRE$CHR_NUM$CHR_LET$VCF_SUF;
    #
    OUT_TXT_1=$OUT_TXT_PRE$CHR_NUM$CHR_LET$OUT_TXT_SUB'.tmp';
    #
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT_1 --positions $SNP_FILE_IN \
      --max-missing $MISS_CUT --keep $SAMP_FILE_1 --kept-sites;
    #
    OUT_TXT_2=$OUT_TXT_PRE$CHR_NUM$CHR_LET$OUT_TXT_SUB'.filt';
    SNP_FILE_2=$OUT_TXT_1'.kept.sites';
    #
    vcftools --gzvcf $VCF_FULL --out $OUT_TXT_2 --positions $SNP_FILE_2 \
      --max-missing $MISS_CUT --keep $SAMP_FILE_2 --kept-sites;
    done;
  done;

```

