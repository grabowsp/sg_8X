# Steps for adding ancestral genotypes to all_samps VCF

## Steps
1. Make file of positions of SNPs with ancestral genotypes
2. Extract REF and ALT alleles for those positions
3. Compare REF and ALT to AP13 genoyptes John included in file
4. Make script calling REF and ALT for ancestral genotypes
5. Make VCF for ancestral (and AP13) genotypes
6. Sort, bgzip, and tabix the VCF files prior to merging
  * see: `~/sg_ha_effort/polyploid_genos/genotype_file_generation/combine_sim8X_geo.md`

## Important Files
* Original all_samps CDS SNPs positions
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/SG_8X_snp.positions.txt.gz`
* Ancestral alleles for CDS SNPs
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/SG_8X_snp.positions_anc.txt.gz`
* Position file for SNPs with ancestral allele calls
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/sg_8x_ancestral_snp_pos.txt`


## Make file with Ancestral genotype positions
* `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/sg_8x_ancestral_snp_pos.txt`
```
#module load python/3.7-anaconda-2019.10
#source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs

gunzip -c SG_8X_snp.positions_anc.txt.gz | \
awk -v FS=',' -v 'OFS'='\t' '/Chr/ {print $1,$2}' > \
sg_8x_ancestral_snp_pos.txt 


```

## Get REF and ALT alleles from VCFs
* 4th column = REF
* 5th column = ALT
```
module load python/3.7-anaconda-2019.10
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs

VCF_FILE=Chr01K.tetrasomic.CDS.allsamps.vcf.gz
SNP_POS_FILE=sg_8x_ancestral_snp_pos.txt
OUT_FILE=Chr01K.ancestral_snp_ref_alt.txt

vcftools --gzvcf $VCF_FILE --positions $SNP_POS_FILE --recode -c | \
awk -v FS='\t' -v OFS='\t' '/^Chr/ {print $1,$2,$3,$4,$5}' > $OUT_FILE

```

## Make chromosome file of ancestral genotypes
```
module load python/3.7-anaconda-2019.10
source activate gen_bioinformatics

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs

gunzip -c SG_8X_snp.positions_anc.txt.gz | \
awk -v FS=',' -v 'OFS'='\t' '/^Chr01K/ {print $1,$2,$4,$5}' > \
Chr01K.sg_8x_ancestral_genos.txt

```

## Compare genotypes in R
```
# module load python/3.7-anaconda-2019.10
# module swap PrgEnv-intel PrgEnv-gnu
# source activate R_tidy

### LOAD PACKAGES ###

library(data.table)

### IMPORT DATA ###

ancestral_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/Chr01K.sg_8x_ancestral_genos.txt'

ancest_info <- fread(ancestral_file, header = F)
colnames(ancest_info) <- c('CHR', 'POS', 'ANCESTRAL', 'AP13')

ra_file <- '/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/Chr01K.ancestral_snp_ref_alt.txt'

ra_info <- fread(ra_file, header = F)
colnames(ra_info) <- c('CHR', 'POS', 'NAME', 'REF', 'ALT')

setdiff(ancest_info$POS, ra_info$POS)
# none

ancest_info <- ancest_info[-which(duplicated(ancest_info$POS))]

if(sum(ancest_info$POS != ra_info$POS) == 0){
  combo_info <- merge(ancest_info, ra_info)
}

sum(combo_info$AP13 != combo_info$REF)
# [1] 0

sum(combo_info$ANCEST != combo_info$AP13)
# [1] 39764

sum(combo_info$ANCEST != combo_info$AP13) / nrow(combo_info)
# [1] 0.06063055

anc_alt_inds <- which(combo_info$ANCEST != combo_info$AP13)

sum(combo_info$ANCEST[anc_alt_inds] != combo_info$ALT[anc_alt_inds])
# 9138

mismatch_ancest_inds <- intersect(which(combo_info$ANCEST != combo_info$ALT), 
  anc_alt_inds)

length(mismatch_ancest_inds)/nrow(combo_info)
# [1] 0.01393326

combo_info_1 <- combo_info[-mismatch_ancest_inds]

# VCF geno = R/A
# tet HOM REF = 4/0:20,0
# tet HOM ALT = 0/4:0,20
# dip HOM REF = 2/0:20,0
# dip HOM ALT = 0/2:0,20

combo_info_1[, ANCESTRAL_TET_GENO := as.character(NA)]
combo_info_1[, ANCESTRAL_DIP_GENO := as.character(NA)]
combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$REF), 
  ANCESTRAL_TET_GENO := '4/0:20,0']
combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$REF),
  ANCESTRAL_DIP_GENO := '2/0:20,0']

combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$ALT),
  ANCESTRAL_TET_GENO := '0/4:20,0']
combo_info_1[which(combo_info_1$ANCESTRAL == combo_info_1$ALT),
  ANCESTRAL_DIP_GENO := '0/2:0,20']

combo_info_1[, QUAL := '.']
combo_info_1[, FILTER := '.']
combo_info_1[, INFO := '.']
combo_info_1[, FORMAT := '']


combo_info_1[, c('ANCESTRAL_TET_GENO', 'REF_TET_GENO')]




# NEXT
# remove the mismatch_ancest_inds
# generate VCF-like ancestral and AP13 genotypes
# generate VCF file with those genotypes



```


