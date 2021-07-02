## Location of Files
* hiFst VCFs
 * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMW_allsampshiFst.v3.disomic.CDS.vcf.gz`
 * `/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/GulfVsMW_allsampshiFst.v3.disomic.CDS.vcf.gz`

## Sample sets
* all 706 samples
 * `/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt`
* Gulf
 * `/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/GULF_all_names.txt`
* Atlantic
 * `/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/ATL_all_names.txt`
* MW_all
 * `/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_all_names.txt`
* MW_introgression_samps
 * `/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_introgressed_names.txt`
* MW_8X_introgression_samps
 * `/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW8X_introgressed_names.txt`
* MW_Low_introgression_samps
 * `/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_no_introgress_names.txt`

## Test with AvM hiFst SNPs
```
bash
source activate bioinformatics_env

OUT_DIR=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/

cd $OUT_DIR

VCF_IN=/home/f2p1/work/grabowsk/data/switchgrass/introgression_v3/AtlanticVsMW_allsampshiFst.v3.disomic.CDS.vcf.gz

MAF_CUT=0.05

vcftools --gzvcf $VCF_IN --chr Chr08N --out AvM_Chr08N_all \
--geno-r2

# 

SAMPS_IN=/home/grabowsky/tools/workflows/sg_8X/sample_sets/samp_set_files/nat_filt_v2_names.txt

OUT_PRE=AvM_Chr08N_natv2

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--maf $MAF_CUT --geno-r2

#

SAMPS_IN=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/GULF_all_names.txt

OUT_PRE=AvM_Chr08N_GULF

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--maf $MAF_CUT --geno-r2

#

SAMPS_IN=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/ATL_all_names.txt

OUT_PRE=AvM_Chr08N_ATL

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--maf $MAF_CUT --geno-r2

# 

SAMPS_IN=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_all_names.txt

OUT_PRE=AvM_Chr08N_MWall

MAF_CUT=0.05

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--maf $MAF_CUT --geno-r2

#

SAMPS_IN=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_introgressed_names.txt

OUT_PRE=AvM_Chr08N_MWintrogressed

MAF_CUT=0.05

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--maf $MAF_CUT --geno-r2

#

SAMPS_IN=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW8X_introgressed_names.txt

OUT_PRE=AvM_Chr08N_MW8X_introgressed

MAF_CUT=0.05

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--geno-r2

#

SAMPS_IN=/home/f2p1/work/grabowsk/data/switchgrass/PlantGroup_presentation_July2021/MW_no_introgress_names.txt

OUT_PRE=AvM_Chr08N_MW_noIntrogression

vcftools --gzvcf $VCF_IN --chr Chr08N --out $OUT_PRE --keep $SAMPS_IN \
--maf $MAF_CUT --geno-r2







```
