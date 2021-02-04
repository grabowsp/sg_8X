# Steps for extracting the SNP positions from the allsamps VCF

* use vcftools
* use --kept-sites flag

# First test

```
VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs

VCF_FILE=Chr01K.tetrasomic.CDS.allsamps.vcf.gz

module load python/3.7-anaconda-2019.10
source activate gen_bioinformatics

cd $VCF_DIR

vcftools --gzvcf $VCF_FILE --kept-sites --out Chr01K

tail -n +2 Chr01K.kept.sites > Chr01K.snp.positions.txt

VCF_FILE=Chr01N.tetrasomic.CDS.allsamps.vcf.gz
vcftools --gzvcf $VCF_FILE --kept-sites -c \
| tail -n +2 > Chr01N.snp.positions.txt

for CHR_NUM in {02..09};
  do
  for CHR_LET in K N;
    do
    echo Chr$CHR_NUM$CHR_LET;
    VCF_FILE=Chr$CHR_NUM$CHR_LET'.tetrasomic.CDS.allsamps.vcf.gz';
    vcftools --gzvcf $VCF_FILE --kept-sites -c \
    | tail -n +2 > Chr$CHR_NUM$CHR_LET'.snp.positions.txt';
    done;
  done;

cat *snp.positions.txt > SG_8X_snp.positions.txt
gzip SG_8X_snp.positions.txt

```
