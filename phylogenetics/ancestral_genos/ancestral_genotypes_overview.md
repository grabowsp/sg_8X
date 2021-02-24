# Overview of adding ancestral genotypes to VCF

## Overview
* I asked John to estimate ancestral genotypes for the CDS SNPs we're using for the 4X/8X analysis
* I pulled out the SNP positions
* John used genespace to estimate ancestral genotypes
* I generated VCF of ancestral genotypes to merge with other VCFs

## Generating SNP positions
* I used VCFtools to extract the SNP positions from the all_samps VCFs
* Script
  * `~/sg_8X/data_processing/get_SNP_positions_allsamps.md`
* Location of SNP position file:
  * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/SG_8X_snp.positions.txt.gz`
  * HAGSC: `/home/f2p1/work/grabowsk/data/switchgrass/SG_8X_snp.positions.txt.gz`

## Estimating Ancestral Genotypes
* Note from John
```
Hi Paul,
I had hoped to re-calculate these data using the new genespace, but it turns out that isn't quite ready yet and needs a couple tuneups. I won't have the new genespace ready for another week or two.
In lieu of the new pipe, I ran the old one, which should be OK, since I split switchgrass into its subgenomes and ran each independently. I copied the ancestral and reference allele calls, as new columns in your 'snp.positions.txt.gz' file and uploaded as SG_8X_snp.positions_anc.txt.gz. I was able to make ancestral allele calls for 10874224 / 16922099 (64%) of sites. Just those sites with ancestral state calls are included. Also included the AP13 allele in there to QC it... make sure that matches your REF call in the VCF.
JL
```
* Location of file:
  * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/allsamps_tet_vcfs/SG_8X_snp.positions_anc.txt.gz`
  * HAGSC: `/home/f2p1/work/grabowsk/data/switchgrass/SG_8X_snp.positions_anc.txt.gz`



