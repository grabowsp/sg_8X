# Look at general patterns of LD in the MW8X

## Steps
* Extract SNPs using VCFtools
 * all Chr08N CDS snps; MAF > 0.1; Prune within 100bp of each other
* Check how many SNPs
* Calculate LD in MW8X

## Extract SNPs in plink


--thin 100
--maf 0.1
--keep MW8X
--kept-sites
```

VCF_DIR=/global/cscratch1/sd/grabowsp/sg_8X_scratch/all_samp_disomic_vcfs/
VCF_PRE=Chr
VCF_SUF=.disomic.CDS.allsamps.vcf.gz

VCF_IN=$VCF_DIR'Chr08N'$VCF_SUF



```


