# Info about extracting and location of header info for the full VCFs

## Location of files
### Full VCF with 1070 Libraries
* Sample line
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tot_VCF_samp_ID.txt`
* Header (without sample line)
  * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tot_VCF_header.txt`


## Sample name line of header for SNP VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files

gunzip -kc Pvirgatum_1070g_Chr01K.snp.sort.norepeats.vcf.gz | head -648 \
| tail -1 > tot_VCF_samp_ID.txt

```

## Header (without sample line) for SNP VCFs
```
gunzip -kc Pvirgatum_1070g_Chr01K.snp.sort.norepeats.vcf.gz | head -647 \
> tot_VCF_header.txt
```

