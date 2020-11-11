# Steps for transferring data from HA to NERSC

## Transfer disomic SNPs
* transfer the VCF and index file
```
# on NERSC

cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files

scp grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr01K.snp.sort.norepeats.vcf.gz .

scp grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr01K.snp.sort.norepeats.vcf.gz.csi .

scp grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr01N.snp.sort.norepeats.vcf.gz* .

scp grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr02K.snp.sort.norepeats.vcf.gz* .

scp grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr02N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr03K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr03N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr04K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr04K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr04N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr05K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr05N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr06K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr06N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr07K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr07N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr08K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr08N.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr09K.snp.sort.norepeats.vcf.gz* .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr09N.snp.sort.norepeats.vcf.gz* .
```

## Transfer tetrasomic VCFs
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr01K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr01N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr02K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr02N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr03K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr03N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr04K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr04N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr05K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr05N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr06K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr06N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr07K.5alleles.snp.sort.norepeats.vcf.gz .

```
