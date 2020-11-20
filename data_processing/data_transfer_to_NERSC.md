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

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr07N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr08K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr08N.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr09K.5alleles.snp.sort.norepeats.vcf.gz .

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr09N.5alleles.snp.sort.norepeats.vcf.gz .

# run this to make sure everything finished transferring
rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_Chr0*.5alleles.snp.sort.norepeats.vcf.gz .

```

## Check that files are complete
### on HA
```
cd /home/f2p1/work/grabowsk/data/switchgrass/check_8X_tet_files

bash

DATA_DIR=/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/

FILE_PRE=Pvirgatum_1070g_
FILE_SUF=.5alleles.snp.sort.norepeats.vcf.gz

for CHR_TYPE in K N;
  do
  for CHR_NUM in {01..09};
    do
    TEST_CHR=Chr$CHR_NUM$CHR_TYPE;
    echo $TEST_CHR;
    TOT_FILE=$DATA_DIR$FILE_PRE$TEST_CHR$FILE_SUF;
    md5sum $TOT_FILE > $TEST_CHR'.tet.md5';
    done
  done
```
### Transfer files to Cori
* on Cori
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs

scp grabowsk@pants.hagsc.org:/home/f2p1/work/grabowsk/data/switchgrass/check_8X_tet_files/*md5 .

```
### Generate checksum files on Cori
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs

FILE_PRE=Pvirgatum_1070g_
FILE_SUF=.5alleles.snp.sort.norepeats.vcf.gz

for CHR_TYPE in K N;
  do
  for CHR_NUM in {01..09};
    do
    TEST_CHR=Chr$CHR_NUM$CHR_TYPE;
    echo $TEST_CHR;
    TOT_FILE=$FILE_PRE$TEST_CHR$FILE_SUF;
    md5sum $TOT_FILE > $TEST_CHR'.tet.cori.md5';
    done
  done

```
### compare checksum files
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/tetrasomic_vcfs

for i in *md5;
do
cat $i;
done

# visual inspection shows that the md5 files made on the different systems
#  are the same
```

## Transfer MNP file from HA
```
cd /global/cscratch1/sd/grabowsp/sg_8X_scratch/orig_sujan_files/MNP_files

rsync -P grabowsk@pants.hagsc.org:/home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_all_combined.MNP.noRepeats.vcf.gz* .


```

### Checksum
* on HA
```
cd /home/f2p1/work/grabowsk/data/switchgrass/check_8X_tet_files

md5sum /home/smamidi_scratch/Pvirgatum_V5/files_share/Pvirgatum_1070g_variants/Pvirgatum_1070g_all_combined.MNP.noRepeats.vcf.gz > GW.MNP.HA.md5
```

