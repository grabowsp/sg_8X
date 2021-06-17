# Set up globus file for Acer to use for accessing data

## Set up directory
* `/global/cfs/cdirs/m342/gsharing/sg_acer`

## Make softlinks for files
* Location of files from Sujan
  * /global/cscratch1/sd/sujan/Pvirgatum_Acer_paul

```
cd /global/cscratch1/sd/sujan/Pvirgatum_Acer_paul
for BASE_FILE in `ls *`;
  do
  ln -s $BASE_FILE /global/cfs/cdirs/m342/gsharing/sg_acer/$BASE_FILE;
  done

```

```
cd /global/cfs/cdirs/m342/gsharing/sg_acer
rm libs_30g

cp /global/cscratch1/sd/sujan/Pvirgatum_Acer_paul/libs_30g .

rm *bz2
rm *bam
rm libs*

scp /global/cscratch1/sd/sujan/Pvirgatum_Acer_paul/* .
```

## Second set of files
### Location of files
* globus directory
  * `/global/cfs/cdirs/m342/gsharing/sg_acer`
  * Sujan's file with prepped fastq files
    * `/global/cscratch1/sd/sujan/Pvirgatum_Acer_paul/`
### Transfer files to glubus directory
```
cd /global/cfs/cdirs/m342/gsharing/sg_acer

scp /global/cscratch1/sd/sujan/Pvirgatum_Acer_paul/*fastq.bz2 .
# waiting to do this until permissions are adjusted on the files
```


