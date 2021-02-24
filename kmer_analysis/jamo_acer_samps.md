# Process for finding the libraries for Acer
* some on NERSC, some on HA, I think

```
cd /global/homes/g/grabowsp/jamo_reports

for LIBNAME in JBDT IJBQ IEYM JBJF IFBD IRYU IELR INFN ITAQ JBFD IICN ITAW IIBP JBEG IIJD IGQW IIJX JBFW IRUX JBEE JBEH ITAT IEYI IEMQ IHZD IEMF IEMM IENN IRUR IENK;
  do
  jamo info library $LIBNAME >> /global/cscratch1/sd/grabowsp/sg_8X_scratch/samp_filtering/acer_pilot_jamo_out.txt;
  done


```
