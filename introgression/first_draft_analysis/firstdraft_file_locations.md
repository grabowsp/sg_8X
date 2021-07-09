# Locations of files used for updated introgression analysis for first draft

## Training sample sets
* Midwest
  * 40 samples
    * HAIB: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/MW_firstdraft_train_filt40.txt`
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MW_firstdraft_train_filt40.txt`
  * 60 sample expanded set - includes 40 samples plus more for bootstrapping
    * HAIB: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/MW_firstdraft_train_all60.txt`
    * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/MW_firstdraft_train_all60.txt`
* Atlantic
  * 40 samples: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/ATL_firstdraft_train_filt40.txt`
  * 60 sample expanded set: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/ATL_firstdraft_train_all60.txt`
* Gulf
  * 40 samples: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/GULF_firstdraft_train_filt40.txt`
  * 60 sample expanded set: `/home/f2p1/work/grabowsk/data/switchgrass/introgression_firstdraft/GULF_firstdraft_train_all60.txt`

## High Fst SNP files
* Directory with raw Fst output files
 * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft`
  * File names end with `...Fst_CHROMOSOME_firstdraft.weir.fst`
### Non-pruned hiFst files:
 * Midwest-vs-Atlantic
  * 167,400 with Fst >= 0.5
   * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt`
  * 8,435 with Fst = 1
   * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt`
 * Midwest-vs-Gulf
  * 149,467 with Fst >= 0.5
   * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMW_hiFst_SNPs_firstdraft_Fst0.5.txt`
  * 3,848 with Fst = 1
   * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/GulfVsMW_hiFst_SNPs_firstdraft_Fst1.0.txt`
 * Atlantic-vs-Gulf
  * 65,318 with Fst >= 0.5
   * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst0.5.txt`
  * 569 with Fst = 1
   * NERSC: `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/AtlanticVsGulf_hiFst_SNPs_firstdraft_Fst1.0.txt`
### LD-pruned and missing-data filtered hiFst SNPs
* NERSC Directory:
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/`
* HAIB Directory:
 * `/home/f1p1/tmp/switchgrass_8X/`
* Atlantic vs Gulf
 * Fst 0.5 LD-pruned file:
  * `AtlanticvsGulf_hiFst_SNPs_firstdraft_Fst0.5.LD.MISS.snps.txt`
 * Fst 1.0 LD-pruned file:
  * `AtlanticvsGulf_hiFst_SNPs_firstdraft_Fst1.0.LD.MISS.snps.txt`
* Atlantic vs Midwest
 * Fst 0.5 LD-pruned file:
  * `AtlanticvsMidwest_hiFst_SNPs_firstdraft_Fst0.5.LD.MISS.snps.txt`
 * Fst 1.0 LD-pruned file:
  * `AtlanticvsMidwest_hiFst_SNPs_firstdraft_Fst1.0.LD.MISS.snps.txt`
* Gulf vs Midwest
 * Fst 0.5 LD-pruned file:
  * `GulfvsMidwest_hiFst_SNPs_firstdraft_Fst0.5.LD.MISS.snps.txt`
 * Fst 1.0 LD-pruned file:
  * `GulfvsMidwest_hiFst_SNPs_firstdraft_Fst1.0.LD.MISS.snps.txt`

## High Fst VCFs
* HAIB Directory:
 * `/home/f1p1/tmp/switchgrass_8X/firstdraft_vcfs/`
* NERSC Directory:
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/vcf_files/`
* Midwest vs Atlantic
 * `AtlanticvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`
  * 90,151 SNPs with Fst > 0.5
 * `AtlanticvsMidwest.hiFst.firstdraft.Fst1.0.disomic.CDS.vcf.gz`
  * 1,714 SNPs with Fst = 1
* Midwest vs Gulf
 * `GulfvsMidwest.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`
  * 82,345 SNPs with Fst > 0.5
 * `GulfvsMidwest.hiFst.firstdraft.Fst1.0.disomic.CDS.vcf.gz`
  * 740 SNPs with Fst = 1
* Atlantic vs Gulf
 * `AtlanticvsGulf.hiFst.firstdraft.Fst0.5.disomic.CDS.vcf.gz`
  * 38,762 SNPs with Fst > 0.5
 * `AtlanticvsGulf.hiFst.firstdraft.Fst1.0.disomic.CDS.vcf.gz`
  * 104 SNPs with Fst = 1

## Allele State Files
* HAIB Directory
 * `/home/f1p1/tmp/switchgrass_8X/introgression_files/`
* NERSC Directory
 * `/global/cscratch1/sd/grabowsp/sg_8X_scratch/introgression_files/firstdraft/`
* Midwest vs Atlantic
 * Fst > 0.5
  * `MW_v_ATL_Fst0.5_allele_states.txt`
  * `MW_v_ATL_Fst0.5_keep_pos.txt`
  * `MW_v_ATL_Fst0.5_ref_freq.txt`
 * Fst = 1.0
  * `MW_v_ATL_Fst1.0_allele_states.txt`
  * `MW_v_ATL_Fst1.0_keep_pos.txt`
  * `MW_v_ATL_Fst1.0_ref_freq.txt`
* Midwest vs Gulf
 * Fst > 0.5
  * `MW_v_GULF_Fst0.5_allele_states.txt`
  * `MW_v_GULF_Fst0.5_keep_pos.txt`
  * `MW_v_GULF_Fst0.5_ref_freq.txt`
 * Fst = 1.0
  * `MW_v_GULF_Fst1.0_allele_states.txt`
  * `MW_v_GULF_Fst1.0_keep_pos.txt`
  * `MW_v_GULF_Fst1.0_ref_freq.txt`
* Atlantic vs Gulf
 * Fst > 0.5
  * `ATL_v_GULF_Fst0.5_allele_states.txt`
  * `ATL_v_GULF_Fst0.5_keep_pos.txt`
  * `ATL_v_GULF_Fst0.5_ref_freq.txt`
 * Fst = 1.0
  * `ATL_v_GULF_Fst1.0_allele_states.txt`
  * `ATL_v_GULF_Fst1.0_keep_pos.txt`
  * `ATL_v_GULF_Fst1.0_ref_freq.txt`



