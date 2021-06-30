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
* Non-pruned hiFst files:
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


