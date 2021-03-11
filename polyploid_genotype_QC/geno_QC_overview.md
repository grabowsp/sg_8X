# Steps to demonstrate benefit of using tetrasomic genotypes for 8X

## Heterozygosity
* UPDATE: I don't think HET is the best metric because the distribution of
readcounts that affects RAW vs ADJ HET is also what is used for calling
ploidy in nQuire - so they are not independent measures
* Calculate heterozygosity for all samples using diploid genotypes
* Calculate heterozygosity for all samples using tetraploid genotypes
* Look at changes in het values in 4X an 8X samples using the different genotypes

## Trees
* Calculate genetic distance using different genotypes - use same SNPs for each
  * Calc euclidean distance using diploid genotypes for all samples
  * Calc euclidean distance using tetrasomic genotypes for all samples
  * Calc euclidean distance using ploidy-appropriate genotypes
* Generate NJ trees using each distance matrix
* Calculate reproducibility of tres for each method

## Pairwise distances
* Divide samples into groups
  * try using both K=3 and further subdivisions
* Calculate pairwise differences within and between groups using disomic,
tetrasomic, and polyploid genotypes
* Compare differences by ploidy

## PCA
* Generate disomic, tetrasomic, and polyploid genlight objects
  * subselect 100k SNPs
  * repeat 100 times
* Run PCA
* Standardize PC1 and PC2
  * choose "standard" samples that are used for deciding orientation of PCs
  * Standardize PCs to be on scale from 0 to 1
* Get PC1 and PC2 value for each sample
* Calculate variance for each sample

