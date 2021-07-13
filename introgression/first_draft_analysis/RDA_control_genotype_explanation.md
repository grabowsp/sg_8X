# Explanation of Process of Generating Control genotypes for RDA

## Background
* hiFst RDA uses SNPs chosed because of hi Fst between training sample sets
* a high number of climate associations when using hiFst SNPs and introgressed 
alleles suggests that the introduced alleles have an important role in
climate adaptation
* However, we need a control set of SNPs to compare to the hiFst SNP results
* The idea is to select a random set of SNPs with the same chromosome- and
MAF distribution in a population as the hiFst SNPs, and compare the
number of climate associations for the hiFst and control SNPs

## Steps
1. Identify SNPs used for hiFst RDA
 * These SNPs should have been filted by MAF and presence of ancestry-
informative allele for the introgressing population
 * Currently, this is performed in the R script used for calling ancestry-
informative SNPs in each sample of a population, running window-based 
analyses, generating the genotypes to use for the RDA, and generating
manhattan-style plots showing the putative introgression patterns in each
sample
2. Get allele frequency of hiFst-RDA SNPs
 * Currenlty, this is included in the SNP info output file when the hiFst RDA
genotype file is generated
3. Extract allele frequencies of ALL SNPs using vcfTools
 * Only for population of interest
  * use --keep flag
 * Only process SNPs with MAF at or above the MAF threshold used for the
hiFst-RDA SNPs
  * --maf flag
 * Currently using --freq flag, which adds extra info about the alleles
  * might want to try --freq2 because it outputs less extra info
4. Generate random, control set of SNPs with same chromosome and MAF 
distribution as the hiFst-RDA snps
 * Using R for this, that outputs the SNP postion info needed for vcftools
5. Generate VCF of population-of-interest and control SNPs
 * Usually need to make a VCF for each chromosome separately and then
concatenate the VCFs
6. Generate genotype matrix for RDA
 * Rows = samples, columns = SNPs
 * Using R for this
 * 2 = HOM_REF, 1 = HET, 0 = HOM_ALT 


