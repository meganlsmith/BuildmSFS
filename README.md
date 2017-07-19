# BuildmSFS
Build the multidimensional SFS from pyrad output files

This script is a modification of Jordan Satler's (https://github.com/jordansatler) downsampling script that works for three to ten populations.

This script will build a multidimensional allele frequency spectrum from an input 
matrix of SNPs. For use for as few as three and as many as ten populations. SNP matrix is 
filtered to all bi-allelic SNPs, that equal or surpass the population 
thresholds. The script will sample a single SNP per locus, and if a threshold that requires
subsampling is used, the user can replicate the observed AFS N times 
due to the subsampling of alleles per SNP. Output is an observed allele 
frequency spectrum for use in fastsimcoal2.

The version named with an appended 'linked' does the same, except it will not sample only one SNP per locus, and will instead use multiple SNPs per locus, creating a SFS that includes linked SNPs. Such an SFS should not be used for model selection in fastsimcoal2. I recommend using the other script, which does not use linked SNPs, if you're planning to perform model selection.

For use with modified pyRAD output.

Note that the infile ('SNP_infile.txt') should be created from pyRAD SNP output using Jordan Satler's script: 'SNPtoAFSready.py'

Note that the traits file needs a line for each allele. If you have an individual in your dataset named 'MLS_1', the SNPtoAFSready script will name each allele 'MLS_1_a' and 'MLS_1_b'. You need a line for both of these names in your traits.txt file.

                                                     
python AFS_FSC_total_MultPops_Done.py traits.txt SNP_infile.txt 
       Threshold Monomorphics.txt/species.loci Replicate 

"""

Please cite: Smith ML, Ruffley MR, Espíndola AE, Tank DC, Sullivan J, Carstens BC. Demographic model selection using random forests and the site frequency spectrum. Molecular Ecology.

Frequently asked questions: 

1. I'm getting an error that one of my keys is not found (e.g. KeyError: ‘11579_a’). 
First, check that you have a line for each allele in your traits file. 
Second, make sure the line endings are in Unicode (UTF-8). Some endings may cause an error.
