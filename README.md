# BuildmSFS
Build the multidimensional SFS from pyrad output files

Introduction:
These scripts are a modification of Jordan Satler's (https://github.com/jordansatler) downsampling script that works for three to ten populations. There are two versions: 

Version 1. AFS_FSC_total_MultPops_Revision.py

This script will build a multidimensional allele frequency spectrum from an input 
matrix of SNPs. For use for as few as three and as many as ten populations. SNP matrix is 
filtered to all bi-allelic SNPs, that equal or surpass the population 
thresholds. The script will sample a single SNP per locus, and if a threshold that requires
subsampling is used, the user can replicate the observed AFS N times 
due to the subsampling of alleles per SNP. Output is an observed allele 
frequency spectrum for use in fastsimcoal2.

Version 2. AFS_FSC_total_MultPops_Revision_linked.py

This script will build a multidimensional allele frequency spectrum from an input 
matrix of SNPs. For use for as few as three and as many as ten populations. SNP matrix is 
filtered to all bi-allelic SNPs, that equal or surpass the population 
thresholds. The script will sample all SNPs per locus, and if a threshold that requires
subsampling is used, the user can replicate the observed AFS N times 
due to the subsampling of alleles per SNP. Output is an observed allele 
frequency spectrum for use in fastsimcoal2. This creates a SFS that uses linked SNPs, which is not recommended
for use in model selection. Instead, if performing model selection, use Version 1.

Input: 
As input, both scripts require two things: 
1.  The user must provide a SNP file. This file should be created from pyRAD SNP output using
    Jordan Satler's script 'SNPtoAFSready.py' (https://github.com/jordansatler).
2.  The user must provide a TRAITS file. This file should be a tab-delimited file. The
    first line should be the headers traits\t species. In the first column, list the names
    of the alleles sampled (as they are listed in your SFS). In the second column, 
    list the population each individual belongs to.  Note that the traits file needs a line for each allele. 
    If you have an individual in your dataset named 'MLS_1', the SNPtoAFSready script will name each allele 
    'MLS_1_a' and 'MLS_1_b'. You need a line for both of these names in your TRAITS file.

Usage: 

The required arguments are as follows, in the specified order

script,Traits,file,Threshold,locus_file,nreps = argv #arguments the user must supply

file is the input descripted in (1) above.

The locus_file is the .loci output from pyRAD.

The Threshold is the user-specified percent of individuals that must be sequenced at a locus for that locus to be retained.

nreps is the number of replicated downsampled SFS the user wants to create.

example usage:

python AFS_FSC_total_MultPops_Done.py traits.txt SNP_infile.txt 50 Monomorphics.txt/species.loci 10 

Output: 
The script will create #nreps downsampled SFS ready for use in FSC2.

If you have any questions or comments, please email me at megansmth67@gmail.com.

Please cite: Smith ML, Ruffley MR, Espíndola AE, Tank DC, Sullivan J, Carstens BC. Demographic model selection using random forests and the site frequency spectrum. Molecular Ecology.

Frequently asked questions: 

1. I'm getting an error that one of my keys is not found (e.g. KeyError: ‘11579_a’). 

First, check that you have a line for each allele in your traits file. 

Second, make sure the line endings are in Unicode (UTF-8). Some endings may cause an error.
