# Codon-Usage-in-Python
Python3  implementation of popular codon usage metrics for manipulating genomic files (.fasta).  I work with hundreds of species in parallel so these scripts to handle batch processing of multiple files - a task that was difficult to accomplish with previously published tools. I will update with more code once I finish the manuscript that I am using these algorithms i.

# Terms
## Codon Usage Bias
The unequal usage synonymous codons within a gene or genome from uniform distribution. 
## Effective Number of Codons (ENC or Nc):
Measures the deviation of synonymous codon usage from unformity. Ranges from 20 (extreme bias, one codon per amino acid) to 61 (no bias, equal usage of synonymous codons for each amino acids)

## Relative Synonymous Codon Usage
The RSCU of a codon is computed as its observed frequency  divided by its expected frequency within a gene or whole transcriptome under the assumption of equal synonymous codon usage. RSCU greater that 1 means that the codon is used more than expected by random chance. [Sharp & Li 1987]. Codons with high RSCU in highly expressed genes are referred to as "optimal codons". For many species the optimal codons are also recognised by the abundant tRNAs [Ikemura 2001]. 

