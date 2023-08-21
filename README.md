# Python tools for Codon Usage Bias Analysis
Python3 command-line programs for calculating popular gene-specific or genome-wide codon usage frequencies, e.g. RSCU, from sequence files (.fasta).  I work with hundreds of species in parallel so these scripts to handle batch processing of multiple files - a task that was difficult to accomplish with previously published tools. 

## Tools:
- <b>Compute_RSCU_gene </b>:  computes the RSCU for each individual  transcript
- <b>Compute_RSCU_tw </b>:  computes the RSCU and absolute codon counts over the entire CDS ('transciptome-wde') in fasta file.
- <b>CodonUsage_per_1000 </b>:  Computes codon usage per 1000 of the whole transcriptome.
- <b>fasta2csv</b>: converts fasta file to two-column csv table (Header | Sequence); 
- <b> aa_usage </b>: computes the Amino acid usage
- <b> fix_fasta.py </b>: corrects the issue of newlines within the same sequence. 

# Terms
## Codon Usage Bias
The unequal usage synonymous codons within a gene or genome from uniform distribution. 
 synonymous codon usage from unformity. Ranges from 20 (extreme bias, one codon per amino acid) to 61 (no bias, equal usage of synonymous codons for each amino acids)

## Relative Synonymous Codon Usage
The RSCU of a codon is computed as its observed frequency  divided by its expected frequency within a gene or whole transcriptome under the assumption of equal synonymous codon usage. RSCU greater that 1 means that the codon is used more than expected by random chance. [Sharp & Li 1987]. Codons with high RSCU in highly expressed genes are referred to as "optimal codons". For many species the optimal codons are also recognised by the abundant tRNAs [Ikemura 2001]. 

