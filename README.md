# Python tools for Codon Usage Bias Analysis
- Python3 command-line programs for calculating popular single-gene or genome-wide codon usage frequencies, e.g. RSCU, from sequence files (.fasta).  I worked with hundreds of species in parallel so these scripts to handle batch processing of multiple files and outputs a CSV format table that is easier to parse - a task that was difficult to accomplish with previously published tools. 
- These tools were validated against the original CodonW software by Peden, 1995

## Tools:
- <b>Compute_RSCU_gene </b>:  computes the RSCU for each individual  Coding Sequence (CDS). The output is a matrix of 59 codons x n genes
- <b>Compute_RSCU_tw </b>:  computes the RSCU and absolute codon counts of a <i> single CDS </i> or over entire set of CDS ('transciptome-wde') in fasta file.
- <b>CodonUsage_per_1000 </b>:  Computes codon usage per 1000 of the whole transcriptome.
- <b>fasta2csv</b>: converts fasta file to two-column csv table (Header | Sequence); 
- <b> aa_usage </b>: computes the Amino acid usage
- <b> fix_fasta.py </b>: corrects the issue of newlines within the same sequence. 

# Terms
## Codon Usage Bias
The unequal usage synonymous codons within a gene or genome i.e. the deviation of synonymous codons from a uniform distribution. 

## Relative Synonymous Codon Usage
- The RSCU of a codon is computed as its observed frequency  divided by its expected frequency within a gene or whole transcriptome under the null hypothesis of equal synonymous codon usage. 
- RSCU greater that 1 means that the codon is used more than expected by random chance. [Sharp & Li 1987].
- Codons with high RSCU in highly expressed genes are referred to as "optimal codons". For many species the optimal codons are selectively recognised by the abundant tRNAs, which is often taken as an indication selection pressures shaping codon usage patterns [Ikemura 1983; Wint et al 2022]. 

