# Python tools for Codon Usage Bias Analysis
- Python3 command-line programs for calculating popular single-gene or genome-wide codon usage frequencies, e.g. RSCU, from sequence files (.fasta).  I worked with hundreds of species in parallel so these scripts to handle batch processing of multiple files and outputs a CSV format table that is easier to parse - a task that was difficult to accomplish with previously published tools. 
- These tools were validated against the original CodonW software by Peden, 1995

# Tools and Their Usage:
All tools require that  python3 version 3.8 or higher is installed. Recommended to install python3 via anaconda https://docs.anaconda.com/anaconda/install/index.html 
### Compute_RSCU_gene :  
- Computes relative synonymous codon usage of each 59 degenerate codons per each coding sequence (CDS) according to Sharp and Li, 1986 PMCID: PMC340524
- Input:  FASTA file of N coding sequences (CDS)
- Output: comma-separated table (csv) of the relative synonymous codon usage for each transcript: i.e. a matrix of N transcripts x 59 RSCU values
  ******************************************************************************************************
How to Use :

<li> Download the <i>Compute_RSCU_gene.pyz</i> binary from the Compute_RSCU_gene github repo into your project folder containing the input FASTA file.</li>
<li>Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.</li> 
<li>Type the following in the terminal, be sure to replace the input and output arguments with your own :</li>
     
     ```
       python Compute_RSCU_gene.pyz -CDS example_cds.fasta -out rscu_results
   - Also run ```python Compute_RSCU_gene.pyz --help```  for help menu.

### Compute_RSCU_tw :  
-  Computes relative synonymous codon usage (RSCU) and absolute counts of the 59 synonymous codons over the entire set (aggregate) of coding sequences('transcriptome-wide'). Implemented  according to  Sharp and Li, 1986  PMCID: PMC340524
- Input: single or multifasta file of coding sequences (CDS)
- Output: a comma-separated table (.csv) file of the 59 RSCU values
  ******************************************************************************************************
How to Use :
1. Download Compute_RSCU_tw.pyz binary from Compute_RSCU_tw repo into your working folder that contains the input fasta file of CDS.
2. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
3. To run the programn, type the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files): ```python Compute_RSCU_tw.pyz -CDS example.fasta -out results```  

### CodonCount: 
- Command-line tool that computes the length normalized codon usage of each 61 sense codons of a coding sequence (CDS), and returns a comma-separated table.
            
	    Relative freq. of codon_i=  (frequency of codon_i)/(total number of codons in the CDSj)
******************************************************************************************************
How to Use :
2. Download the *CodonCount.pyz* file in CodonCount github repo into your working folder with the input fasta file(s). 
3. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. To run the programn, type the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files):
	```console
	python CodonCount.pyz -CDS example.fasta -out example_output
 Also run ```python CodonCount.pyz --help```  for help menu.

 ### CodonUsage_per_1000:  
- Computes codon usage per 1000 of the whole transcriptome.
2. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. To run the programn, type the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files):
	```console
	python CodonUsage_per_1000.pyz -CDS all_CDS.fasta -out  results_cu
 Also run ```python CodonUsage_per_1000.pyz --help```  for help menu.
### fasta2csv : 
- Converts fasta file to two-column csv table (Header | Sequence); 
### aa_usage :
- Computes the Amino acid usage
### fix_fasta.py: 
- Corrects the issue of newlines within the same sequence. 

# Terms
## Codon Usage Bias
The unequal usage synonymous codons within a gene or genome i.e. the deviation of synonymous codons from a uniform distribution. 

## Relative Synonymous Codon Usage
<li> The RSCU of a codon is computed as its observed frequency  divided by its expected frequency within a gene or whole transcriptome under the null hypothesis of equal synonymous codon usage. </li>
<li> RSCU greater that 1 means that the codon is used more than expected by random chance. [Sharp & Li 1987]. </li>
<li>Codons with high RSCU in highly expressed genes are referred to as "optimal codons". For many species the optimal codons are selectively recognised by the abundant tRNAs, which is often taken as an indication selection pressures shaping codon usage patterns [Ikemura 1983; Wint et al 2022]. </li>

