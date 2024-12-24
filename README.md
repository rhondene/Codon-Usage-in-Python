# Python tools for Codon Usage Bias Analysis  - 
- This repository hosts different Python3 command-line programs for calculating popular codon usage and amino acid usage frequency statistics from FASTA sequence files (.fasta).
- To use these tools , all that needs to be done is to install python3 and then download the executable binary file <b>(.pyz)</b> from the tool's folder in this repo.	
- <b> Motivation </b>: I worked with hundreds of genomes so I wrote these scripts to handle batch processing of multiple genomes/ input files and outputs a CSV formatted table that is easier to parse and amenable to statistical analysis like PCA - a task that I found tedious  because  previously published tools would output the conventional wide-form codon usage table that needed extra processing.
  	- <b> Genomics Publications that use these tools <b>:
  		<ul>
    <li><a href="https://www.mdpi.com/2073-4425/15/12/1602" target="_blank">Transfer RNA Levels Are Tuned to Support Differentiation During Drosophila Neurogenesis</a></li>
	</ul>	
 <ul>
    <li><a href="https://academic.oup.com/mbe/article/39/2/msab372/6513383" target="_blank">Kingdom-Wide Analysis of Fungal Protein-Coding and tRNA Genes Reveals Conserved Patterns of Adaptive Evolution</a></li>
	</ul>
   
- These Codon Usage tools were validated against the original CodonW software by Peden, 1995
  
   <b> <font color='d#f59e42'> ***[Interactive Jupyter Notebook version Coming Soon in Summer 2024!]*** </font></b>
   
# Software Setup:
- All tools require python3 version 3.8 or higher is installed and pandas version 2.0 or higher.
- Recommended to install python3 via anaconda https://docs.anaconda.com/anaconda/install/index.html 

- <font color='green'> See </font> the ```test_data``` folder for examples of the outputs of each tool on the same input fasta file ('NB_CDS.fasta')
# ------How to Use Each Tool ------

### Compute_RSCU_gene :  
- Computes relative synonymous codon usage of each 59 degenerate codons per each coding sequence (CDS) according to Sharp and Li, 1986 PMCID: PMC340524
- `Input`:  FASTA file of N coding sequences (CDS)
- `Output`: comma-separated table (csv) of the relative synonymous codon usage for each transcript: i.e. a matrix of N transcripts x 59 RSCU values
  ******************************************************************************************************
<b>How to Use</b> :

1. Download the `Compute_RSCU_gene.pyz` binary from the Compute_RSCU_gene github repo into your project folder containing the input FASTA file.
2. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.</li> 
3. Type the following in the terminal, be sure to replace the names of the input and output arguments with your own :</li>

```console
	python Compute_RSCU_gene.pyz -CDS example_cds.fasta -out rscu_results
 ```
    
   - Also run ```python Compute_RSCU_gene.pyz --help```  for help menu.

### Compute_RSCU_tw :  
-  Computes relative synonymous codon usage (RSCU) and absolute counts of the 59 synonymous codons over the entire set (aggregate) of coding sequences('transcriptome-wide'). Implemented  according to  Sharp and Li, 1986  PMCID: PMC340524
- `Input`: single or multifasta file of coding sequences (CDS)
- `Output`: a comma-separated table (.csv) file of the 59 RSCU values
  ******************************************************************************************************
**How to Use** :
1. Download Compute_RSCU_tw.pyz binary from Compute_RSCU_tw repo into your working folder that contains the input fasta file of CDS.
2. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
3. To run the programn, type the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files):

   ```console
   	python Compute_RSCU_tw.pyz -CDS example.fasta -out results
   ```  

### CodonCount: 
Computes the length normalized codon frequency of each 61 sense codons of a coding sequence (CDS), and returns CSV .
            
	    Relative Frequency of Codon_i=  (frequency of codon_i)/(total number of codons in the CDSj)
******************************************************************************************************
**How to Use** :
1.  Download the `CodonCount.pyz` file in CodonCount github repo into your working folder with the input fasta file(s). 
2. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
3.  To run the programn, type the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files):
	```console
	python CodonCount.pyz -CDS example.fasta -out example_output
 	``` 
 Also run ```python CodonCount.pyz --help```  for help menu.

 ### CodonUsage_per_1000:   
 Computes codon usage per 1000 of the whole transcriptome.
1. Download the `CodonUsage_per_1000.pyz` file in CodonUsage_per_1000 github repo into your working folder with the input fasta file(s).
2. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. To run the programn, type the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files):
	```console
	python CodonUsage_per_1000.pyz -CDS all_CDS.fasta -out  results_cu
  	```
 Also run ```python CodonUsage_per_1000.pyz --help```  for help menu.
### fasta2csv : 
- Converts fasta file to two-column csv table (Header | Sequence); 
### aa_usage :
- Computes the Expected and Observed Amino acid usage according to methods outlined in  and https://pubmed.ncbi.nlm.nih.gov/5767777/ the https://qubeshub.org/publications/979/serve/1/3067?el=1&download=1 
- To run, download the script in your project folder and type in the terminal

  ``` python aa_usage.py -CDS YOUR_CDS.fasta -out OUTPUT_NAME```
  
### fix_fasta.py: 
- Corrects the issue of newlines within the same sequence. 

# Glossary Codon Usage Metrics

## Codon Usage Bias
The unequal usage of synonymous codons within a gene or genome i.e. the deviation of synonymous codons from a uniform distribution due to a combination of natural selection, neutral mutational bias and genetic drift. 

## Relative Synonymous Codon Usage
<li> The RSCU of a codon is computed as its observed frequency  divided by its expected frequency within a gene or whole transcriptome under the null hypothesis of equal synonymous codon usage. </li>
<li> RSCU greater that 1 means that the codon is used more than expected by random chance. [Sharp & Li 1987]. </li>
<li>Codons with high RSCU in highly expressed genes are referred to as "optimal codons". For many species the optimal codons are selectively recognised by the abundant tRNAs, which is often taken as an indication selection pressures shaping codon usage patterns [Ikemura 1983; Wint et al 2022]. </li>

## Amino Acid Frequency:
- If a particular amino acid is in some way adaptive, then it should occur more frequently than expected by chance.
- This can easily be tested by calculating the expected frequencies of amino acids and comparing to observed. The codons and observed frequencies of particular amino acids are given in the table.
- The frequencies of DNA bases in nature are 22.0% uracil, 30.3% adenine, 21.7% cytosine, and 26.1% guanine. The expected frequency of a particular codon can then be calculated by multiplying the frequencies of each DNA base comprising the codon. The expected frequency of the amino acid can then be calculated by adding the frequencies of each codon that codes for that amino acid.
- As an example, the RNA codons for tyrosine are UAU and UAC, so the random expectation for its frequency is (0.220)(0.303)(0.220) + (0.220)(0.303)(0.217) = 0.0292. Since 3 of the 64 codons are nonsense or stop codons, this frequency for each amino acid is multiplied by a correction factor of 1.057.
