#Author: Rhondene Wint (rwint@ucmerced.edu)

- Purpose: Computes relative synonymous codon usage of each 59 degenerate codons per each coding sequence (CDS)
            according to Sharp and Li, 1986
- Input:  FASTA file of N coding sequences (CDS)
- Output: comma-separated table (csv) of the relative synonymous codon usage for each transcript: i.e. a matrix of N transcripts x 59 RSCU values
******************************************************************************************************
How to Use :
1. Ensure that python3 (version 3.8 or higher) and pandas v2.0 or higher are installed. 
	Recommended to install python via anaconda https://docs.anaconda.com/anaconda/install/index.html
2. download the Compute_RSCU_gene.pyz binary from github repo into your working folder containing the input FASTA file.
3. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. Type the following in the terminal, be sure to replace the input and output arguments with your own :
   - `python Compute_RSCU_gene.pyz -CDS example_cds.fasta -out rscu_results `
   - Also run `python  ./Compute_RSCU_gene --help`  for help menu.
