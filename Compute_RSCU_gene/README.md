#Author: Rhondene Wint (rwint@ucmerced.edu)

Purpose: Computes relative synonymous codon usage of each 59 degenerate codons per each coding sequence (CDS)
            according to Sharp and Li, 1986
Input:  FASTA file of N coding sequences (CDS)
Output: comma-separated table (csv) of the relative synonymous codon usage for each transcript: i.e. a matrix of N transcripts x 59 RSCU values
******************************************************************************************************
How to Use :
1. Ensure that python3 (version 3.5 or higher) is installed. 
	Recommended to install python via anaconda https://docs.anaconda.com/anaconda/install/index.html
2. download and unzipped the Compute_RSCU_gene folder from github repo into your working folder containing the input FASTA file.
3. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. Run the command: `python ./Compute_RSCU_gene --help'  for how to add arguments.`