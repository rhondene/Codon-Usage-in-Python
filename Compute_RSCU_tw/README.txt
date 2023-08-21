## Author: Rhondene Wint,PhD (rwint@ucmerced.edu)

Purpose: Computes relative synonymous codon usage (RSCU) of the ALL codons over the entire set of coding sequnces(CDS)/transcriptome-wide:
            implemented  according to  Sharp and Li, 1986  PMCID: PMC340524
Input: a fasta file of coding sequences (CDS)
Output: a comma-separated table (.csv) file of the RSCU values
******************************************************************************************************
## How to Use :
1. Ensure that python3 (version 3.5 or higher) is installed. 
	Recommended to install python via anaconda https://docs.anaconda.com/anaconda/install/index.html
2. Download and unzipped the Compute_RSCU_tw folder from github repo into your working folder that contains the input fasta file of CDS.
3. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. Run the command: 'python ./Compute_RSCU_tw --help'  for how to add arguments.
