#### Author: Rhondene Wint,PhD (rwint@ucmerced.edu)

Purpose: Computes relative synonymous codon usage (RSCU) and absolute counts of the 59 synonymous codons over the entire set (aggregate) of coding sequences(CDS). Implemented  according to  Sharp and Li, 1986  PMCID: PMC340524
- Input: a fasta file of one or more coding sequences (CDS)
- Output: a comma-separated table (.csv) file of the 59 RSCU values
******************************************************************************************************
## How to Use :
1. Ensure that python3 (version 3.5 or higher) is installed. 
	Recommended to install python via anaconda https://docs.anaconda.com/anaconda/install/index.html
2. Download and unzipped the Compute_RSCU_tw folder from github repo into your working folder that contains the input fasta file of CDS.
3. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. Run the command: `python ./Compute_RSCU_tw --help`  for how to add arguments.