## About CodonCount


Purpose: Command-line tool that computes the length normalized codon frequency of each 61 sense codons of a coding sequence (CDS), and returns a comma-separated table.
            
	    Relative freq. of codon_i=  (frequency of codon_i)/(total number of codons in the CDSj)
******************************************************************************************************
How to Use :
1. Ensure that python3 (version 3.8 or higher) and pandas v2.0 or higher are installed.  
	Recommended to install python3 via anaconda https://docs.anaconda.com/anaconda/install/index.html 
2. Download the *CodonCount.pyz* file in this github repo into your working folder with the input fasta file(s). 
3. Open a terminal window (bash, gitbash, powershell, etc) in the same working folder.
4. To run the programn, tpe the command below in the terminal shell (be sure to replace arguments with the actual name the input and output files):
	```console 
	python3 CodonCount.pyz -CDS example.fasta -out example_output
	```
	
