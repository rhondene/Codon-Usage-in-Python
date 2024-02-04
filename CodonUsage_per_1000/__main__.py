#!/usr/bin/python 

# Author: Rhondene Wint
import argparse
import pandas as pd
import numpy as np
import fix_fasta
from pandas.errors import SettingWithCopyWarning
import warnings
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)


#codon_aa_table
codon_aa = {
	"UUU":"Phe", "UUC":"Phe",		  
	"UCU":"Ser4", "UCC":"Ser4", "UCA":"Ser4", "UCG":"Ser4",
	"AGU":"Ser2", "AGC":"Ser2",
	"CUU":"Leu4", "CUC":"Leu4", "CUA":"Leu4", "CUG":"Leu4",
	"UUA":"Leu2", "UUG":"Leu2",
	
	"UAU":"Tyr", "UAC":"Tyr", "UAA":"STOP", "UAG":"STOP",
	"UGU":"Cys", "UGC":"Cys", "UGA":"STOP", "UGG":"Trp",
	"CGU":"Arg4", "CGC":"Arg4", "CGA":"Arg4", "CGG":"Arg4",
	"AGA":"Arg2", "AGG":"Arg2",
	"CCU":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
	"CAU":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
	
	"AUU":"Ile", "AUC":"Ile", "AUA":"Ile", "AUG":"Met",
	"ACU":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
	"AAU":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
   
	"GUU":"Val", "GUC":"Val", "GUA":"Val", "GUG":"Val",
	"GCU":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
	"GAU":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
	"GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"}
	

##compute codon count
def get_cod_freq(headers:list, seqs:list)->pd.DataFrame:
	""" Computes the codon usage per 1000 over the entire set of CDS
		headers: list of headers
		seqs: list of CDS
		Returns a 59-dim dataframe of total absolute codon frequencies
	"""
	###build dictionary to accumulate codon counts; omit  stop codons
	codon_count=dict() 
	STOP_CODONS = ['UAA', 'UAG', 'UGA']
	codon_count = {codon: 0 for codon in codon_aa if codon not in STOP_CODONS }

	for i,cds in enumerate(seqs):
		cds = cds.upper().replace('T','U')
		if len(cds)%3 !=0:
			ID = headers[i].split(' ')[0]
			print(f"WARNING: Skipping {ID} Length of CDS is not a multiple of 3.")
			continue
		
		##count codons in cds
		for i in range(0, len(cds), 3):
			codon = cds[i:i+3]
			if codon in codon_count:
				codon_count[codon] += 1
	
	df_codcount=pd.DataFrame(list(codon_count.items()) )
	df_codcount.columns=['Codon', 'Frequency']
	df_codcount['Amino_Acid'] = [codon_aa[codon] for codon in df_codcount['Codon'].values] ## add amino acid column
	per_1000 = df_codcount['Frequency'].sum()/1000
	df_codcount['Frequency'] = df_codcount['Frequency']/per_1000
	
	return df_codcount[['Codon','Amino_Acid','Frequency']]	#reorder the columns
	
if __name__=='__main__':

	about = 'Computes codon usage per 1000 of a whole transcriptome. Rules: \n Skips CDS if the length of CDS is not a multiple of 3 \n Skips Codons with unknown bases.  \nWritten by Rhondene Wint, rwint@ucmerced.edu.'
	epi_note = 'To contact the author about problems or errors,make a pull request at https://github.com/rhondene/Codon-Usage-in-Python'
	parser = argparse.ArgumentParser(description=about,epilog=epi_note)
	parser.add_argument('-CDS', help='Path to fasta file with species coding sequences', type=str, required=True, metavar='')
	parser.add_argument('-out', help='Path of destination folder for output csv file', type=str, default='./file.cu.csv', metavar='')	

	args=parser.parse_args()
	
	headers,seqs=fix_fasta.fix_fasta(args.CDS)##formats fasta into csv of sequences
	CU = get_cod_freq(headers,seqs)	 ##computes codon usage per 1000

	#save the file
	CU.to_csv('{}.cu1000.csv'.format(args.out), index=False)


	
	
	
	
	
	