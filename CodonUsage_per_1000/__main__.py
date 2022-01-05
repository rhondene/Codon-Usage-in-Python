#!/usr/bin/python 

# Author: Rhondene Wint
import argparse
import pandas as pd
import numpy as np
import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning) #suppress copy warning

import fix_fasta

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
	
	
### preprocess the fasta file


def preproc(fasta_file):
	"""formats fasta sequences to a list"""
	
	#flybase fasta file has internal newline in the same seqeunce 
	seqs=fix_fasta.fix_fasta(fasta_file)[1] #contains list of sequences
	return seqs
##compute codon count

def get_cod_freq(seqs):
	""" seqs: list of CDS
		Returns a 59-dim dataframe of total absolute codon frequencies
	"""
	
	codon_count=dict() 

	for codon in list(codon_aa.keys()):
		codon_count[codon]=0  ##dictionary to accumulate codon count
	#count the codons in each CDS	
	for cds in seqs:
		cds = cds.upper().replace('T','U')
		codons = []
		##make a list of codons
		for c in range(0,len(cds),3):
			if len(cds)%3 ==0:
				cod=cds[c:c+3]
				if 'N' not in cod:	##ignore N and seqs not multiple of 3
					codons.append(cod)
			else:
				continue

		for c in list(codon_count.keys()):
			codon_count[c]+= codons.count(c)
	
	df_codcount=pd.DataFrame(list(codon_count.items()) )
	df_codcount.columns=['Codon', 'Frequency']
	df_codcount['Amino_Acid'] = [codon_aa[codon] for codon in df_codcount['Codon'].values] ## add amino acid column
	per_1000 = df_codcount['Frequency'].sum()/1000
	df_codcount['Frequency'] = df_codcount['Frequency']/per_1000
	
	return df_codcount[['Codon','Amino_Acid','Frequency']]	#reorder the columns
	
if __name__=='__main__':

	about = 'Computes codon usage per 1000 of a whole transcriptome. Written by Rhondene Wint, rwint@ucmerced.edu.'
	epi_note = 'To contact the author about problems or errors,make a pull request at https://github.com/rhondene/Codon-Usage-in-Python'
	parser = argparse.ArgumentParser(description=about,epilog=epi_note)
	parser.add_argument('-CDS', help='Path to fasta file with species coding sequences', type=str, required=True, metavar='')
	parser.add_argument('-out', help='Path of destination folder for output csv file', type=str, default='./file.cu.csv', metavar='')	

	args=parser.parse_args()
	
	seqs=preproc(args.CDS)##formats fasta into csv of sequences
	CU = get_cod_freq(seqs)	 ##computes codon usage per 1000

	#save the file
	CU.to_csv('{}.cu.csv'.format(args.out), index=False)


	
	
	
	
	
	