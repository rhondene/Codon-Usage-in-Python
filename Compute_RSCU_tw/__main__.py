"""Compute Relative Synonymous Codon Usage (RSCU) of  a transcriptome(Sharp and Li 1986)"""
import argparse
import pandas as pd
import warnings
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
import fix_fasta

## dictionary that maps codons to amino acids
"""break up 6-codon family into 2 and 4 fold (Ser (S), L(Leu), R (Arg))"""
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



def get_cod_freq(headers:list, seqs:list):
    """ 
        headers: list of headers
        seqs: list of CDS
    	Returns a 64-dim dataframe of total absolute codon frequencies
    """
    
    codon_count=dict()
	#ignore 1-fold Met and Trp, along with stop codons
    non_deg=['AUG', "UAA","UAG", "UGA", "UGG" ]
    codon_count=dict() 
    codon_count = {codon: 0 for codon in codon_aa if codon not in non_deg }
	
    for i,cds in enumerate(seqs):
        if len(cds)%3 !=0:
            ID = headers[i].split(' ')[0]
            print(f"WARNING: Skipping {ID} Length of CDS is not a multiple of 3.")
            continue
        cds = cds.upper().replace('T','U')
        ##count codons in cds
        for i in range(0, len(cds), 3):
            codon = cds[i:i+3]
            if codon in codon_count:
                codon_count[codon] += 1
    
    df_codcount=pd.DataFrame(list(codon_count.items()) )
    df_codcount.columns=['Codon', 'Obs_Freq']
    df_codcount['Amino_Acid'] = [codon_aa[codon] for codon in df_codcount['Codon'].values] ## add amino acid column
    
    return df_codcount

def compute_rscu_weights(df_codcount:pd.DataFrame)->pd.DataFrame:
    """ Caclculates Relative Synonymous codon usage (RSCU) wij = RSCUij/ RSCU i,max
    Returns: 59-dim dataframe of RSCU values for each codon """
    aa_groups = df_codcount.groupby('Amino_Acid')
    aa =  df_codcount['Amino_Acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['RSCU'] = d['Obs_Freq'].values/d['Obs_Freq'].mean() #obs/expected freq 
        d['Relative_Adaptive_Weights'] = d['RSCU'].values/d['RSCU'].max() 
        d['optimal'] = [True if rscu==d['RSCU'].max() else False for rscu in d['RSCU'].values] #marks optimal codon
        df_list.append(d)
    return pd.concat(df_list)


if __name__=='__main__':
	about = 'Computes transcriptome-wide Relative synonymous codon usage and absolute codon counts. Written by Rhondene Wint, rwint@ucmerced.edu.'
	epi_note = 'To contact the author about problems or errors,make a pull request at https://github.com/rhondene/Codon-Usage-in-Python'
	parser = argparse.ArgumentParser(description=about,epilog=epi_note)
	parser.add_argument('-CDS', help='Path to fasta file with species coding sequences', type=str, required=True, metavar='')
	parser.add_argument('-out', help='Path of destination folder for output file (text file)', type=str, default='./file_out.rscu', metavar='')	

	args=parser.parse_args()


	headers,seqs=fix_fasta.fix_fasta(args.CDS) ##formats fasta into csv of sequences
	df_codcount = get_cod_freq(headers,seqs)	  ##computes absolute codon frequencies
	rscu = compute_rscu_weights(df_codcount)  ##computes RSCU and adaptive weights

	#save the file
	rscu.to_csv('{}.rscu'.format(args.out), index=False)

   
                      
       
    
