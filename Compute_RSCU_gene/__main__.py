
#!/usr/bin/python

import argparse
import os
import pandas as pd
import warnings
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

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
    "GGU":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly",}
    

def get_cod_freq(gene,ID):
    "computes absolute codon counts in the input gene coding sequence"
    codon_count=dict() 
    
    #ignore 1-fold Met and Trp, along with stop codons
    non_deg=['AUG', "UAA","UAG", "UGA", "UGG" ]
    
    for codon in list(codon_aa.keys()):
        if codon not in non_deg :
            codon_count[codon]=0
    #make a list of codons
    codons=[]
    for c in range(0,len(gene),3):
        cod=gene[c:c+3]
        if 'N' not in cod:  ##ignore N
            codons.append(cod)
        else:
            continue
    for c in list(codon_count.keys()):
        codon_count[c]+= codons.count(c)

    df_codcnt=pd.DataFrame(list(codon_count.items()) )
    df_codcnt.columns=['Codon', 'Obs_Freq']
    df_codcnt['Amino_acid'] = [codon_aa[codon] for codon in df_codcnt['Codon'].values]
    df_codcnt['Length']= len(gene)
    df_codcnt['SeqID']= ID
    return df_codcnt
    
#compute relative usagae from codon frequency
def compute_rscu_weights(df_codcnt):
    """ wij = RSCUij/ RSCU i,max"""
    aa_groups = df_codcnt.groupby('Amino_acid')
    aa =  df_codcnt['Amino_acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['RSCU'] = d['Obs_Freq'].values/d['Obs_Freq'].mean() #obs/expected freq,
        d['AA-Codon']=d['Amino_acid']+'-'+d['Codon'] 
        df_list.append(d)
        rscu = pd.concat(df_list).fillna(0)  #some genomes may not use any amino acids
    return rscu # rscu of a gene

       
if __name__=='__main__':

    about = 'Computes relative synonymous codon usage per coding sequeunce. Written by Rhondene Wint, rwint@ucmerced.edu.'
    epi_note = 'To contact the author about problems or errors,make a pull request at https://github.com/rhondene/Codon-Usage-in-Python'
    parser = argparse.ArgumentParser(description=about,epilog=epi_note)
    parser.add_argument('-CDS', help='Path to fasta file with species coding sequences', type=str, required=True, metavar='')
    parser.add_argument('-out', help='Path of destination folder for output file', type=str, default='./file_out.rscu', metavar='') 

    args=parser.parse_args()


    headers,seqs=fix_fasta.fix_fasta(args.CDS)##preprocess fasta to a paired list of headers and sequences
    df_list=[]
    for i in range(len(seqs)):
        gene = seqs[i].replace('T','U').upper()
        ID = headers[i]
        if len(gene)%3 !=0:
            print('Gene {} not multiple of 3'.format(ID))
            continue
        
        df_codcnt = get_cod_freq(gene,ID)
        df_codcnt.to_csv('df_count2.csv',index=False)
        rscu = compute_rscu_weights(df_codcnt)
        df_list.append(rscu) ## append RSCU matrix of each gene

    ##collect and stack rscu matrices for cds 
    omit = ['Codon','Obs_Freq','Amino_acid','Length', 'SeqID']
    comb = []
    for rscu in df_list:  #rname rscu_list
        r3 = rscu.set_index('AA-Codon').drop(omit, axis=1).T
        r3.reset_index(drop=True)
        
            ##add a gene information column
        r3['SeqID']= rscu['SeqID'].values[0]
        r3['Length']=rscu['Length'].values[0] 
      
        comb.append(r3)


    pd.concat(comb,axis=0).reset_index(drop=True).to_csv(args.out+'_rscu.csv',index=False)
