
#!/usr/bin/python

import sys
import argparse
import pandas as pd


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
	


def read_fasta(fasta):
    headers=[]
    seqs=[]
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line.rstrip().split(';')[0].split(' ')[0])
            else:
                seqs.append(line.rstrip())
        assert len(headers)==len(seqs),' Number of IDs and sequences do not match'

    return (headers,seqs)
    


def get_cod_freq(gene):
    "computes absolute codon counts in the input gene coding sequence"
    ID = gene['ID'].values[0]
    cds = gene['CDS'].values[0].upper().replace('T','U')
    codon_count=dict() 
    
    #ignore 1-fold Met and Trp, along with stop codons
    for codon in list(codon_aa.keys()):
        if codon not in ['AUG', "UAA","UAG", "UGA", "UGG" ]:
            codon_count[codon]=0
    ##count codons 59 codons in each sequence, 
    codons = []
    
   ##count codons in cds
    codons = []
    for c in range(0,len(cds),3):   #O(len cds)
        cod=cds[c:c+3]
        try:
            codon_count[cod]+=1
        except KeyError:
            continue
        
    #store the fractional freq of each codon in the codon dictionary
    total_cod=len(cds)/3 #total number of codons in the cds
    for c in list(codon_count.keys()):    #O(len codondict)
        codon_count[c]/=total_cod
    
    df_codcnt=pd.DataFrame(list(codon_count.items()) )
    df_codcnt.columns=['Codon', 'Obs_Freq']
    df_codcnt['Amino_acid'] = [codon_aa[codon] for codon in df_codcnt['Codon'].values]
    df_codcnt['Length']= len(cds)
    df_codcnt['GeneID']= ID
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
        d['Relative_Adaptive_Weights'] = d['RSCU'].values/d['RSCU'].max() 
        df_list.append(d)
        rscu = pd.concat(df_list).fillna(0)  #some genomes may not use any amino acids
    return rscu # rscu of a gene
	
cols2= ['Phe-UUU', 'Phe-UUC', 'Ser4-UCU', 'Ser4-UCC', 'Ser4-UCA',
       'Ser4-UCG', 'Ser2-AGU', 'Ser2-AGC', 'Leu4-CUU', 'Leu4-CUC',
       'Leu4-CUA', 'Leu4-CUG', 'Leu2-UUA', 'Leu2-UUG', 'Tyr-UAU',
       'Tyr-UAC', 'Cys-UGU', 'Cys-UGC', 'Arg4-CGU', 'Arg4-CGC',
       'Arg4-CGA', 'Arg4-CGG', 'Arg2-AGA', 'Arg2-AGG', 'Pro-CCU',
       'Pro-CCC', 'Pro-CCA', 'Pro-CCG', 'His-CAU', 'His-CAC', 'Gln-CAA',
       'Gln-CAG', 'Ile-AUU', 'Ile-AUC', 'Ile-AUA', 'Thr-ACU', 'Thr-ACC',
       'Thr-ACA', 'Thr-ACG', 'Asn-AAU', 'Asn-AAC', 'Lys-AAA', 'Lys-AAG',
       'Val-GUU', 'Val-GUC', 'Val-GUA', 'Val-GUG', 'Ala-GCU', 'Ala-GCC',
       'Ala-GCA', 'Ala-GCG', 'Asp-GAU', 'Asp-GAC', 'Glu-GAA', 'Glu-GAG',
       'Gly-GGU', 'Gly-GGC', 'Gly-GGA', 'Gly-GGG']
	   
	   

if __name__=='__main__':

    #print a welcome message with my contacts
    msg = """Welcome to PyRSCU! All rights reserved @2020. \n To send your questions and suggestions, you may make a pull request at github OR email rwint@ucmerced.edu """

    parser = argparse.ArgumentParser(description="Computes the relative synonymous codon usage (RSCU) for each coding sequence. Written by Rhondene Wint (rwint@ucmerced.edu).")
    parser.add_argument('fasta', help='Name or path of of input fasta file.')
    seq_df = pd.read_csv(,sep=',' )
    print(seq_df.columns)
    df_list =[]
    #iterate over each CDS/gene in table

    for i in range(seq_df.shape[0]):
        gene = seq_df.iloc[[i]]
        if len(gene['CDS'].values[0])%3 !=0:
            print('Gene {} not multiple of 3'.format(gene['ID'].values[0]))
            continue
        
        df_codcnt = get_cod_freq(gene)
        rscu = compute_rscu_weights(df_codcnt)
        df_list.append(rscu) ## append RSCU matrix of each gene

    ##collect and stack rscu matrices for cds 
    omit = ['Obs_Freq','Amino_acid','Length','Relative_Adaptive_Weights', 'GeneID']
    comb = []
    for rscu in df_list:  #rname rscu_list
        r3 = rscu.set_index('Codon').drop(omit, axis=1).T
        r3.reset_index(drop=True)
        r3.columns=cols2  ##codon with amino acid
        
            ##add a gene information column
        r3['GeneID']= rscu['GeneID'].values[0]
        r3['Length']=rscu['Length'].values[0]; 
      
        comb.append(r3)

    pd.concat(comb,axis=0).reset_index(drop=True).to_csv(output+'_rscu.csv',index=False)
