##Author: Rhondene Wint
##Compute the 59 relative synonymous codon usage (RSCU) from a fasta file or list of fasta file
## I wrote this because the present codon usage tools could not handle multiple files or
## the output was in the form of standard codon usage table that needs extra parsing



# map codons to amino acids
## break up 6-codon family into 2 and 4 fold subsets(Ser (S), L(Leu), R (Arg))

codon_aa_2 = {
           
    "UCU":"S4", "UCC":"S4", "UCA":"S4", "UCG":"S4",
    "AGU":"S2", "AGC":"S2",
    "CUU":"L4", "CUC":"L4", "CUA":"L4", "CUG":"L4",
    "UUA":"L2", "UUG":"L2",
	"CGU":"R4", "CGC":"R4", "CGA":"R4", "CGG":"R4",
    "AGA":"R2", "AGG":"R2",
    "UUU":"F", "UUC":"F",  
    "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
    "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
    
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
   
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

####--------HELPER FUNCTIONS---------####

## preprocess coding sequence fasta file into a list

def get_seqs(filename):
	"""filename: path to file
	NOTE: ensure that newlines occur only between header lines and sequences,
	and not within sequences. If you're uncertain use my fix_fasta function beforehand"""
    
	seqs = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('>') is False:
                seqs.append(line.rstrip())
    return seqs
	

## compute codon count
def get_cod_freq(seqs):
    """ seqs: list of sequences of each Coding sequence """
    
    codon_count=dict() 
    
    for codon in list(codon_aa_2.keys()):
        if codon not in ['AUG', "UAA","UAG", "UGA", "UGG" ]:
            codon_count[codon]=0
    ##count codons 59 codons in each sequence, 
    for cds in seqs:
        cds = cds.upper().replace('T','U')
        codons = []
        ##make a list of codons
        for c in range(0,len(cds),3):
            if len(cds)%3 ==0:
                cod=cds[c:c+3]
                if 'N' not in cod:  ##ignore N and seqs not multiple of 3
                    codons.append(cod)
            else:
                continue

        for c in list(codon_count.keys()):
            codon_count[c]+= codons.count(c)
    
    df_rscu=pd.DataFrame(list(codon_count.items()) )
    df_rscu.columns=['Codon', 'Observed_Freq']
    df_rscu['Amino_Acid'] = [codon_aa_2[codon] for codon in df_rscu['Codon'].values]
    
    return df_rscu	## with absolute codon frequencies

## computes RSCU
def compute_rscu_weights(df_rscu):
    """ wij = RSCUij/ RSCU i,max"""
    aa_groups = df_rscu.groupby('Amino_Acid')
    aa =  df_rscu['Amino_Acid'].unique()  #make a list of all amino acids to iterate over
    df_list = []
    for a in aa:
        d=aa_groups.get_group(a)
        d['RSCU'] = d['Obs_Freq'].values/d['Obs_Freq'].mean() #obs/expected freq 
        d['Relative_Adaptive_Weights'] = d['RSCU'].values/d['RSCU'].max() 
        d['optimal'] = [True if rscu==d['RSCU'].max() else False for rscu in d['RSCU'].values] #marks optimal codon
        #d['Species'] = species
        df_list.append(d)
    return pd.concat(df_list)
	

for species in names:  #names is a list of filenames of fasta files
    seqs = get_seqs(species)  ##formats fasta into list of sequences
    df_rscu = get_cod_freq(seqs)  ##computes absolute codon frequencies
    rscu = compute_rscu_weights(df_rscu)  ##computes RSCU and adaptive weights
##saves final table with RSCU, adaptive weights and codon frequencies	
    rscu.to_csv('../RSCU_genomes/{}_rscu.csv'.format(species),index=False, sep=',')
	

