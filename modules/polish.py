import time
import feather
import numpy as np
import pandas as pd
import more_itertools as mit
from Bio import SeqIO

nucleotide = {0:'A', 1:'T', 2:'C', 3:'G'}


def stitch(draft, result, path):
    d_position = []
    df = feather.read_dataframe(result)
    record = next(SeqIO.parse(draft, "fasta"))

    deletion = df[df['predict'] == 4].position.values #deletion 的df
    insertion = df[df['predict'] < 4] #只剩insertion
    
    d_temp = [list(group) for group in mit.consecutive_groups(deletion)]  
    for key in d_temp:
        if len(key) <= 6:
            d_position = np.concatenate((d_position, key), axis=0)

    i_position = df[df['predict'] < 4].position.values 
    seq = []
    seq.append('>{}_homopolish\n'.format(record.id))
    
    pos = np.concatenate((i_position, d_position), axis=0)
    sort_pos = np.sort(np.unique(pos))
    
    start = 0
    record = list(record.seq)

    for i in sort_pos.astype(int): 
        seq.append(''.join(record[start:i]))
        if i in i_position:                 
            match = insertion[insertion['position'] == i]
            if i not in d_position:            
                seq.append(record[i])              
            for index, row in match.iterrows():
                seq.append(nucleotide.get(row.predict))
        start = i+1
    seq.append(''.join(record[start:]))

    name = path + '/polished.fasta'
    polished = open(name,'w')
    polished.write(''.join(seq)) 
    polished.write('\n') 
    return name