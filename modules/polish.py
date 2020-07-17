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

    seq = []
    seq.append('>{}_homopolish\n'.format(record.id))

    for i in range(len(record)): #只對insertion/deletion 處理
        if i not in d_position: #不管deletion -> label = 4 
            if i in insertion['position'].values: #處理 insertion
                match = insertion[insertion['position'] == i]     
                seq.append(record[i])  
                for index, row in match.iterrows():                       
                    seq.append(nucleotide.get(row.predict))               
            else: #不動
                seq.append(record[i])
    name = path + '/polished.fasta'
    polished = open(name,'w')
    polished.write(''.join(seq)) 
    polished.write('\n') 
    return name