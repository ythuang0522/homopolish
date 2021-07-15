import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

def preprocessing(df):
    new_draft = []
    df['A'] = df['A'].astype(int) / df['coverage'].astype(int)
    df['T'] = df['T'].astype(int) / df['coverage'].astype(int)
    df['C'] = df['C'].astype(int) / df['coverage'].astype(int)
    df['G'] = df['G'].astype(int) / df['coverage'].astype(int)
    df['gap'] = df['gap'].astype(int) / df['coverage'].astype(int)
    df['Ins_A'] = df['Ins_A'].astype(int) / df['coverage'].astype(int)
    df['Ins_T'] = df['Ins_T'].astype(int) / df['coverage'].astype(int)
    df['Ins_C'] = df['Ins_C'].astype(int) / df['coverage'].astype(int)
    df['Ins_G'] = df['Ins_G'].astype(int) / df['coverage'].astype(int)    
    df.loc[df['homopolymer'] > 10, 'homopolymer'] = 10 
    for i in range(0,11):
        if i not in df.homopolymer.value_counts():
            name = 'homopolymer_' + str(i)
            df[name] = 0
    df = pd.DataFrame(df.drop(['draft'], axis=1))
    df = pd.get_dummies(df, columns=['homopolymer'])
    df = pd.DataFrame(df.drop(['position'], axis=1))
    scaler = MinMaxScaler()
    df['coverage'] = scaler.fit_transform(df[['coverage']])
    return df
    
def haplotype(df):
    haplotype = np.zeros((df.shape[0],), dtype=int)
    ins_pos = df[df.draft == "-"].position.values
   
    side = 10
    dis = np.diff(ins_pos)
    idxs = np.where((dis < 10) & (dis > 1))[0]
    

    for i in idxs:        
        head = df[df['position'] == ins_pos[i]]
        head_total = head[['Ins_A','Ins_T','Ins_C','Ins_G']].sum(axis=1).values[0]
        tail = df[df['position'] == ins_pos[i+1]]
        tail_total = tail[['Ins_A','Ins_T','Ins_C','Ins_G']].sum(axis=1).values[0]
        
        head_idx=head.index.values[0]
        tail_idx=tail.index.values[0]
        
        coverage = head.coverage.values[0]
        if (head_total + tail_total) == coverage:   
            if head_total > tail_total:
                haplotype[head_idx] = 1
                haplotype[tail_idx] = 2
            elif head_total < tail_total:
                haplotype[head_idx] = 2
                haplotype[tail_idx] = 1
            else:
                haplotype[head_idx] = 1
                haplotype[tail_idx] = 1
    df['haplotype'] = haplotype
    #print(df['haplotype'])
    return df
'''
def haplotype(df):
    haplotype = np.zeros((df.shape[0],), dtype=int)
    ins_index = df[df.draft == '-'].index.values
    side = 10
    dis = np.diff(ins_index)
    idxs = np.where((dis < 10) & (dis > 1))[0]
    

    for i in idxs:        
        head = df[df.index == ins_index[i]]
        head_total = head[['Ins_A','Ins_T','Ins_C','Ins_G']].sum(axis=1).values[0]
        tail = df[df.index == ins_index[i+1]]
        tail_total = tail[['Ins_A','Ins_T','Ins_C','Ins_G']].sum(axis=1).values[0]
        coverage = head.coverage.values[0]
        if (head_total + tail_total) == coverage:   
            if head_total > tail_total:
                haplotype[ins_index[i]] = 1
                haplotype[ins_index[i+1]] = 2
            elif head_total < tail_total:
                haplotype[ins_index[i]] = 2
                haplotype[ins_index[i+1]] = 1
            else:
                haplotype[ins_index[i]] = 1
                haplotype[ins_index[i+1]] = 1
    df['haplotype'] = haplotype
    return df
'''