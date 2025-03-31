import os
import sys
import numpy as np
import time
from Bio import SeqIO
from modules.utils.TextColor import TextColor
from modules.utils.FileManager import FileManager

LONG_DELETION_LENGTH = 50

def pileup(paf, genome_size):
    """[summary]
    
    Arguments:
        filename {paf file name} -- db alignment result
        Genome_size {integer} -- TGS assembly only genome
    
    Returns:
        allele_count: genome each base [0-3]: ATCG, [4]: deletion
        coverage: genome each base coverage 
        Ins: genome each base [0-2]: 1st Ins, 2nd Ins, 3rd Ins; [0-3]: ATCG
    """
    ins_len = 7
    arr = np.zeros((genome_size, 5), dtype=np.int64)
    coverage = np.zeros(genome_size, dtype=np.int64)
    ins = np.zeros((genome_size, ins_len, 4), dtype=np.int64)
    
    over_ins = []  
    with open(paf, 'r') as f:        
        for line in f:
            line = line.split()            
            t_start = line[7] #reference
            if line[11] != '0': #mapping quality != 0
                cigar = line[-1]
                start_pos = int(t_start)                
                flag = 0
                longdel_count = 0
                longdel_status = 0
                for i in cigar: # ex. =ATC*c
                    
                    if i == 'A' or i == 'a':  # A:0, T:1, C:2, G:3
                        base = 0
                    elif i == 'T' or i == 't':
                        base = 1
                    elif i == 'C' or i == 'c':
                        base = 2
                    elif i == 'G' or i == 'g':
                        base = 3
                    
                    if i == '=': #match:1
                        flag = 1
                    elif i == '*': #mismatch:2
                        flag = 2
                        mismatch = 0
                    elif i == '+': #insertion
                        flag = 3
                        ins_pos = 0
                        over_ins = []
                    elif i == '-': #deletion
                        flag = 4
                        longdel_status = 0
                    
                    elif flag == 1:
                        longdel_count = 0
                        arr[start_pos][base] += 1 #arr1[100][0] = 1
                        coverage[start_pos] += 1
                        start_pos += 1 
                    elif flag == 2: 
                        # *gc
                        # -01
                        # 01
                        longdel_count = 0
                        if mismatch != 1:
                            mismatch += 1
                        else:
                            arr[start_pos][base] += 1 #Mismatch position
                            coverage[start_pos] += 1 
                            start_pos += 1
                        
                    elif flag == 3: 
                        #+AAAAA
                        #-0123
                        #01234
                        longdel_count = 0
                        if ins_pos < ins_len:
                            ins[start_pos-1][ins_pos][base] += 1
                            ins_pos += 1
                            over_ins.append(base)
                        elif ins_pos == ins_len:
                            for x,y in zip(range(ins_len), over_ins):
                                ins[start_pos-1][x][y] -= 1
                            over_ins = []
                            
                            
                    elif flag == 4:                    
                        if longdel_status == 0:
                            longdel_count += 1
                        if longdel_count > LONG_DELETION_LENGTH and longdel_status == 0:
                            for i in range(1,LONG_DELETION_LENGTH + 1):
                                arr[start_pos-i][4] -= 1
                                coverage[start_pos-i] -= 1                                                  
                            longdel_status = 1
                            longdel_count = 0
                        elif longdel_status != 1:
                            arr[start_pos][4] += 1
                            coverage[start_pos] += 1
                        start_pos+=1
        return arr, coverage, ins

def make_output_dir(type, output_dir, contig_id=None):
    if type=='contig':
        contig_output_dir = output_dir + '/' + contig_id
        contig_output_dir = FileManager.handle_output_directory(contig_output_dir)
        return contig_output_dir
    else:
        output_dir_debug = output_dir + '/' + type
        output_dir_debug = FileManager.handle_output_directory(output_dir_debug)
        return output_dir_debug

def align(draft, minimap_args, threads, db, path, reference=None):

    t_start=time.time()

    if reference:
        paf = '{}/truth.paf'.format(path)
        minimap2_cmd= 'minimap2 -cx asm5 --cs=long -t {thread} {draft} {reference} > {paf}'.format(thread=threads, draft=draft, reference=reference,paf=paf)
        
        out = path + '/truth_ANI.txt'
        Ani_cmd = 'fastANI -q {draft} -r {reference} -o {out}'.format(draft=draft,reference=reference,out=out)
        os.system(Ani_cmd)
    else:
        paf = '{}/contig.paf'.format(path)
        minimap2_cmd = 'minimap2 -cx {asm} --cs=long -t {thread} {draft} {db} > {paf}'\
            .format(asm=minimap_args, thread=threads, draft=draft, db=db, paf=paf)

    os.system(minimap2_cmd)
    t_end = time.time()
    print(t_end-t_start)

    return paf

