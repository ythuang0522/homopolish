import numpy as np
cimport numpy as np

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
    cdef int ins_len = 7
    cdef np.ndarray arr = np.zeros((genome_size, 5), dtype=np.int)
    cdef np.ndarray coverage = np.zeros(genome_size, dtype=np.int)
    cdef np.ndarray ins = np.zeros((genome_size, ins_len, 4), dtype=np.int)
    cdef list over_ins
    cdef int flag = 0
    cdef int base = 0 # A/0 T/1 C/2 G/3
    cdef int longdel_count = 0
    cdef int longdel_status = 0
    cdef int mismatch = 0
    cdef int ins_pos = 0 
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
                for i in cigar: 
                    
                    if i == 'A' or i == 'a':  
                        base = 0
                    elif i == 'T' or i == 't':
                        base = 1
                    elif i == 'C' or i == 'c':
                        base = 2
                    elif i == 'G' or i == 'g':
                        base = 3
                    
                    if i == '=': #match
                        flag = 1
                    elif i == '*': #mismatch
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
                        arr[start_pos][base] += 1 
                        coverage[start_pos] += 1
                        start_pos += 1 
                    elif flag == 2: 
                        longdel_count = 0
                        if mismatch != 1:
                            mismatch += 1
                        else:
                            arr[start_pos][base] += 1 
                            coverage[start_pos] += 1 
                            start_pos += 1
                        
                    elif flag == 3: 
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
