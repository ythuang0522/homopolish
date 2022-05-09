import sys
import numpy as np
import time
import modules.download as dl
from modules.VAtypeClass import FixSNP
import gzip
import modules.polish_interface as mlp
from Bio import SeqIO
import modules.alignment as ma
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import pysam
import pandas as pds
from modules.utils.FileManager import FileManager
import multiprocessing
from modules import ani
from modules.utils.TextColor import TextColor
from modules.VAtypeClass import FixSNP
import os
import modules.getCSV as CSV




def fixFlagFn(fixData,fixPosAry,totalCovergae):
   #use expected value: (fixPos/coverage)/length
   num_fixPos = len(fixPosAry)
   genLen = len(fixData.seq)   
   fixFlag = True
   fix_expected_val = (num_fixPos/totalCovergae)/genLen
   print(num_fixPos)
   print(totalCovergae)
   print(genLen)
   print(fix_expected_val*100000)
   if((fix_expected_val*100000)<0.000001):
     fixFlag = False
     
   return fixFlag



def getPos(fixData):
    fixData,flag = getGenLength(fixData)
    fileName = fixData.draft_genome_file.split('/')[-1].split('.')[0]
    
    #download homogenome file
    dlHomoFile(fixData,fileName)  
    
    #Homogenomes array  
    print('load Homo')    
    H_misAry,H_AllAry = getPileUpAry(fixData,"Homo/"+fileName,"Homo/"+fileName+"/All_homologous_sequences.fna.gz")   
    
    
    #Reads array   
    print('load Read')
    if(fixData.bamFile != ""):
        R_misAry_bam,R_AllAry_bam,totalCovergae = MismatchPileup_read_bam(fixData.bamFile,len(fixData.seq))
    else:
        bamFile = getBamPileUp(fileName,16,fixData.draft_genome_file,fixData.reads_file)
        print(bamFile)
        R_misAry_bam,R_AllAry_bam,totalCovergae = MismatchPileup_read_bam(bamFile,len(fixData.seq))
    


    
    if(fixData.spPattern != ""):
        fixary = getMisPosVal(R_misAry_bam,R_AllAry_bam,H_AllAry,fixData.seq) 
    else:
        fixary = fix_pos_homo_read_pattern(fixData,H_AllAry,R_misAry_bam,R_AllAry_bam,fixData.spPattern)
   
    
    #need to fix?
    fixFlag = fixFlagFn(fixData,fixary,totalCovergae)
    if(fixFlag == False):
      print("mismatch less than threshold!")
      return
    
    if(fixData.get_fixCSV_Flag == True):
      CSV.getFixPosCSV(fileName,fixary)
     
    if(fixData.get_EorCSV_Flag == True and fixData.get_fixCSV_Flag == True):
       T_misAry,T_AllAry = getPileUpAry(fixData,"fixDraft/toolsFile",fixData.true_genome_file)#True misAry
       FixEorPosAry = CSV.getFixEorPosAry(fileName+"_ErrorPos",fixData,T_AllAry,T_misAry,R_AllAry_bam,H_AllAry,fixary)
       FixMissPosAry = CSV.getFixMisPosAry(fileName+"Pattern_MissPos",fixData,T_AllAry,T_misAry,R_AllAry_bam,H_AllAry,fixary)      
    
    fixProcess(fixary,H_AllAry,fixData,fileName)
    del_file(fixData)

     
def del_file(fixData:FixSNP):
   for file in glob.glob(fixData.rmFilePath+"/*"):
     if(os.path.isfile(file)):
        os.remove(file)

def getGenLength(fixData):
  flag = True
  for contig in SeqIO.parse(fixData.draft_genome_file, 'fasta'):
    if(len(contig.seq)>fixData.genomeLen):
      fixData.contig_id = contig.name#get fa Contig_id
      fixData.seq = contig.seq#get fa seq
      flag = False      
  return fixData,flag

def getBamPileUp(fileName,threads,fasta,fastq):
    if(os.path.exists('Read') == False):
        os.mkdir('Read')
    if(os.path.exists('Read/'+fileName) == False):
        os.mkdir('Read/'+fileName)
    os.system('minimap2 -ax asm5 --cs=long -t {thread} {draft} {reference} > Read/{sam}/reads.sam'.format(thread=threads, draft=fasta, reference=fastq,sam=fileName))
    os.system('samtools view -S -b Read/{f}/reads.sam > Read/{f}/reads.bam'.format(f=fileName))
    os.system('samtools sort Read/{f}/reads.bam -o Read/{f}/reads_sorted.bam'.format(f=fileName))
    os.system('samtools index Read/{f}/reads_sorted.bam'.format(f=fileName))
    bam = 'Read/{f}/reads_sorted.bam'.format(f=fileName)
    return bam
    
    
def getPileUpAry(fixData:FixSNP,pafPath,asemberlyFile):

   if(os.path.exists('fixDraft') == False):
        os.mkdir('fixDraft')
   if(os.path.exists('fixDraft/toolsFile') == False):
        os.mkdir('fixDraft/toolsFile')
   ma.align(fixData.draft_genome_file,"asm5",10,"",pafPath,asemberlyFile)
   misAry,posData_ary = MismatchPileup(pafPath+"/truth.paf",len(fixData.seq))#get SNP position ATCG array
   
   return misAry,posData_ary   
   
   

def dlHomoFile(fixData:FixSNP,fileName):
   if(os.path.exists('Homo') == False):
        os.mkdir('Homo')
   if(os.path.exists('Homo/'+fileName) == False):
        os.mkdir('Homo/'+fileName)
   ncbi_id =  mlp.mash_select_closely_related(fixData.sketch_path,False,10,fixData.output_dir,fixData.mash_threshold,fixData.dl_contig_nums,fixData.draft_genome_file,fixData.contig_id)
   url_list =  dl.parser_url(ncbi_id)
   dl_path = dl.download("Homo/"+fileName,ncbi_id,url_list,fixData.draft_genome_file,99,5)




def getSibVal(S_PosAry,S_percentRate):  # return Homogenome pos value
    Homo_val = ""
    decision_Fix = False

    Homo_A = S_PosAry[0]+S_PosAry[5]
    Homo_T = S_PosAry[1]+S_PosAry[6]
    Homo_C = S_PosAry[2]+S_PosAry[7]
    Homo_G = S_PosAry[3]+S_PosAry[8]
    Homo_Ins_Del = S_PosAry[4]

    S_total = Homo_A+Homo_T+Homo_C+Homo_G
    if(Homo_A >= (S_total*S_percentRate)):
      decision_Fix = True 
    if(Homo_T >= (S_total*S_percentRate)):
      decision_Fix = True 
    if(Homo_C >= (S_total*S_percentRate)):
      decision_Fix = True 
    if(Homo_G >= (S_total*S_percentRate)):
      decision_Fix = True 

    if(decision_Fix):
      Homo_ATCG = [Homo_A,Homo_T,Homo_C,Homo_G,Homo_Ins_Del] #是否連ins del判斷比較好? 
      #Sib_ATCG = [Sib_A,Sib_T,Sib_C,Sib_G]        
      pos =  Homo_ATCG.index(max(Homo_ATCG))  
      Homo_val = getATCG_pos_use(pos)
      if(pos == 0 and Homo_A == 0):
      	Homo_val = ""

    return Homo_val

def pattern(ST,local_actg):
    patt = True
    if(ST == 'CCAGC'):
        if(local_actg[0:5] == "CCAAG"):
            value = 'C'
        elif(local_actg[1:6] == "CCGAC" or local_actg[1:6] == "CCAAC" or local_actg[1:6] == "CCAAG"):
            value = 'G'
        elif(local_actg[2:7] == "CCGGC" or local_actg[2:7] == "CCGAC"):
            value = 'A'
        elif(local_actg[2:7] == "GTCGG" or local_actg[2:7] == "GCCGG"):
            value = 'T'
        elif(local_actg[3:8] == "GTCGG" or local_actg[3:8] == "GTTGG" or local_actg[3:8] == "TTTGG"):
            value = 'C'
        elif(local_actg[4:9] == "TTTGG"):
            value = 'G'
        else:
            patt = False
            value = ''    
    elif(ST == 'GCAGC'):
        if(local_actg[1:6] == "GCGAC"):
            value = 'G'
        elif(local_actg[2:7] == "GCGAC" or local_actg[2:7] == "GCGGC"):
            value = 'A'
        elif(local_actg[2:7] == "GTCGC" or local_actg[2:7] == "GCCGC"):
            value = 'T'
        elif(local_actg[3:8] == "GTCGC"):
            value = 'C'
        else:
            patt = False
            value = ''
    else:
        patt = False
        value = ''
    return value,patt

def fix_pos_homo_read_pattern(fixData,S_AllAry,R_MisAry,R_Ary_bam,ST):
    pos_ary = []
    fix_pass = False
    count = 0
    for reads_all in range(0,len(fixData.seq)):
        if(R_MisAry[reads_all][0]!=0):
            S_val=""
            S_percentRate = 1
            diff_D_T = False
            R_precentFlag = True
            misPos = R_MisAry[reads_all]

            R_Pos_ATCG = R_Ary_bam[reads_all]
            R_ATCG_Total = R_Pos_ATCG[0][0]+R_Pos_ATCG[1][0]+R_Pos_ATCG[2][0]+R_Pos_ATCG[3][0]+R_Pos_ATCG[4][0]+R_Pos_ATCG[5][0]+R_Pos_ATCG[6][0]+R_Pos_ATCG[7][0]

            S_val = getSibVal(S_AllAry[misPos[0]],1)

            if((fixData.seq[misPos[0]] != S_val) and S_val !=""):
               diff_D_T = True
            if((R_Ary_bam[reads_all][0][0]+R_Ary_bam[reads_all][4][0])>(R_ATCG_Total*0.95) ):
              R_precentFlag = False
            elif((R_Ary_bam[reads_all][1][0]+R_Ary_bam[reads_all][5][0])>(R_ATCG_Total*0.95)):
              R_precentFlag = False
            elif((R_Ary_bam[reads_all][2][0]+R_Ary_bam[reads_all][6][0])>(R_ATCG_Total*0.95) ):
              R_precentFlag = False
            elif((R_Ary_bam[reads_all][3][0]+R_Ary_bam[reads_all][7][0])>(R_ATCG_Total*0.95) ):
              R_precentFlag = False
 
            if(S_val!="" and diff_D_T  and R_precentFlag):
               ary = [misPos[0],S_val]
               count += 1
               pos_ary.append(ary)             
            elif(reads_all>=4 or reads_all <=(len(fixData.seq)-5)):    
                q_all = 0
                n_sum = 0
                low_sum = 0
                local_actg = ""

                for i in range(reads_all-4,reads_all+5):
                    local_actg += fixData.seq[i]
                for j in range(8):
                    n_sum += R_Ary_bam[reads_all][j][0]
                    q_all += R_Ary_bam[reads_all][j][1]
                    low_sum += R_Ary_bam[reads_all][j][2]
                if(q_all/n_sum < 15):
                    P_val,patt = pattern(ST,local_actg)
                    if(patt):
                        H_val = getSibVal(S_AllAry[reads_all],0.7)
                        if(H_val == ''):
                            ary = [reads_all,P_val]
                            pos_ary.append(ary)
                        else:
                            ary = [reads_all,H_val]
                            pos_ary.append(ary)
        else:
            q_all = 0
            n_sum = 0
            low_sum = 0
            local_actg = ""
            
            if(reads_all<4 or reads_all>(len(fixData.seq)-5)):
             continue
            for i in range(reads_all-4,reads_all+5):
                local_actg += fixData.seq[i]
            for j in range(8):
                n_sum += R_Ary_bam[reads_all][j][0]
                q_all += R_Ary_bam[reads_all][j][1]
                low_sum += R_Ary_bam[reads_all][j][2]    
            if(n_sum == 0):
             continue
            if(q_all/n_sum < 10):
                P_val,patt = pattern(ST,local_actg)
                if(patt):
                    H_val = getSibVal(S_AllAry[reads_all],0.7)
                    if(H_val == ''):
                        ary = [reads_all,P_val]
                        pos_ary.append(ary)
                    else:
                        ary = [reads_all,H_val]
                        pos_ary.append(ary)
    return pos_ary

def getMisPosVal(R_MisAry,R_AllAry,S_AllAry,DraftSeq):
      posAry = []       
      
      for misPos in R_MisAry:
        if(misPos[0] == 0):
         continue

        S_val=""
        S_percentRate = 1
        diff_D_T = False
        R_precentFlag = True
     
        R_Pos_ATCG = R_AllAry[misPos[0]]
        R_ATCG_Total = R_Pos_ATCG[0][0]+R_Pos_ATCG[1][0]+R_Pos_ATCG[2][0]+R_Pos_ATCG[3][0]+R_Pos_ATCG[4][0]+R_Pos_ATCG[5][0]+R_Pos_ATCG[6][0]+R_Pos_ATCG[7][0]
 
        S_val = getSibVal(S_AllAry[misPos[0]],1)
              
  
        if((DraftSeq[misPos[0]] != S_val) and S_val !=""):
           diff_D_T = True

        if((R_AllAry[misPos[0]][0][0]+R_AllAry[misPos[0]][4][0])>(R_ATCG_Total*0.95) ):
          R_precentFlag = False   
        elif((R_AllAry[misPos[0]][1][0]+R_AllAry[misPos[0]][5][0])>(R_ATCG_Total*0.95)):
          R_precentFlag = False 
        elif((R_AllAry[misPos[0]][2][0]+R_AllAry[misPos[0]][6][0])>(R_ATCG_Total*0.95) ):
          R_precentFlag = False 
        elif((R_AllAry[misPos[0]][3][0]+R_AllAry[misPos[0]][7][0])>(R_ATCG_Total*0.95) ):
          R_precentFlag = False  
 
          
        if(S_val!="" and diff_D_T  and R_precentFlag):                 
           ary = [misPos[0],S_val]
           posAry.append(ary) 
         
      return posAry
      
def getATCG_pos_use(index:int):
   if(index == 0):
   	return "A"
   elif(index == 1):
   	return "T"
   elif(index == 2):
   	return "C"
   elif(index == 3):
   	return "G"
   else :
   	return ""

def write_for_new_fasta(contig, output_dir_debug,fileName):
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO: RUN-ID: " + contig.id + "\n" + TextColor.END)
    print(contig.id)
    print(output_dir_debug)
    # create a directory for each contig
    contig_output_dir = mlp.make_output_dir("contig", output_dir_debug, contig.id)
    # new fasta for each contig
    contig_name = contig_output_dir + '/' + fileName + '.fasta'
    SeqIO.write(contig, contig_name, "fasta")
    return contig_name, contig_output_dir

def fixGem(fasta,posAry):
  for ary in  posAry:
     if(ary[1] != ""):
       str1 = str(fasta[:int(ary[0])])
       str2 = str(fasta[int(ary[0])+1:])
       fasta = str1+str(ary[1])+str2

  if(not isinstance(fasta, str)):
      return str(fasta)
  return fasta


def fixProcess(fixAry,S_arr,fixData,fileName):        #fix the draft
    if(len(fixAry)>1):
     #fix genome
     fixSeq = fixGem(fixData.seq,fixAry)
     
     #save fasta
     record = SeqRecord(
         Seq(fixSeq),
         id=fixData.contig_id
     )
     if(os.path.exists('modpolish') == False):
        os.mkdir('modpolish')
     contig_name, contig_output_dir=write_for_new_fasta(record,"modpolish",fileName)

def MismatchPileup(file_name, genome_size):
    LONG_DELETION_LENGTH = 50
    misAry = np.array([np.array([0,np.zeros(8)])for i in range(genome_size)])
    ins_len = 7
    arr = np.zeros((genome_size, 9), dtype=np.int)
    coverage = np.zeros(genome_size, dtype=np.int)
    ins = np.zeros((genome_size, ins_len, 4), dtype=np.int)
    over_ins = [] 
    with open(file_name, 'r') as f:        
        for line in f:
            line = line.split()            
            t_start = line[7] #reference
           
            if line[11] != '0': #mapping quality != 0
               cigar = line[-1]
               start_pos = int(t_start)                
               flag = 0
               longdel_count = 0
               longdel_status = 0
               strain = line[4]#取正反股
               if(start_pos < genome_size):                 
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
                        coverage[start_pos] += 1
    
                        if(strain == "+"):
                          arr[start_pos][base] += 1  
                        if(strain == "-"):
                         pos = base +5
                         arr[start_pos][pos] += 1  
                        
                        start_pos += 1 
                    elif flag == 2: 
                        # *gc
                        # -01
                        # 01
                        longdel_count = 0
                        if mismatch != 1:
                          mismatch += 1
                          misAry[start_pos][0] = start_pos
                          
                          if(strain == "+"):
                            misAry[start_pos][1][base] += 1
                          if(strain == "-"): 
                            pos  = base + 4
                            misAry[start_pos][1][pos] += 1

                        else:
                          if(strain == "+"):
                            arr[start_pos][base] += 1
                          if(strain == "-"): 
                            pos  = base + 5
                            arr[start_pos][pos] += 1

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
  
        return misAry,arr
def MismatchPileup_read_bam(bam,genome_size):
    LONG_DELETION_LENGTH = 50
    misAry = np.array([np.array([0,np.zeros(8)])for i in range(genome_size)])
    ins_len = 7
    arr = np.zeros((genome_size, 9,3), dtype=np.int)
    coverage = np.zeros(genome_size, dtype=np.int)
    ins = np.zeros((genome_size, ins_len, 4), dtype=np.int)
    over_ins = []
    totalCovergae = 0
    bamData = pysam.AlignmentFile(bam,'rb')
    for each_read in bamData.fetch(until_eof=True):
        x = each_read.mapping_quality
        if(x == 0):
            continue
        for cs in range(len(each_read.tags)-1,0,-1):
            if(each_read.tags[cs][0] == 'cs'):
                cigar = each_read.tags[cs][1]
                break
        start_pos = each_read.pos
        base_quality = each_read.query_alignment_qualities
        pos_counter = 0
        flag = 0
        longdel_count = 0
        longdel_status = 0


        if(each_read.is_reverse == True):
            reverse = 4
        else:
            reverse = 0
        for i in cigar: # ex. =ATC*c
            #print('in_cigar')
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
                arr[start_pos][reverse + base][0] += 1
               
                arr[start_pos][reverse + base][1] += base_quality[pos_counter]
                if base_quality[pos_counter] < 10:
                    arr[start_pos][reverse + base][2] += 1
                coverage[start_pos] += 1
                pos_counter += 1
                start_pos += 1
                totalCovergae +=1
          
            elif flag == 2:
                # *gc
                # -01
                # 01
                longdel_count = 0
                if mismatch != 1:
                    mismatch += 1
                    misAry[start_pos][0] = start_pos
                    misAry[start_pos][1][reverse+base] += 1
                    
                else:
                    arr[start_pos][reverse+base][0] += 1
                    arr[start_pos][reverse + base][1] += base_quality[pos_counter]
                    if base_quality[pos_counter] < 10:
                        arr[start_pos][reverse + base][2] += 1
                    coverage[start_pos] += 1
                    start_pos += 1
                    pos_counter += 1
                    mismatch = 0
                    totalCovergae +=1
                    
            elif flag == 3:
                #+AAAAA
                #-0123
                #01234
                longdel_count = 0
               
                pos_counter += 1
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
                        arr[start_pos-i][8][0] -= 1
                        coverage[start_pos-i] -= 1
                        totalCovergae -=1
                    longdel_status = 1
                    longdel_count = 0
                elif longdel_status != 1:
                    arr[start_pos][8][0] += 1
                    coverage[start_pos] += 1
                    totalCovergae +=1
                start_pos+=1    
    return misAry,arr,totalCovergae


