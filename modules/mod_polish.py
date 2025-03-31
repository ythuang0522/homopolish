import sys
import os
import numpy as np
import time
import gzip
import pysam
import pandas as pds
import multiprocessing
import glob
import modules.download as dl
import modules.polish_interface as mlp
import modules.alignment as ma
import modules.getCSV as CSV
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from modules.utils.FileManager import FileManager
from modules.utils.TextColor import TextColor
from modules import ani
from modules.VAtypeClass import FixSNP


#save fasta when modpolish not work
def saveNofixSeq(fixData):
  filePath = fixData.output_dir+"debug/"+fixData.contig_id
  record = SeqRecord(
         fixData.seq,
         id=fixData.contig_id,
         description = fixData.gen_desc
  )
  contig_output_filePath=write_for_new_fasta(record,filePath,fixData.contig_id+"_modpolish")
  return contig_output_filePath

def saveDraftContig(fixData):
  filePath = fixData.output_dir+"debug/"+fixData.contig_id
  record = SeqRecord(
    fixData.seq,
    id=fixData.contig_id,
    description = fixData.gen_desc
  )
  contig_output_filePath=write_for_new_fasta(record,filePath,fixData.contig_id+"_contig")
  return contig_output_filePath

def starModpolsh(fixData,debug_mod):
  modpolish_filePath = ""
  contig_fix_ary = []
  fileName = fixData.draft_genome_file.split('/')[-1].split('.')[0]	
  
 

  #file path end must exist "/"
  if(fixData.output_dir[-1] != "/"):
       fixData.output_dir = fixData.output_dir+"/"


  #create filePath
  if not os.path.isdir(fixData.output_dir+"debug"):
       os.makedirs(fixData.output_dir+"debug")
  
  

  for contig in SeqIO.parse(fixData.draft_genome_file, 'fasta'):
      fixData.contig_id = contig.name#get fa Contig_id
      fixData.seq = contig.seq#get fa seq
      fixData.gen_desc = contig.description
      
      #create each contig debug filePath
      if not os.path.isdir(fixData.output_dir+"debug/"+ fixData.contig_id):
        os.makedirs(fixData.output_dir+"debug/"+ fixData.contig_id)

      #target contig 
      fixData.contig_path = saveDraftContig(fixData)

      filePath,mod_fix_flag = getPos(fixData,debug_mod)
      contig_fix_ary.append(mod_fix_flag)
      
      if(filePath != ""):
        modpolish_filePath = modpolish_filePath +" "+filePath
  
  #check each contig without modpolish
  if(not any(contig_fix_ary)):
     timestr = time.strftime("[%Y/%m/%d %H:%M]")
     sys.stderr.write(TextColor.PURPLE + str(timestr) + " INFO:Methyl position less than threshold!"+ "\n" + TextColor.END)  
  else:
     fileName = fileName+"_modpolish"
     
  catAllfasta(fixData,modpolish_filePath,fileName)
  
  if(debug_mod):
     shutil.rmtree(fixData.output_dir+"debug")


def fixFlagFn(fixData,fixPosAry,totalCovergae,fileName):
   filePath = fixData.output_dir+"debug/"+fixData.contig_id
   #use expected value: (fixPos/coverage)/length
   num_fixPos = len(fixPosAry)
   genLen = len(fixData.seq)   
   fixFlag = True
   average_cov  = totalCovergae/genLen
   fix_expected_val = (num_fixPos/average_cov)/genLen
   
   #saveFn csv
   Fndata = []
   Fndata.append([average_cov,genLen,num_fixPos,fix_expected_val])

   df = pds.DataFrame(Fndata)
   df = df.reset_index()
   if(len(Fndata)!=0):
     df.columns =["","avg_cov","genLen","num_fixPos","fix_expected_val"]
   df.to_csv(filePath+"/"+fileName+"_fnFlag.csv",encoding="ascii",index=False)
   
   
   
   if(fix_expected_val<0.000001):
     fixFlag = False
     
   return True



def getSibFile(fixData):
   file_path = ""
   flag = True

   #use siblings from user side
   if(fixData.sib_files != ""):
     for sibFile in fixData.sib_files:
        file_path = file_path +" "+sibFile
     
     print(file_path)
     db_path = fixData.output_dir+"debug/"+fixData.contig_id+"/All_homologous_sequences.fna.gz"
     os.system('cat {} > {}'.format(file_path, db_path))
   
   else:
    #download siblings homogenome file from NCBI
    flag = dlHomoFile(fixData)  
   
   return flag


def getPos(fixData,debug_mod):
    

    fileName = fixData.draft_genome_file.split('/')[-1].split('.')[0]
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star modpolish with sequence length: "+ str(len(fixData.seq))  + "\n" + TextColor.END)
    
    #download siblings files
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star download Homo sibling files"+ "\n" + TextColor.END)
    
    sib_flag = getSibFile(fixData)
    if(sib_flag == False):
      return saveNofixSeq(fixData),True
    #siblings array  
    sibPath = fixData.output_dir+"debug/"+fixData.contig_id
    H_misAry,H_AllAry = getPileUpAry(fixData,sibPath,sibPath+"/All_homologous_sequences.fna.gz")   
   
    #Reads array   
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star get Reads File data"+ "\n" + TextColor.END)
  
    if(fixData.bamFile != ""):
        R_misAry_bam,R_AllAry_bam,totalCovergae = MismatchPileup_read_bam(fixData.bamFile,len(fixData.seq))
    else:
        bamFile = getBamPileUp(fixData,fileName,fixData.thread,fixData.contig_path,fixData.reads_file)
        R_misAry_bam,R_AllAry_bam,totalCovergae = MismatchPileup_read_bam(bamFile,len(fixData.seq))

    #use special pattern
    if(fixData.spPattern == ""):
        fixary = getMisPosVal(R_misAry_bam,R_AllAry_bam,H_AllAry,fixData.seq) 
    else:
        fixary = fix_pos_homo_read_pattern(fixData,H_AllAry,R_misAry_bam,R_AllAry_bam,fixData.spPattern)
    
    #need to fix?
    #fixFlag = fixFlagFn(fixData,fixary,totalCovergae,fileName)
    #if(fixFlag == False):
    #return saveNofixSeq(fixData),False

    #star fix
    modpolish_filePath= fixProcess(fixary,H_AllAry,fixData,fileName)
    CSV.getFixPosCSV(fixData.contig_id,fixary,fixData.output_dir)

    return modpolish_filePath,True



  


def catAllfasta(fixData,genFilePath,fileName):
    subFilePath = fixData.draft_genome_file.split('/')[-1].split('.')[0]
    timestr =time.strftime("[%Y/%m/%d %H:%M]")
    # create a directory for each contig
    contig_output_dir = mlp.make_output_dir("contig", fixData.output_dir, subFilePath)
    outPutFilePath = contig_output_dir+"/"+fileName+".fa"
    os.system('cat {} > {}'.format(genFilePath,outPutFilePath)) 
    
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO: output dir: " + contig_output_dir + "\n" + TextColor.END)




def getBamPileUp(fixData:FixSNP,fileName,threads,fasta,fastq):
    filePath = fixData.output_dir+'debug/'+fixData.contig_id
    bam = filePath + '/reads_sorted.bam'
    os.system('minimap2 -ax asm5 --cs=long -t {thread} {draft} {reference} | '
              'samtools view -@ {thread} -S -b - | '
              'samtools sort -@ {thread} -o {bam} -'.format(thread=threads, draft=fasta, reference=fastq, bam=bam))
    os.system('samtools index {}'.format(bam))
    bam = filePath+'/reads_sorted.bam'
    return bam
    
    
def getPileUpAry(fixData:FixSNP,pafPath,asemberlyFile):
   filePath = fixData.output_dir+'debug/'+fixData.contig_id
   ma.align(fixData.contig_path,"asm5",fixData.thread,"",filePath,asemberlyFile)
   misAry,posData_ary = MismatchPileup(filePath+"/truth.paf",len(fixData.seq))#get SNP position ATCG array
   
   return misAry,posData_ary   
   
   

def dlHomoFile(fixData:FixSNP):
   fixFlag = True
   filePath = fixData.output_dir+'debug/'+fixData.contig_id
   ncbi_id =  mlp.mash_select_closely_related(fixData.sketch_path,False,fixData.thread,filePath,fixData.mash_threshold,fixData.dl_contig_nums,fixData.contig_path,fixData.contig_id)
   
   
   if(len(ncbi_id)<5):
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.PURPLE + str(timestr) + "Closely-related genomes less than 5, not to modpolish...\n" + TextColor.END) 
    fixFlag = False
   
   url_list =  dl.parser_url(ncbi_id)
   dl_path = dl.download(filePath,ncbi_id,url_list,fixData.contig_path,99,5)

   return fixFlag



def getSibVal(S_PosAry,S_percentRate):  # return Homogenome siblings pos value
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
                local_actg = ""

                for i in range(reads_all-4,reads_all+5):
                    local_actg += fixData.seq[i]
                for j in range(8):
                    n_sum += R_Ary_bam[reads_all][j][0]
                    q_all += R_Ary_bam[reads_all][j][1]
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
            local_actg = ""
            
            if(reads_all<4 or reads_all>(len(fixData.seq)-5)):
             continue
            for i in range(reads_all-4,reads_all+5):
                local_actg += fixData.seq[i]
            for j in range(8):
                n_sum += R_Ary_bam[reads_all][j][0]
                q_all += R_Ary_bam[reads_all][j][1] 
            if(n_sum == 0):
             continue
            if(q_all/n_sum < 10):
                P_val,patt = pattern(ST,local_actg)
                if(patt):
                    H_val = getSibVal(S_AllAry[reads_all],0.7)
                    if(H_val == ''):
                        ary = [reads_all,P_val]
                        pos_ary.append(ary)
                    #else:
                    elif(H_val != fixData.seq[reads_all]):
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
    
    contig_name = output_dir_debug + '/' + fileName + '.fasta'
    SeqIO.write(contig, contig_name, "fasta")
    return contig_name

def fixGem(fasta,posAry):
  for ary in  posAry:
     if(ary[1] != ""):
       str1 = str(fasta[:int(ary[0])])
       str2 = str(fasta[int(ary[0])+1:])
       fasta = str1+str(ary[1])+str2

  if(not isinstance(fasta, str)):
      return str(fasta)
  return fasta


def getSeq(fixAry,fixData):
  if(len(fixAry)>1):
    fixSeq = fixGem(fixData.seq,fixAry)   
     #save fasta
    record = SeqRecord(
         Seq(fixSeq),
         id=fixData.contig_id,
         description = fixData.gen_desc
    )
    return record
  else:
    record = SeqRecord(
         fixData.seq,
         id=fixData.contig_id,
         description = fixData.gen_desc
    )
    return record


def fixProcess(fixAry,S_arr,fixData,fileName):#fix the draft
    
    record = getSeq(fixAry,fixData)
    contig_output_filePath=write_for_new_fasta(record,fixData.output_dir+"debug",fixData.contig_id+"_modpolish")
    return contig_output_filePath

def MismatchPileup(file_name, genome_size):

    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO:star Homo files pileup with sequence length: "+ str(genome_size)  + "\n" + TextColor.END)
     

    LONG_DELETION_LENGTH = 50
    misAry = np.array([np.array([0,np.zeros(8)])for i in range(genome_size)])
    ins_len = 7
    arr = np.zeros((genome_size, 9), dtype=np.int64)
    coverage = np.zeros(genome_size, dtype=np.int64)
    ins = np.zeros((genome_size, ins_len, 4), dtype=np.int64)
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
                    if(start_pos>=genome_size):
                     break

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
    arr = np.zeros((genome_size, 9,3), dtype=np.int64)
    coverage = np.zeros(genome_size, dtype=np.int64)
    ins = np.zeros((genome_size, ins_len, 4), dtype=np.int64)
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
            if(start_pos>=genome_size):
               break

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
