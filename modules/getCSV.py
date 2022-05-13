import pandas as pds
import modules.mod_polish as mp


def getFixPosCSV(fileName,FixPosAry):
   if(len(FixPosAry)==0):
     print("number of fix-position is zero")
     return
   finalAry = []
   for posData in FixPosAry:
    finalAry.append(posData[0])
    
   df = pds.DataFrame(finalAry)
   df = df.reset_index()
   df.columns =["","pos"]
   df.to_csv("modpolish/"+fileName+"_fix_position.csv",encoding="ascii",index=False)



def getT_val(T_AllAry,pos):
    T_val = ""

    T_ATCG = T_AllAry[pos]
    if(T_ATCG[0] !=0 or T_ATCG[5] !=0 ):
        T_val = "A"
    elif(T_ATCG[1] !=0 or T_ATCG[6] !=0):
        T_val = "T"
    elif(T_ATCG[2]!=0 or T_ATCG[7] !=0):
        T_val = "C"
    elif(T_ATCG[3] !=0 or T_ATCG[8] !=0):
        T_val = "G" 
    return T_val


def getFixMisPosCSV(fileName,fixData,T_AllAry,T_misAry,R_AllAry,S_AllAry,FixPosAry):
    MissPosAry = []
    fixPos_ary = []
       
    for posData in FixPosAry:
     if(posData[0]!=0):
       fixPos_ary.append(posData[0])
 
    for T_pos in T_misAry:

     if(T_pos[0] == 0):
      continue  
      
     T_val = getT_val(T_AllAry,T_pos[0]) 
      
     if(T_pos[0] not in fixPos_ary):
       #MissPosAry.append([T_pos[0],fixData.seq[T_pos[0]],T_val,R_AllAry[T_pos[0]][0],R_AllAry[T_pos[0]][1],R_AllAry[T_pos[0]][2],R_AllAry[T_pos[0]][3],R_AllAry[T_pos[0]][5],R_AllAry[T_pos[0]][6],R_AllAry[T_pos[0]][7],R_AllAry[T_pos[0]][8],S_AllAry[T_pos[0]][0],S_AllAry[T_pos[0]][1],S_AllAry[T_pos[0]][2],S_AllAry[T_pos[0]][3],S_AllAry[T_pos[0]][4]])
       
       subPattern = mp.getSpecialPattern(T_pos[0],fixData.seq)
       MissPosAry.append([T_pos[0],fixData.seq[T_pos[0]],subPattern]) #getSpecialPatternCSV
       
    df = pds.DataFrame(MissPosAry)
    df = df.reset_index()
    df.columns =["","pos","draft_Val","True_Val","Read_A+","Read_T+","Read_C+","Read_G+","Read_A-","Read_T-","Read_C-","Read_G-","Sib_A","Sib_T","Sib_C","Si_G","SibInsDel"]
    df.to_csv("modpolish/"+fileName+".csv",encoding="ascii",index=False)
  
    return MissPosAry
    
    
def getFixEorPosCSV(fileName,fixData,T_AllAry,T_misAry,R_AllAry,S_AllAry,FixPosAry): 
    ErrorPosAry = []
   
    for F_pos in FixPosAry: 
      T_val = getT_val(T_AllAry,F_pos[0])       
      if(F_pos[1] != T_val and T_val!="" ):
       ErrorPosAry.append([F_pos[0],fixData.seq[F_pos[0]],T_val,R_AllAry[F_pos[0]][0],R_AllAry[F_pos[0]][1],R_AllAry[F_pos[0]][2],R_AllAry[F_pos[0]][3],R_AllAry[F_pos[0]][5],R_AllAry[F_pos[0]][6],R_AllAry[F_pos[0]][7],R_AllAry[F_pos[0]][8],S_AllAry[F_pos[0]][0],S_AllAry[F_pos[0]][1],S_AllAry[F_pos[0]][2],S_AllAry[F_pos[0]][3],S_AllAry[F_pos[0]][4]])
    
    if(len(ErrorPosAry)!=0):
     df = pds.DataFrame(ErrorPosAry)
     df = df.reset_index()
     df.columns =["","pos","draft_Val","True_Val","Read_A+","Read_T+","Read_C+","Read_G+","Read_A-","Read_T-","Read_C-","Read_G-","Sib_A","Sib_T","Sib_C","Si_G","SibInsDel"]
     df.to_csv("Qscore/cntCSV/"+fileName+".csv",encoding="ascii",index=False)
    
    return ErrorPosAry
