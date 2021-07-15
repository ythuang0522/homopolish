import os
import sys
from modules.utils.TextColor import TextColor

def writeDbInTxt(db,path):
    temp = path + '/db.txt'
    with open(temp, "w") as outfile:
        outfile.write("\n".join(db))
    return temp
    
def putInList(path,url,id,All_db):
    path = path + '/'
    if 'ftp' in url:
        filename = '{}{}.fna.gz'.format(path, id)
        if filename not in All_db:
            All_db.append(filename)
    else:
        filename ='{}{}.fasta'.format(path, id)
        if filename not in All_db:
            All_db.append(filename)
    return All_db
def computeAni(draft,db_txt,path):
    out = path + '/ANI.txt'
    Ani_cmd = 'fastANI -q {draft} --rl {db_txt} -o {out}'.format(draft=draft,db_txt=db_txt,out=out)
    os.system(Ani_cmd)
    return out

def parseAni(db_txt):
    remove_file=[]
    count = 0
    with open(db_txt,'r') as f:
        
        for line in f:
            count += 1
            line = line.split()
            Ani = float(line[2])
            first = float(line[3])
            second = float(line[4])
            
            if second - first > 5 or Ani <= 99:
                remove_file.append(line[1])
    
    if count != len(remove_file):
        for i in remove_file:
            os.remove(i)
        
