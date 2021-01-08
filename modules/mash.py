import os
import sys
from modules.utils.TextColor import TextColor

def screen(contig_name, sketch_path, threads, output_dir, mash_threshold, contig_id=None, download_contig_nums=None,):
    mash_out = '{}/{}.sort.tab'.format(output_dir, contig_id)
    mash_cmd = 'mash screen -p {thread} -i {ratio} {db} {draft} > temp.tab'\
        .format(thread=threads, ratio=mash_threshold, db=sketch_path, draft=contig_name)
    os.system(mash_cmd)

    if contig_id=="meta":
        sort_cmd = 'sort -gr temp.tab > {out}'.format(out=mash_out)
        os.system(sort_cmd)
    else:
        sort_cmd = 'sort -gr temp.tab > temp.sort.tab'
        os.system(sort_cmd)
        head_cmd = 'head -n {num} temp.sort.tab > {out}'.format(num=download_contig_nums, out=mash_out)
        os.system(head_cmd)
        os.remove('temp.sort.tab')

    os.remove('temp.tab')

    return mash_out


def dist(contig_name, sketch_path, threads, output_dir, mash_threshold ,download_contig_nums, contig_id):
    ratios=1- float (mash_threshold)
    mash_out = '{}/{}.sort.tab'.format(output_dir, contig_id)
    mash_cmd = 'mash dist -p {thread} -d {ratio} {db} {draft} > {output_dir}/temp.tab'\
            .format(thread=threads, ratio=str (ratios),db=sketch_path, draft=contig_name, output_dir=output_dir)
    sort_cmd = 'sort -gk3 {output_dir}/temp.tab > {output_dir}/temp.sort.tab'.format(output_dir=output_dir)
    head_cmd = 'head -n {num} {output_dir}/temp.sort.tab > {out}'.format(num=download_contig_nums, out=mash_out, output_dir=output_dir)
    os.system(mash_cmd)
    os.system(sort_cmd)
    os.system(head_cmd)
    os.remove('{output_dir}/temp.tab'.format(output_dir=output_dir))
    os.remove('{output_dir}/temp.sort.tab'.format(output_dir=output_dir))
    return mash_out

    
def get_ncbi_id(mashfile, mash_screen=None):
    ncbi_id = []
    with open(mashfile,'r') as f:
        for line in f:
            line = line.split('\t')
            if mash_screen:
                ncbi_id.append(line[4]) #Use mash screen
            else:
                ncbi_id.append(line[0]) #Use mash dist
    return ncbi_id

def meta_get_ncbi_id(mashfile, download_contig_nums):
    download_2dlist = []

    with open(mashfile, 'r') as f:
        for line in f:
            line = line.split('\t')
            test = line[5].split(" ")
            genus = ""

            chromosome_type = False
            if "_genomic.fna.gz" in line[4]:  # chromosome
                chromosome_type = True
                if line[5][0] != "[":
                    if test[2]=='sp.':
                        continue
                    genus = test[1] + " " + test[2]
                else:
                    if test[4]=='sp.':
                        continue
                    genus = test[3] + " " + test[4]


            else:  # plasmid
                genus = test[0] + " " + test[1]
                continue

            if genus[-1] == ',':
                genus = genus[:-1]

            chk = False
            for i, x in enumerate(download_2dlist):
                if genus in x:
                    chk = True
                    if chromosome_type == True and download_2dlist[i][1]<int(download_contig_nums):
                        download_2dlist[i][1] += 1
                        download_2dlist[i].append(line[4])
                    elif chromosome_type == False and download_2dlist[i][2]<int(download_contig_nums):
                        download_2dlist[i][2] += 1
                        download_2dlist[i].append(line[4])

            if chk == False:
                download_2dlist.append([genus, 1, 0])
                download_2dlist[-1].append(line[4])


    ncbi_id = []
    for i, x in enumerate(download_2dlist):
        if x[1] >= 5:
            print(x[0]+" "+str(x[1]))
            for i in range(3, len(x)):
                ncbi_id.append(x[i])

    return ncbi_id
