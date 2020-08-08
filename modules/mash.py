import os
import sys

def screen(contig_name, sketch_path, threads, output_dir, mash_threshold, download_contig_nums, contig_id):
    
    mash_out = '{}/{}.sort.tab'.format(output_dir, contig_id)
    mash_cmd = 'mash screen -p {thread} -i {ratio} {db} {draft} > temp.tab'\
            .format(thread=threads, ratio=mash_threshold, db=sketch_path, draft=contig_name)    
    sort_cmd = 'sort -gr temp.tab > temp.sort.tab'
    head_cmd = 'head -n {num} temp.sort.tab > {out}'.format(num=download_contig_nums, out=mash_out)
    os.system(mash_cmd)
    os.system(sort_cmd)
    os.system(head_cmd)
    os.remove('temp.tab')
    os.remove('temp.sort.tab')
    return mash_out


def dist(contig_name, sketch_path, threads, output_dir,mash_threshold ,download_contig_nums, contig_id):
    
    ratios=1- float (mash_threshold)
    mash_out = '{}/{}.sort.tab'.format(output_dir, contig_id)
    mash_cmd = 'mash dist -p {thread} -d {ratio} {db} {draft} > temp.tab'\
            .format(thread=threads, ratio=str (ratios),db=sketch_path, draft=contig_name)   
    sort_cmd = 'sort -gk3 temp.tab > temp.sort.tab'
    head_cmd = 'head -n {num} temp.sort.tab > {out}'.format(num=download_contig_nums, out=mash_out)
    os.system(mash_cmd)
    os.system(sort_cmd)
    os.system(head_cmd)
    os.remove('temp.tab')
    os.remove('temp.sort.tab')
    return mash_out
    
    
def get_ncbi_id(mashfile, mash_screen=None):
    ncbi_id = []
    with open(mashfile,'r') as f:
        for line in f:
            line = line.split()
            if mash_screen:
                ncbi_id.append(line[4]) #Use mash screen
            else:
                ncbi_id.append(line[0]) #Use mash dist
    return ncbi_id