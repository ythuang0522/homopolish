import os
import sys
import urllib.request
import wget
import requests
import multiprocessing
from modules.utils.FileManager import FileManager


def checkInternetRequests(url='http://www.google.com/', timeout=3):
    try:
        r = requests.head(url, timeout=timeout)
        return True
    except requests.ConnectionError as ex:
        print(ex)
        return False

def run_process(id, url, path):
    path = path + '/'
    if 'ftp' in url:
        try:
            urllib.request.urlretrieve(url, path)
        except urllib.error.HTTPError as e:
            print(e)
    else:
        try:
            r = requests.get(url, allow_redirects=True)
            r.raise_for_status()
        except
            requests.exceptions.RequestException as e:
            print(e)
        open('{}{}.fasta'.format(path, id), 'wb').write(r.content)


def parser_url(ncbi_id):
    url_list = []
    for filename in ncbi_id:
        if '_genomic.fna.gz' in filename: #download chromosome
            id = filename.split('_genomic.fna.gz')[0]
            gcf = id.split('_')[0]
            second = id.split('_')[1]
            number = second.split('.')[0]
            letter = '/'.join(number[i:i+3] for i in range(0, len(number), 3))                   
            url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{gcf}/{letter}/{id}/{filename}".format(gcf=gcf, letter=letter, id=id, filename=filename) 
            url_list.append(url)
        else: #download plasmid                  
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id}&rettype=fasta'.format(id=filename)
            url_list.append(url)
    return url_list

def parser_genus(genus):
    
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

    response = requests.get(url)
    html_doc = response.text

    line = html_doc.split("\n")
    ncbi_id = []
    url_list = []

    for i in range(2, len(line) - 1):
        element = line[i].split("\t")
        if genus in element[7]:
            filename = element[19].split('/')[-1]
            ftp = element[19] + "/" + filename + "_genomic.fna.gz"
            ncbi_id.append(element[0])
            url_list.append(ftp)
    return ncbi_id, url_list
    
def download(path, ncbi_id, url_list): 

    db_dir = path + '/homologous_sequences/'
    db_dir = FileManager.handle_output_directory(db_dir)
    max_pool_size = 3 #API rate limit exceeded, can't go higher
    cpus = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(cpus if cpus < max_pool_size else max_pool_size)
    for id, url in zip(ncbi_id, url_list):
        pool.apply_async(run_process, args=(id, url, db_dir))
    pool.close()
    pool.join()

    file_path = db_dir + '/*'
    db_path = path + '/All_homologous_sequences.fna.gz'
    os.system('cat {} > {}'.format(file_path, db_path))
    print('')
    return db_path
