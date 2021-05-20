import os
import sys
import urllib.request 
import pycurl
import wget
import requests
import multiprocessing
from modules.utils.FileManager import FileManager
from subprocess import call
import shutil
import time
from modules.utils.TextColor import TextColor
from contextlib import closing

def checkInternetRequests(url, timeout=3):
    try:
        req = urllib.request.Request(url)
        response = urllib.request.urlopen(req)
        return True
    except urllib.URLError as ex:
        print(ex)
        return False

def download_NCBI(filename,url,):
    with open(filename, 'wb') as f:
                curl = pycurl.Curl()
                curl.setopt(pycurl.URL, url)
                curl.setopt(pycurl.WRITEDATA, f)
                curl.perform()
                curl.close()
                
    with open(os.devnull, 'w') as gzip:
            result = call(['gzip', '-t', filename], stdout=gzip, stderr=gzip)
            if result != 0: 
                #print("file is corrupted")
                raise OSError

def run_process(id, url, path):
    path = path + '/'
    filename = '{}{}.fna.gz'.format(path, id)
    if 'ftp' in url:
        try:
            
            ok = checkInternetRequests(url,3)
            if ok :
                print("Downloaded " + id)
            else :
                print('Connect is bad')
            
            #urllib.request.urlretrieve(url, filename)
            
            '''
            with closing(urllib.request.urlopen(url)) as r:
                with open(filename, 'wb') as f:
                    curl = pycurl.Curl()
                    curl.setopt(pycurl.URL, url)
                    curl.setopt(pycurl.WRITEDATA, f)
                    curl.perform()
                    curl.close()
                    #shutil.copyfileobj(r, f)
            '''
            download_NCBI(filename,url)
            
        except pycurl.error as e :
            time.sleep(3)
            run_process(id,url,path)
            print(e)
        except urllib.error.ContentTooShortError as e:
            print(e)
            run_process(id,url,path)
        
        except OSError as e:
            print("{} is corrupted".format(filename))
            #sys.exit(0)
            run_process(id,url,path)
            #sys.exit(0)
        
    else:
        try:
            '''
            r = requests.get(url, allow_redirects=True, verify = False)
            r.raise_for_status()
            '''
            
            print("Downloaded " + id)
           
            download_NCBI(filename,url)

            
        except requests.exceptions.RequestException as e:
            print(e)
            #sys.exit(0)
        except pycurl.error as e :
            time.sleep(3)
            run_process(id,url,path)
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

# def parser_genus_species(genus, ):
   
#    url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

#    response = requests.get(url)
#    html_doc = response.text

#    line = html_doc.split("\n")
#    ncbi_id = []
#    url_list = []

#    for i in range(2, len(line) - 1):
#        element = line[i].split("\t")
#        if genus in element[7]:
#            filename = element[19].split('/')[-1]
#            ftp = element[19] + "/" + filename + "_genomic.fna.gz"
#            ncbi_id.append(element[0])
#            url_list.append(ftp)
#    return ncbi_id, url_list

def parser_genus_species(genus_species, download_contig_nums=None):
    
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"

    response = requests.get(url)
    html_doc = response.text

    line = html_doc.split("\n")
    ncbi_id = []
    url_list = []
    allData = []

    if "_" in genus_species:
        genus_species = genus_species.replace("_", " ")

    for i in range(2, len(line) - 1):
        element = line[i].split("\t")

        if genus_species in element[7]:
            
            filename = element[19].split('/')[-1]
            #print(filename)
            ftp = element[19] + "/" + filename + "_genomic.fna.gz"

            if len(ncbi_id) < int(download_contig_nums):
                ncbi_id.append(element[0])
                url_list.append(ftp)
                allData.append(element[7])
            else:
                break

    # add homologous sequences (same genus) quantity up to download_contig_num
    if len(ncbi_id) < int(download_contig_nums):
        genus_species = genus_species.split(" ")

        for i in range(2, len(line) - 1):
            element = line[i].split("\t")

            if genus_species[1] not in element[7] and genus_species[0] in element[7]:
                filename = element[19].split('/')[-1]
                ftp = element[19] + "/" + filename + "_genomic.fna.gz"

                if len(ncbi_id) < int(download_contig_nums):
                    ncbi_id.append(element[0])
                    url_list.append(ftp)
                    allData.append(element[7])
                else:
                    break

    if len(ncbi_id) < 5:  # Would'nt polish if closely-related genomes less than 5
        sys.stderr.write(TextColor.PURPLE + "Closely-related genomes less than 5, not to polish...\n" + TextColor.END)
        return

    return ncbi_id, url_list 

def download(path, ncbi_id, url_list): 

    db_dir = path + '/homologous_sequences/'
    
    db_dir = FileManager.handle_output_directory(db_dir)
    max_pool_size = 3 #API rate limit exceeded, can't go higher
    cpus = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(cpus if cpus < max_pool_size else max_pool_size)
    for id, url in zip(ncbi_id, url_list):
        try:
            pool.apply_async(run_process, args=(id, url, db_dir))
        except multiprocessing.BufferTooShort as e:
            print(e)
        except multiprocessing.ProcessError as e:
            print(e)
        except multiprocessing.TimeoutError:
            print(e)
        except multiprocessing.AuthenticationError:
            print(e)
    pool.close()
    pool.join()

    file_path = db_dir + '/*'
    db_path = path + '/All_homologous_sequences.fna.gz'
    os.system('cat {} > {}'.format(file_path, db_path))
    print('')
    return db_path
