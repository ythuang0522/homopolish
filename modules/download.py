import os
import sys
import urllib.request
import requests
import multiprocessing
from modules.utils.FileManager import FileManager
import time
import socket


#socket.setdefaulttimeout(30)


def checkInternetRequests(url='http://www.google.com/', timeout=3):
    try:
        r = requests.head(url, timeout=timeout)
        return True
    except requests.ConnectionError as ex:
        print(ex)
        return False


def run_process(id, url, path):
    path = path + '/'
    print("Downloaded " + id)

    if 'ftp' in url:
        filename = '{}{}'.format(path, id)
        try:
            urllib.request.urlretrieve(url, filename)
            pass

        except urllib.error.URLError as e:
            #print(e)
            for i in range(1, 5):
                
                if os.path.isfile(filename):
                    break
                else:
                    print("URL ERROR: time out !!!!!!!!")
                    print("sleep 10 sec...")
                    time.sleep(10)
                    print("Download again............"+str(i))
                    urllib.request.urlretrieve(url,filename)  
            

        except urllib.error.ContentTooShortError as e:
            print ('Network conditions is not good. Reloading...')
            run_process(id, url, path)

    else:
        try:
            r = requests.get(url, allow_redirects=True)
            r.raise_for_status()
            open('{}{}.fasta'.format(path, id), 'wb').write(r.content)
            
            pass

        except requests.exceptions.RequestException as e:
            print(e)
            print("sleep 10 sec...")
            time.sleep(10)
            return (id, url)


 
def parser_url(ncbi_id):  #GCF_002060415.1_ASM206041v1_genomic.fna.gz
    url_list = []
    for filename in ncbi_id:
        if '_genomic.fna.gz' in filename: #download chromosome
            id = filename.split('_genomic.fna.gz')[0]  #GCF_002060415.1_ASM206041v1
            gcf = id.split('_')[0]                      #GCF
            second = id.split('_')[1]                   #002060415.1
            number = second.split('.')[0]               #002060415
            letter = '/'.join(number[i:i+3] for i in range(0, len(number), 3))         #002/060/415/
            url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{gcf}/{letter}/{id}/{filename}".format(gcf=gcf, letter=letter, id=id, filename=filename) 
            url_list.append(url)
        else: #download plasmid
            url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id}&rettype=fasta'.format(id=filename)
            url_list.append(url)
    return url_list

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

    cnt = 0

    results = [(ncbi_id[i], url_list[i]) for i in range(0, len(ncbi_id))]
    loss = []
    with multiprocessing.Pool(cpus if cpus < max_pool_size else max_pool_size) as pool:
        for pair in results:
            cnt += 1
            pool.apply_async(run_process, args=(pair[0], pair[1], db_dir))
        pool.close()
        pool.join()

    results = list(filter(None, results))
    

    file_path = db_dir + '/*'
    db_path = path + '/All_homologous_sequences.fasta.gz'

    os.system('cat {} > {}'.format(file_path, db_path))
    return db_path
