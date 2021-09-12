import os
import sys
import time
import shutil
from Bio import SeqIO
from modules import mash
from modules import download
from modules import alignment
from modules import align2df
from modules import prediction
from modules import polish
import numpy as np
from modules.utils.TextColor import TextColor
from modules.utils.FileManager import FileManager


def print_system_log(stage):
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) +" INFO: Stage: "+ stage + "\n" + TextColor.END)

def print_stage_time(stage, time):
    sys.stderr.write(TextColor.PURPLE + "TIME "+ stage + ": " + str(time) + "\n" + TextColor.END)
        
def get_elapsed_time_string(start_time, end_time):
    """
    Get a string representing the elapsed time given a start and end time.
    :param start_time: Start time (time.time())
    :param end_time: End time (time.time())
    :return:
    """
    elapsed = end_time - start_time
    # hours = int(elapsed / 60**2)
    mins = int(elapsed % 60**2 / 60)
    secs = int(elapsed % 60**2 % 60)
    time_string = "{} MINS {} SECS.".format(mins, secs)

    return time_string

def mash_select_closely_related(sketch_path, mash_screen, threads, output_dir, mash_threshold,
                           download_contig_nums, contig_name, contig_id):
    if mash_screen:
        mash_file = mash.screen(contig_name, sketch_path, threads, output_dir, mash_threshold, contig_id, download_contig_nums)
    else:
        mash_file = mash.dist(contig_name, sketch_path, threads, output_dir, mash_threshold, download_contig_nums,
                              contig_id)

    ncbi_id = mash.get_ncbi_id(mash_file, mash_screen)

    return ncbi_id

def homologous_retrieval(paf_name, genome_size, contig_output_dir, contig_id, contig_name, ref_paf=None):
    print_system_log('Homologous retrieval')
    homologous_retrieval_start_time = time.time()

    arr, coverage, ins = alignment.pileup(paf_name, genome_size)
    npz = '{}.npz'.format(contig_output_dir + '/' + contig_id)
    np.savez(npz, arr=arr, coverage=coverage, ins=ins)

    if ref_paf:
        arr, coverage, ins = alignment.pileup(ref_paf, genome_size)
        ref_npz = '{}/truth.npz'.format(contig_output_dir)
        np.savez(ref_npz, arr=arr, coverage=coverage, ins=ins)
        dataframe_path = align2df.todf(contig_name, npz, contig_output_dir, ref_npz)
    else:
        dataframe_path = align2df.todf(contig_name, npz, contig_output_dir)
    
    homologous_retrieval_end_time = time.time()
    homologous_retrieval_time = get_elapsed_time_string(homologous_retrieval_start_time, homologous_retrieval_end_time)
    print_stage_time('Homologous retrieval', homologous_retrieval_time)
    return dataframe_path



def make_output_dir(type, output_dir, contig_id=None):
    if type=='contig':
        contig_output_dir = output_dir + '/' + contig_id
        contig_output_dir = FileManager.handle_output_directory(contig_output_dir)
        return contig_output_dir
    else:
        output_dir_debug = output_dir + '/' + type
        output_dir_debug = FileManager.handle_output_directory(output_dir_debug)
        return output_dir_debug



def homopolish(contig_name, minimap_args, threads, db_path, model_path, contig_output_dir, dataframe):

    print_system_log('Prediction')
    predict_start_time = time.time()
    result = prediction.predict(dataframe, model_path, threads, contig_output_dir)
    predict_end_time = time.time()
    predict_time = get_elapsed_time_string(predict_start_time, predict_end_time)
    print_stage_time('Prediction', predict_time)

    print_system_log('Polish')
    polish_start_time = time.time()
    finish = polish.stitch(contig_name, result, contig_output_dir)
    polish_end_time = time.time()
    polish_time = get_elapsed_time_string(polish_start_time, polish_end_time)
    print_stage_time('Polish', polish_time)

    return finish



def write_for_new_fasta(contig, output_dir_debug):
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO: RUN-ID: " + contig.id + "\n" + TextColor.END)
    print(contig.id)
    print(output_dir_debug)
    # create a directory for each contig
    contig_output_dir = make_output_dir("contig", output_dir_debug, contig.id)
    # new fasta for each contig
    contig_name = contig_output_dir + '/' + contig.id + '.fasta'
    SeqIO.write(contig, contig_name, "fasta")
    return contig_name, contig_output_dir



def check_homopolish(paf, contig_name, contig_output_dir, contig, minimap_args, threads, db_path, model_path):
    if os.stat(paf).st_size != 0:
        record = SeqIO.read(contig_name, "fasta")
        genome_size = len(record)
                
        dataframe = homologous_retrieval(paf, genome_size, contig_output_dir, contig.id, contig_name)

        if dataframe==False:
            return contig_name
        #run homopolish
        finish = homopolish(contig_name, minimap_args, threads, db_path, model_path, contig_output_dir, dataframe)
        
    else:
        sys.stderr.write(TextColor.PURPLE + "This contig's npz file is empty.\n" + TextColor.END)
        return contig_name
    return finish
    


def download_action(ncbi_id, homologous_output_dir,contig_name =None ):
    download_start_time = time.time()
    print_system_log('Download closely-related genomes')
    url_list = download.parser_url(ncbi_id)
    sys.stderr.write(TextColor.GREEN + " INFO: " + str(len(url_list)) + " homologous sequence need to download: \n" + TextColor.END)
    db_path = download.download(homologous_output_dir, ncbi_id, url_list,contig_name)       
    download_end_time = time.time()
    download_time = get_elapsed_time_string(download_start_time, download_end_time)
    print_stage_time('Download closely-related genomes time', download_time)
    return db_path



def meta_polish(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta):
    
    # create a directory
    homologous_output_dir = make_output_dir("homologous", output_dir_debug)

    # use mash screen to get closely related
    print_system_log('Select closely-related genomes')
    mash_start_time = time.time()
    mash_file = mash.screen(assembly, sketch_path, threads, homologous_output_dir, mash_threshold, "meta")
    ncbi_id = mash.meta_get_ncbi_id(mash_file, download_contig_nums)
    mash_end_time = time.time()
    mash_time = get_elapsed_time_string(mash_start_time, mash_end_time)
    print_stage_time('Select closely_related genomes time', mash_time)
   
    if len(ncbi_id) < 5:  # Would'nt polish if closely-related genomes less than 5
        sys.stderr.write(TextColor.PURPLE + "Closely-related genomes less than 5, not to polish...\n" + TextColor.END)
        out.append(assembly)
        os.system('cat {} > {}/{}_homopolished.fasta'.format(' '.join(out), output_dir, assembly_name))
        return

    # download homologous
    download_path = download_action(ncbi_id, homologous_output_dir)



    # alignment
    # align_start_time = time.time()
    paf = alignment.align(assembly, minimap_args, threads, download_path, output_dir_debug)

    # if os.stat(paf).st_size == 0:
    #     sys.stderr.write(TextColor.PURPLE + "Minimap2 can't align, not to polish...\n" + TextColor.END)
    #     return

    # align_end_time = time.time()
    # alignment_time = get_elapsed_time_string(align_start_time, align_end_time)
    # print_stage_time('Minimap2 alignment time', alignment_time)


    contig_start_time = time.time()

    # Each contig polish
    for contig in SeqIO.parse(assembly, 'fasta'):
        if '/' in contig.id:
            contig.id = contig.id.replace('/', '_')
        
        contig_name, contig_output_dir = write_for_new_fasta(contig, output_dir_debug)

        paf_name = contig_output_dir + '/' + contig.id + '.paf'
        cut_cmd = "more {} | grep '{}' > {}".format(paf, contig.id, paf_name)
        os.system(cut_cmd)

        #check homopolish and run homopolish
        out.append(check_homopolish(paf_name, contig_name, contig_output_dir, contig, minimap_args, threads, download_path, model_path))

    contig_end_time = time.time()
    contig_time = get_elapsed_time_string(contig_start_time, contig_end_time)
    print_stage_time('All contig time', contig_time)
    return out



def genus_species_polish(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta):
    
    # create a directory
    homologous_output_dir = make_output_dir("homologous", output_dir_debug)
        
    # Download closely related genome by given genus_species
    print_system_log('Select closely-related genomes and download')
    print("Genus: "+genus_species)
    collect_start_time = time.time()
    ncbi_id, url_list = download.parser_genus_species(genus_species, download_contig_nums)

    # download homologous
    download_path = download_action(ncbi_id, homologous_output_dir,assembly)
    

    # Each contig alignment and polish
    for contig in SeqIO.parse(assembly, 'fasta'):
        if '/' in contig.id:
            contig.id = contig.id.replace('/', '_')
            
        contig_name, contig_output_dir = write_for_new_fasta(contig, output_dir_debug)
             
            
        # alignment
        paf = alignment.align(contig_name, minimap_args, threads, download_path, contig_output_dir)
        
        #check homopolish and run homopolish
        out.append(check_homopolish(paf, contig_name, contig_output_dir, contig, minimap_args, threads, download_path, model_path))

    return out



def without_genus(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta):

    # Each contig to mash and alignment and polish
    for contig in SeqIO.parse(assembly, 'fasta'):
        if '/' in contig.id:
            contig.id = contig.id.replace('/', '_')
            
        contig_name, contig_output_dir = write_for_new_fasta(contig, output_dir_debug)

        print_system_log('Select closely-related genomes')
        select_start_time = time.time()
        ncbi_id = mash_select_closely_related(sketch_path, mash_screen, threads, contig_output_dir, mash_threshold,
                       download_contig_nums, contig_name, contig.id)
        select_end_time = time.time()
        select_time = get_elapsed_time_string(select_start_time, select_end_time)
        print_stage_time('Select closely-related genomes', select_time)
        
        if len(ncbi_id) < 5:  # Would'nt polish if closely-related genomes less than 5
            sys.stderr.write(TextColor.PURPLE + "This contig " + contig.id + " closely-related genome is less than 5, not to polish...\n" + TextColor.END)
            out.append(contig_name)
            continue
        else:
            # download homologous
            download_path = download_action(ncbi_id, contig_output_dir,contig_name)
            # alignment
            paf = alignment.align(contig_name, minimap_args, threads, download_path, contig_output_dir)
        
            #check homopolish and run homopolish
            out.append(check_homopolish(paf, contig_name, contig_output_dir, contig, minimap_args, threads, download_path, model_path))
        
    return out
 


def local_DB(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta, local_DB_path):

    # Each contig to alignment and polish
    for contig in SeqIO.parse(assembly, 'fasta'):
        if '/' in contig.id:
            contig.id = contig.id.replace('/', '_')
            

        contig_name, contig_output_dir = write_for_new_fasta(contig, output_dir_debug)

        #file_path = local_DB_path + '*'
        file_path = local_DB_path
        db_path = contig_output_dir + '/All_homologous_sequences.fna.gz'
        os.system('cat {} > {}'.format(file_path, db_path))
        print('')
        
        # alignment
        paf = alignment.align(contig_name, minimap_args, threads, db_path, contig_output_dir)
        
        #check homopolish and run homopolish
        out.append(check_homopolish(paf, contig_name, contig_output_dir, contig, minimap_args, threads, db_path, model_path))

    return out



def polish_genome(mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta, local_DB_path):
    output_dir = FileManager.handle_output_directory(output_dir)

    # create a directory
    output_dir_debug = make_output_dir("debug", output_dir)

    assembly_name = assembly.rsplit('/', 1)[-1]
    assembly_name = assembly_name.split('.')[0]

    total_start_time = time.time()

    out = []
    if local_DB_path: # use local DB
        out = local_DB(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta, local_DB_path)
    
    elif meta:  # metagenome used screen
        out = meta_polish(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta)

    elif genus_species:   # given genus_species
        out = genus_species_polish(out, assembly_name, output_dir_debug,mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta)

    else:  # single genome without given genus_species -> screen or dist
        out = without_genus(out, assembly_name, output_dir_debug, mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta)
        
    os.system('cat {} > {}/{}_homopolished.fasta'.format(' '.join(out), output_dir, assembly_name))
    total_end_time = time.time()
    total_time = get_elapsed_time_string(total_start_time, total_end_time)
    print_stage_time('Total', total_time)

    if debug:
        try:
            shutil.rmtree(output_dir_debug)
        except OSError as e:
            print(e)
        else:
            return True



def make_train_data(mash_screen, assembly, reference, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug):
    output_dir = FileManager.handle_output_directory(output_dir)
    contig_output_dir_debug = make_output_dir("debug", output_dir)
    
    assembly_name = assembly.rsplit('/',1)[-1]
    assembly_name = assembly_name.split('.')[0]
    
    total_start_time = time.time()
    for contig in SeqIO.parse(assembly, 'fasta'):
        timestr = time.strftime("[%Y/%m/%d %H:%M]")
        sys.stderr.write(TextColor.GREEN + str(timestr) +" INFO: RUN-ID: "+ contig.id + "\n" + TextColor.END)
        contig_output_dir = make_output_dir("contig", contig_output_dir_debug, contig.id)

        contig_name = contig_output_dir + '/' + contig.id +'.fasta'
        SeqIO.write(contig, contig_name, "fasta")

        print_system_log('Select closely-related genomes and download')
        collect_start_time = time.time()
        #db_path = mash_select_closely_related(sketch_path, mash_screen, threads, contig_output_dir, mash_threshold, download_contig_nums, contig_name, contig.id)
        ncbi_id = mash_select_closely_related(sketch_path, mash_screen, threads, contig_output_dir, mash_threshold, download_contig_nums, contig_name, contig.id)
        '''
        if len(ncbi_id) < 5:
            sys.stderr.write(TextColor.PURPLE + "This contig " + contig.id + " closely-related genome is less than 5, not to polish...\n" + TextColor.END)
            out.append(contig_name)
            continue
        '''

        collect_end_time = time.time()
        collect_time = get_elapsed_time_string(collect_start_time, collect_end_time)
        #print_stage_time('Select closely-related genomes and download', collect_time)

        print_system_log('Download closely-related genomes')
        url_list = download.parser_url(ncbi_id)
        sys.stderr.write(TextColor.GREEN + " INFO: " + str(len(url_list)) + " homologous sequence need to download: \n" + TextColor.END)
        db_path = download.download(contig_output_dir, ncbi_id, url_list)


        seq_paf = alignment.align(contig_name, minimap_args, threads, db_path, contig_output_dir)
        ref_paf = alignment.align(contig_name, minimap_args, threads, db_path, contig_output_dir, reference)


        if os.stat(seq_paf).st_size != 0 and os.stat(ref_paf).st_size != 0:
            record = SeqIO.read(contig_name, "fasta")
            genome_size = len(record)

            dataframe_path = homologous_retrieval(seq_paf, genome_size, contig_output_dir, contig.id, contig_name, ref_paf)
                                                                                                                                                                                                                                                                                                                                                                                                                        
        else:
            sys.stderr.write(TextColor.PURPLE + contig.id + " minimap2 can't align......\n" + TextColor.END)

        shutil.move(dataframe_path, output_dir+'/'+assembly_name+'.feather')
