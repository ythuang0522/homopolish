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

def select_closely_related(sketch_path, genus_species, mash_screen, threads, output_dir, mash_threshold,
                           download_contig_nums, contig_name=None, contig_id=None):
    if sketch_path:
        if mash_screen:
            mash_file = mash.screen(contig_name, sketch_path, threads, output_dir, mash_threshold, download_contig_nums,
                                    contig_id)
        else:
            mash_file = mash.dist(contig_name, sketch_path, threads, output_dir, mash_threshold, download_contig_nums,
                                  contig_id)
        ncbi_id = mash.get_ncbi_id(mash_file)
        if len(ncbi_id) < 5:  # Would'nt polish if closely-related genomes less than 5
            return False
        url_list = download.parser_url(ncbi_id)

    if genus_species:
        ncbi_id, url_list = download.parser_genus_species(genus_species, download_contig_nums)


    db_path = download.download(output_dir, ncbi_id, url_list)

    return db_path

def homologous_retrieval(assembly, minimap_args, threads, sequence, output_dir, reference=None):

    seq_npz = alignment.align(assembly, minimap_args, threads, sequence, output_dir)
    if reference:
        ref_npz = alignment.align(assembly, minimap_args, threads, reference, output_dir)
        df = align2df.todf(assembly, seq_npz, output_dir, ref_npz)
    else:
        print(seq_npz)
        df = align2df.todf(assembly, seq_npz, output_dir)
    return df

def polish_genome(mash_screen, assembly, model_path, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug, meta):
    output_dir = FileManager.handle_output_directory(output_dir)
    contig_output_dir_debug = output_dir + '/debug'
    contig_output_dir_debug = FileManager.handle_output_directory(contig_output_dir_debug)
    assembly_name = assembly.rsplit('/', 1)[-1]
    assembly_name = assembly_name.split('.')[0]

    if meta or genus_species:
        meta_output_dir = contig_output_dir_debug + '/homologous'
        meta_output_dir = FileManager.handle_output_directory(meta_output_dir)
        collect_start_time = time.time()
        db_path = select_closely_related(sketch_path, genus_species, mash_screen, threads, meta_output_dir, mash_threshold, download_contig_nums)
        collect_end_time = time.time()

        out = []
        total_start_time = time.time()
        for contig in SeqIO.parse(assembly, 'fasta'):
            timestr = time.strftime("[%Y/%m/%d %H:%M]")
            sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO: RUN-ID: " + contig.id + "\n" + TextColor.END)
            contig_output_dir = contig_output_dir_debug + '/' + contig.id
            contig_output_dir = FileManager.handle_output_directory(contig_output_dir)
            contig_name = contig_output_dir + '/' + contig.id + '.fasta'
            SeqIO.write(contig, contig_name, "fasta")

            print_system_log('Homologous retrieval')
            homologous_start_time = time.time()
            dataframe = homologous_retrieval(contig_name, minimap_args, threads, db_path, contig_output_dir)
            homologous_end_time = time.time()

            print_system_log('Prediction')
            predict_start_time = time.time()
            result = prediction.predict(dataframe, model_path, threads, contig_output_dir)
            predict_end_time = time.time()

            print_system_log('Polish')
            polish_start_time = time.time()
            finish = polish.stitch(contig_name, result, contig_output_dir)
            polish_end_time = time.time()

            # calculating time
            homologous_time = get_elapsed_time_string(homologous_start_time, homologous_end_time)
            predict_time = get_elapsed_time_string(predict_start_time, predict_end_time)
            polish_time = get_elapsed_time_string(polish_start_time, polish_end_time)

            # print stage time

            print_stage_time('Homologous retrieval', homologous_time)
            print_stage_time('Prediction', predict_time)
            print_stage_time('Polish', polish_time)
            out.append(finish)

        collect_time = get_elapsed_time_string(collect_start_time, collect_end_time)
        print_stage_time('Select closely-related genomes', collect_time)
    else:
        out = []
        total_start_time = time.time()
        for contig in SeqIO.parse(assembly, 'fasta'):
            timestr = time.strftime("[%Y/%m/%d %H:%M]")
            sys.stderr.write(TextColor.GREEN + str(timestr) + " INFO: RUN-ID: " + contig.id + "\n" + TextColor.END)
            contig_output_dir = contig_output_dir_debug + '/' + contig.id
            contig_output_dir = FileManager.handle_output_directory(contig_output_dir)
            contig_name = contig_output_dir + '/' + contig.id + '.fasta'
            SeqIO.write(contig, contig_name, "fasta")

            print_system_log('Select closely-related genomes')
            collect_start_time = time.time()
            db_path = select_closely_related(sketch_path, genus_species, mash_screen, threads, contig_output_dir,
                                             mash_threshold, download_contig_nums, contig_name, contig.id)
            collect_end_time = time.time()
            if db_path == False:  # contig which didn't polish
                out.append(contig_name)
                continue

            print_system_log('Homologous retrieval')
            homologous_start_time = time.time()
            dataframe = homologous_retrieval(contig_name, minimap_args, threads, db_path, contig_output_dir)
            homologous_end_time = time.time()

            print_system_log('Prediction')
            predict_start_time = time.time()
            result = prediction.predict(dataframe, model_path, threads, contig_output_dir)
            predict_end_time = time.time()

            print_system_log('Polish')
            polish_start_time = time.time()
            finish = polish.stitch(contig_name, result, contig_output_dir)
            polish_end_time = time.time()

            # calculating time
            collect_time = get_elapsed_time_string(collect_start_time, collect_end_time)
            homologous_time = get_elapsed_time_string(homologous_start_time, homologous_end_time)
            predict_time = get_elapsed_time_string(predict_start_time, predict_end_time)
            polish_time = get_elapsed_time_string(polish_start_time, polish_end_time)

            # print stage time
            print_stage_time('Select closely-related genomes', collect_time)
            print_stage_time('Homologous retrieval', homologous_time)
            print_stage_time('Prediction', predict_time)
            print_stage_time('Polish', polish_time)
            out.append(finish)

    os.system('cat {} > {}/{}_homopolished.fasta'.format(' '.join(out), output_dir, assembly_name))

    total_end_time = time.time()
    total_time = get_elapsed_time_string(total_start_time, total_end_time)
    print_stage_time('Total', total_time)

    if debug:
        try:
            shutil.rmtree(contig_output_dir_debug)
        except OSError as e:
            print(e)
        else:
            return True

def make_train_data(mash_screen, assembly, reference, sketch_path, genus_species, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug):
    output_dir = FileManager.handle_output_directory(output_dir)
    contig_output_dir_debug = output_dir + '/debug'
    contig_output_dir_debug = FileManager.handle_output_directory(contig_output_dir_debug)
    assembly_name = assembly.rsplit('/',1)[-1]
    assembly_name = assembly_name.split('.')[0]
    
    total_start_time = time.time()
    for contig in SeqIO.parse(assembly, 'fasta'):
        timestr = time.strftime("[%Y/%m/%d %H:%M]")
        sys.stderr.write(TextColor.GREEN + str(timestr) +" INFO: RUN-ID: "+ contig.id + "\n" + TextColor.END)
        contig_output_dir = contig_output_dir_debug + '/' + contig.id
        contig_output_dir = FileManager.handle_output_directory(contig_output_dir)
        contig_name = contig_output_dir + '/' + contig.id +'.fasta'
        SeqIO.write(contig, contig_name, "fasta")

        print_system_log('Select closely-related genomes')
        collect_start_time = time.time()
        db_path = select_closely_related(sketch_path, genus_species, mash_screen, threads, contig_output_dir, mash_threshold, download_contig_nums, contig_name, contig.id)
        collect_end_time = time.time()

        print_system_log('Homologous retrieval')
        homologous_start_time = time.time()
        dataframe_path = homologous_retrieval(contig_name, minimap_args, threads, db_path, contig_output_dir, reference)
        homologous_end_time = time.time()
        
        #calculating time
        collect_time = get_elapsed_time_string(collect_start_time, collect_end_time)
        homologous_time = get_elapsed_time_string(homologous_start_time, homologous_end_time)
        
        #print stage time       
        print_stage_time('Select closely-related genomes', collect_time)
        print_stage_time('Homologous retrieval', homologous_time)
        shutil.move(dataframe_path, output_dir)

    total_end_time = time.time()
    total_time = get_elapsed_time_string(total_start_time, total_end_time)
    print_stage_time('Total', total_time)

    if debug:
        try:
            shutil.rmtree(contig_output_dir_debug)
        except OSError as e:
            print(e)
        else:
            return True