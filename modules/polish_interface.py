import os
import sys
import time
import shutil
from Bio import SeqIO
from modules import mash
from modules import download
from modules import alignment
from modules import align2df
from modules import predict
from modules import polish
from modules.TextColor import TextColor
from modules.FileManager import FileManager


def print_system_log(stage):
    timestr = time.strftime("[%Y/%m/%d %H:%M]")
    sys.stderr.write(TextColor.GREEN + str(timestr) +" INFO: Stage: "+ stage + "\n" + TextColor.END)

def print_stage_time(stage, time):
    sys.stderr.write(TextColor.GREEN + "INFO: "+ stage + " TIME: " + str(time) + "\n" + TextColor.END)
        
def get_elapsed_time_string(start_time, end_time):
    """
    Get a string representing the elapsed time given a start and end time.
    :param start_time: Start time (time.time())
    :param end_time: End time (time.time())
    :return:
    """
    elapsed = end_time - start_time
    hours = int(elapsed / 60**2)
    mins = int(elapsed % 60**2 / 60)
    secs = int(elapsed % 60**2 % 60)
    time_string = "{} HOURS {} MINS {} SECS.".format(hours, mins, secs)

    return time_string

def polish_genome(FLAGS,assembly, model_path, sketch_path, genus, threads, output_dir, minimap_args, mash_threshold, download_contig_nums, debug):    
    
    out = []
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
        
        if sketch_path:
            mash_start_time = time.time()
            if FLAGS.mash_screen:
                print_system_log('MASH SCREEN')
                mash_file = mash.screen(contig_name, sketch_path, threads, contig_output_dir, mash_threshold, download_contig_nums, contig.id)                
            else :
                print_system_log('MASH DIST')
                mash_file = mash.dist(contig_name, sketch_path, threads, contig_output_dir, mash_threshold , download_contig_nums, contig.id)

            mash_end_time = time.time()
            ncbi_id = mash.get_ncbi_id(mash_file)  
            if len(ncbi_id) < 5: #Would'nt polish if closely-related genomes less than 5
                out.append(contig_name)
                continue        
            url_list = download.parser_url(ncbi_id)

        if genus:
            ncbi_id, url_list = download.parser_genus(genus)      
        
        download_start_time = time.time()
        print_system_log('DOWNLOAD CONTIGS')              
        db = download.download(contig_output_dir, ncbi_id, url_list)
        download_end_time = time.time()           


        pileup_start_time = time.time()
        print("\n")
        print_system_log('PILE UP')
        db_npz = alignment.align(contig_name, minimap_args, threads, db, contig_output_dir)
        if db_npz == False:
            continue
        pileup_end_time = time.time()            


        align2df_start_time = time.time()           
        print_system_log('TO DATAFRAME')
        df = align2df.todf(contig_name, db_npz, contig_output_dir)
        align2df_end_time = time.time()
        
    
        predict_start_time = time.time()
        print_system_log('PREDICT')
        df = contig_output_dir + '/' + contig.id + '.feather'
        result = predict.predict(df, model_path, threads, contig_output_dir)
        predict_end_time = time.time()


        polish_start_time = time.time()
        print_system_log('POLISH')
        finish = polish.stitch(contig_name, result, contig_output_dir)
        polish_end_time = time.time()

        #calculating time
        if sketch_path:
            mash_time = get_elapsed_time_string(mash_start_time, mash_end_time)
            if FLAGS.mash_screen:
                print_stage_time('MASH SCREEN', mash_time)
            else:
                print_stage_time('MASH DIST', mash_time)
                
        download_time = get_elapsed_time_string(download_start_time, download_end_time)
        pileup_time = get_elapsed_time_string(pileup_start_time, pileup_end_time)
        align2df_time = get_elapsed_time_string(align2df_start_time, align2df_end_time)
        predict_time = get_elapsed_time_string(predict_start_time, predict_end_time)
        polish_time = get_elapsed_time_string(polish_start_time, polish_end_time)
        
        #print stage time       
        print_stage_time('DOWNLOAD', download_time)
        print_stage_time('PILEUP', pileup_time)
        print_stage_time('TO DATAFRAME', align2df_time)
        print_stage_time('PREDICT', predict_time)
        print_stage_time('POLISH', polish_time)
        out.append(finish)

    
    os.system('cat {} > {}/{}_homopolished.fasta'.format(' '.join(out), output_dir, assembly_name))

    if debug:
        try:
            shutil.rmtree(contig_output_dir_debug)
        except OSError as e:
            print(e)
        else:
            return True

    total_end_time = time.time()
    total_time = get_elapsed_time_string(total_start_time, total_end_time)
    print_stage_time('Total', total_time)