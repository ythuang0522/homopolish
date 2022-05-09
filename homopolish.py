import argparse
from version import __version__
from modules.arguments import *
from modules.polish_interface import polish_genome
from modules.polish_interface import make_train_data
from modules.train_interface import train_model
import os
from os import path
from modules.VAtypeClass import FixSNP
from modules.mod_polish import getPos


def main():
    parser = argparse.ArgumentParser(description="Homopolish fixes systematic errors in ONT-based assemblies. \n",
#                                                 "1) polish: Run the polishing pipeline.\n"
#                                                 "2) train: Train your own SVM model.\n"
#                                                "3) make_train_data: Make training data with reference genome.",
                                                     
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v",
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )
    subparsers = parser.add_subparsers(dest='sub_command')
    parser_polish = subparsers.add_parser('polish', help="Polish your genomes.")
    add_polish_arguments(parser_polish)
    add_common_arguments(parser_polish)
    
    parser_train = subparsers.add_parser('train', help="Train your model.")
    add_train_arguments(parser_train)


    parser_make_train_data = subparsers.add_parser('make_train_data', help="Prepare training data with truth genome.") 
    add_train_data_arguments(parser_make_train_data)    
    add_common_arguments(parser_make_train_data)

    
    mod_polish = subparsers.add_parser('mod_polish',help = "polish genome by reads")
    add_modpolish_arguments(mod_polish)
    
 
    FLAGS, unparsed = parser.parse_known_args()

    if FLAGS.sub_command == 'polish':
        this_directory = path.abspath(path.dirname(__file__))
        __pkg_path__ = os.path.join(this_directory,FLAGS.model_path)
        polish_genome(FLAGS.mash_screen, FLAGS.assembly, __pkg_path__, FLAGS.sketch_path, FLAGS.genus, FLAGS.threads, \
                FLAGS.output_dir, FLAGS.minimap_args, FLAGS.mash_threshold, FLAGS.download_contig_nums, FLAGS.debug, FLAGS.meta, FLAGS.local_DB_path,FLAGS.coverage,FLAGS.distance)

    elif FLAGS.sub_command == 'train':
        train_model(FLAGS.dataframe_dir, FLAGS.output_dir, FLAGS.output_prefix, FLAGS.threads,FLAGS.pacbio)

    elif FLAGS.sub_command == 'make_train_data':
        make_train_data(FLAGS.mash_screen, FLAGS.assembly, FLAGS.reference, FLAGS.sketch_path, FLAGS.genus, FLAGS.threads, \
                FLAGS.output_dir, FLAGS.minimap_args, FLAGS.mash_threshold, FLAGS.download_contig_nums, FLAGS.debug,FLAGS.coverage,FLAGS.distance)
    
    
    elif FLAGS.sub_command == 'mod_polish':
        if(FLAGS.fastq == "" and FLAGS.bam == ""):
            print("need fastq or bam file!")
            return
        
        fixData = FixSNP()
        fixData.draft_genome_file = FLAGS.fasta
        fixData.reads_file = FLAGS.fastq
        fixData.get_fixCSV_Flag = FLAGS.outFixCSV 
        fixData.spPattern = FLAGS.pattern
        fixData.bamFile = FLAGS.bam    
        getPos(fixData)
    
    elif FLAGS.sub_command == 'mod_polish_posData':
        if(FLAGS.fastq == "" and FLAGS.bam == ""):
            print("need fastq or bam file!")
            return
        fixData = FixSNP()
        fixData.draft_genome_file = FLAGS.fasta
        fixData.reads_file = FLAGS.fastq
        fixData.get_fixCSV_Flag = FLAGS.outFixCSV 
        fixData.get_MissCSV_Flag = FLAGS.outMissCSV
        fixData.get_EorCSV_Flag = FLAGS.outErrorCSV
        fixData.spPattern = FLAGS.pattern
        fixData.bamFile = FLAGS.bam    
        getPos(fixData)
    
    
    elif FLAGS.version is True:
        print("Homopolish VERSION: ", __version__)

    
    else:
        parser.print_help()
    
if __name__ == "__main__":
    main()
