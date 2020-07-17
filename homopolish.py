import argparse
from version import __version__
from modules.polish_interface import polish_genome

def add_polish_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-a",
        "--assembly",
        type=str,
        required=True,
        help="[REQUIRED] Path to a assembly genome."
    )
    parser.add_argument(
        "-m",
        "--model_path",
        type=str,
        required=True,
        help="[REQUIRED] Path to a trained model (pkl file). Please see our github page to see options."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-s",
        "--sketch_path",
        type=str,
        required=False,
        help="Path to a mash sketch file."
    )
    group.add_argument(
        "-g",
        "--genus",
        type=str,
        required=False,
        help="Genus name"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use. [1]"
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=False,
        default='./output/',
        help="Path to the output directory. [output]"
    )
    parser.add_argument(
        "--minimap_args",
        type=str,
        required=False,
        default='asm5',
        help="Minimap2 -x argument. [asm5]"
    )
    parser.add_argument(
        "--mash_threshold",
        type=str,
        required=False,
        default='0.95',
        help="Mash screen output threshold. [0.95]"
    )
    parser.add_argument(
        "--download_contig_nums",
        type=str,
        required=False,
        default='20',
        help="How much contig to download from NCBI. [20]"
    )
    

    return parser

def main():
    parser = argparse.ArgumentParser(description="Homopolish is a SVM based polisher for polishing ONT-based assemblies. \n"
                                                 "1) polish: Run the polishing pipeline.\n"
                                                 "2) train: still in development",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v",
        "--version",
        default=False,
        action='store_true',
        help="Show version."
    )
    subparsers = parser.add_subparsers(dest='sub_command')
    parser_polish = subparsers.add_parser('polish', help="Run the polishing pipeline.")
    add_polish_arguments(parser_polish)

    FLAGS, unparsed = parser.parse_known_args()
    if FLAGS.sub_command == 'polish':
        
        polish_genome(FLAGS.assembly, FLAGS.model_path, FLAGS.sketch_path, FLAGS.genus, FLAGS.threads, \
                FLAGS.output_dir, FLAGS.minimap_args, FLAGS.mash_threshold, FLAGS.download_contig_nums)

    elif FLAGS.version is True:
        print("Homopolish VERSION: ", __version__)

    
if __name__ == "__main__":
    main()