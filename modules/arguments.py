def add_common_arguments(parser):
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
        help="Mash output threshold. [0.95]"
    )
    parser.add_argument(
        "--download_contig_nums",
        type=str,
        required=False,
        default='20',
        help="How much contig to download from NCBI. [20]"
    )
    parser.add_argument(
        "-d",
        "--debug",
        required=False,
        action = "store_false",
        help="Keep the information of every contig after mash, such as homologous sequences and its identity infomation. [no]"
    )
    parser.add_argument(
        "--mash_screen",
        required=False,
        action='store_true',
        default=False,
        help="Use mash screen. [mash dist]"
    )

    return parser
def add_polish_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-m",
        "--model_path",
        type=str,
        required=True,
        help="[REQUIRED] Path to a trained model (pkl file). Please see our github page to see options."
    )    
def add_train_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-d",
        "--dataframe_dir",
        type=str,
        required=True,
        help="[REQUIRED] Path to a directory for alignment dataframe."
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
        "-p",
        "--output_prefix",
        type=str,
        required=False,
        default='train',
        help="Prefix for the train model. [train]"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        default=1,
        help="Number of threads to use. [1]"
    )
    return parser

def add_train_data_arguments(parser):
    """
    Add arguments to a parser for sub-command "polish"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        required=True,
        help="[REQUIRED] True reference aligned to assembly genome. Include labels in output."
    )
    return parser