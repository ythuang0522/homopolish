[![DOI](https://zenodo.org/badge/280070417.svg)](https://zenodo.org/badge/latestdoi/280070417)


# Introduction
Homopolish is a Nanopore-specific polisher able to generate a high-quality genome for virus, bacteria, and fungus. Nanopore systematic errors are  corrected by retrieving homologs from closely-related genomes and a trained ML model. In conjunction with Racon or Medaka, the genomes can reach Q40-90 (>99.99%) accuracy on R9.4 flowcells (Guppy 3.x) (see [Reference](#reference))

# Installation
Homopolish is recommendated to install and run within a conda environment

	git clone https://github.com/ythuang0522/homopolish.git
	cd homopolish
	conda env create -f environment.yml
	conda activate homopolish

# Download virus, bacteria, or fungi sketches
Homopolish retrieves homologous sequences by scanning microbial genomes compressed in (Mash) sketches. Three sketches of virus (74Mb), bacteria (719Mb), and fungi (74Mb) can be downloaded from the following addresses using wget or curl.

	Virus: http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/virus.msh.gz
	Bacteria: http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz
	Fungi: http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/fungi.msh.gz

Then unzip the downloaded skeches.

```
gunzip bacteria.msh.gz
```
    
# Quick usage

Homopolish should be run with a pre-trained model (R9.4.pkl or R10.3.pkl) and one sketch (virus, bacteria, or fungi). Note that Homopolish should be run after Racon or Medaka as it focuses on removing systematic indel errors only. For instance, if your Medaka-polished genome (yourgenome.fasta) is bacteria and sequenced by R9.4 flowcell, please type
```
python3 homopolish.py polish -a yourgenome.fasta -s bacteria.msh -m R9.4.pkl -o youroutput
```
You also can set ```-g``` to specify particular genus and species names in [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt) without mash searching.
```
python3 homopolish.py polish -a yourgenome.fasta -g genusname_speciesname -m R9.4.pkl -o youroutput
```
If you wanna use your own local genomes instead of NCBI, specify the path to your local database via ```-l```.
```
python3 homopolish.py polish -a yourgenome.fasta -l path_to_your_genomes.fasta -m R9.4.pkl -o youroutput
```
# Other Options and usage

Run ```python3 homopolish.py polish --help``` to view all the options:
```
usage: homopolish.py polish [-h] -m MODEL_PATH -a ASSEMBLY
                            (-s SKETCH_PATH | -g GENUS) [-t THREADS]
                            [-o OUTPUT_DIR] [--minimap_args MINIMAP_ARGS]
                            [--mash_threshold MASH_THRESHOLD]
                            [--download_contig_nums DOWNLOAD_CONTIG_NUMS] [-d]
                            [--mash_screen]

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL_PATH, --model_path MODEL_PATH
                        [REQUIRED] Path to a trained model (pkl file). Please
                        see our github page to see options.
  -a ASSEMBLY, --assembly ASSEMBLY
                        [REQUIRED] Path to a assembly genome.
  -s SKETCH_PATH, --sketch_path SKETCH_PATH
                        Path to a mash sketch file.
  -g GENUS, --genus GENUS
                        Genus name
  -l LOCAL_DB_PATH, --local_DB_path LOCAL_DB_PATH
                        Path to your local DB (ex: cat closely-related_genomes1.fasta closely-related_genomes2.fasta> DB.fasta)
  -t THREADS, --threads THREADS
                        Number of threads to use. [1]
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to the output directory. [output]
  --minimap_args MINIMAP_ARGS
                        Minimap2 -x argument. [asm5]
  --mash_threshold MASH_THRESHOLD
                        Mash output threshold. [0.95]
  --download_contig_nums DOWNLOAD_CONTIG_NUMS
                        How much contig to download from NCBI. [20]
  -d, --debug           Keep the information of every contig after mash, such
                        as homologous sequences and its identity infomation.
                        [no]
  --mash_screen         Use mash screen. [mash dist]
```
# Output Files

**Ordinary output:** 

All program's output files will be saved in the folder named youroutput(```-o youroutput``` ), and you will only get one file named ```yourgenome_homopolished.fasta```.

**Debugging mode:** 

If you use the parameter ```-d```, directory content in a tree-like format is below.
* ```homologous_sequences``` contains other homologous species
*  ```All_homologous_sequences.fna.gz``` which concatenate all file in ```homologous_sequences```

```
├── yourgenome_homopolished.fasta
└── debug
    ├── contig_1_segment0
    │   ├── All_homologous_sequences.fna.gz 
    │   ├── contig_1_segment0.fasta
    │   ├── contig_1_segment0.feather
    │   ├── contig_1_segment0.npz
    │   ├── contig_1_segment0.paf
    │   ├── contig_1_segment0.sort.tab
    │   ├── homologous_sequences
    │   │   ├── GCF_000775955.1_ASM77595v1_genomic.fna.gz
    │   │   └── .......
    │   ├── polished.fasta
    │   └── result.feather
    ├── contig_2_segment0
    │   └── ......
    └── ......
```

# Reference

Comparison of genome accuracy polished by Racon, Medaka, MarginPolish, HELEN, and Homopolish on R9.4 and on R10.3 Zymo Microbial Community Standard. See also the [data](https://github.com/ythuang0522/homopolish/tree/master/data) folder. Median Q scores computed by [fastmer](https://github.com/jts/assembly_accuracy/blob/master/fastmer.py). 
![Accuracy of Homopolish](https://www.biorxiv.org/content/biorxiv/early/2020/09/20/2020.09.19.304949/F1.large.jpg)
![Accuracy of Homopolish](https://www.biorxiv.org/content/biorxiv/early/2020/09/20/2020.09.19.304949/F2.large.jpg)

If you use homopolish, please cite

Huang, Y.-T., Liu, P.-Y., and Shih, P.-W. [Homopolish: a method for the revmoal of systematic errors in nanopore sequencing by homologous polishing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02282-6), Genome Biology, 2021.




