# Introduction
Homopolish is a Nanopore-specific polisher able to generate a high-quality genome (Q40-90, aka >99.99% accuracy) for virus, bacteria, and fungus. The strength of Homopolish lies in correcting Nanopore systematic errors by training a polishing model using homologous sequences retreived from closely-related genomes.

# Installation
Homopolish is recommendated to install and run within a conda environment

	git clone https://github.com/ythuang0522/homopolish.git
	conda env create -f environment.yml
	conda activate homopolish

# Download virus, bacteria, or fungi sketches
Homopolish retrieves homologous sequences by scanning microbial genomes compressed in (Mash) sketches. Three sketches of virus, bacteria, and fungi can be downloaded from the following addresses.

	Virus: wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/virus.msh.gz
	Bacteria: wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz
	Fungi: wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/fungi.msh.gz

Then unzip the downloaded skeches.

```
gunzip bacteria.msh.gz
```
    
# Quick usage

Homopolish should be run with a pre-trained model (R9.4.pkl or R10.3.pkl) and one sketch (virus, bacteria, or fungi). We note that Homopolish achieves best results when run after Racon or Medaka as it focuses on removing systematic indel errors only. For instance, if your Medaka-polished genome (yourgenome.fasta) is bacteria and sequenced by R9.4 flowcell, please type
```
python3 homopolish.py polish -a yourgenome.fasta -s bacteria.msh -m R9.4.pkl -o youroutput
```
You also can use ```-g``` to set the given Genus name to polish your genome.
```
python3 homopolish.py polish -a yourgenome.fasta -s bacteria.msh -g genusname -m R9.4.pkl -o youroutput
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

**Debuging mode:** 

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

