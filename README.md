# Introduction
Homopolish is a Nanopore-specific polisher able to generate a high-quality genome (Q40-50) for virus, bacteria, and fungus.

# Installation
Homopolish is recommendated to install and run within a conda environment

	git clone https://github.com/ythuang0522/homopolish.git
	conda env create -f environment.yml
	conda activate homopolish

# Download virus, bacteria, or fungi sketches
Homopolish makes uses of compressed genomes (Mash sketches) for identifying most-similar genomic sequences. The pre-sketches (virus/bacteria/fungi) can be downloaded.

	Virus: wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/virus.msh.gz
	Bacteria: wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz
	Fungi: wget http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/fungi.msh.gz
**unzip:**
```
gunzip FileName.gz
```
    
    
# Quick usage

Homopolish is recommended to run in conjunction with one of the three sketches.

**polish:**
```
python3 homopolish.py polish -a yourgenome.fasta -s bacteria.msh -m R10.3.pkl -o youroutput
```
This command uses the ```R10.3.pkl``` model and the ```bacteria.msh``` sketch.
> Other model option: ```R9.4.pkl```
> Other sketchs option: ```virus.msh```, ```fungi.msh```

**train model:**
```

```


# Options and usage

Run ```python3 homopolish.py polish --help``` to view this program's most commonly used options:
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
