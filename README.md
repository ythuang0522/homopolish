[![DOI](https://zenodo.org/badge/280070417.svg)](https://zenodo.org/badge/latestdoi/280070417)
[![DOI](https://img.shields.io/badge/Download-Bioconda-green)](https://anaconda.org/bioconda/homopolish)


## Introduction
Homopolish is a genome polisher originally developed for Nanopore and subsequently extended for PacBio CLR. It generates a high-quality genome (>Q50) for virus, bacteria, and fungus. Nanopore/PacBio systematic errors are corrected by retreiving homologs from closely-related genomes and polished by an SVM. When paired with Racon and Medaka, the genome quality can reach Q50-90 (>99.999%) on Nanopore R9.4/10.3 flowcells (Guppy >3.4). For PacBio CLR, Homopolish also improves the majority of Flye-assembled genomes to Q90 (see [Accuracy](#accuracy)).

---
- [Installation](#installation)
- [Quick usage for polishing ONT homopolymer errors](#quick-usage-for-polishing-homopolymer-errors)
- [Quick usage for polishing ONT modification-mediated errors](#quick-usage-for-polishing-modification-errors)
- [Accuracy](#accuracy)
- [Citation](#citation)
- [Contact](#contact)
---
## Installation
Homopolish is recommendated to install and run within a conda environment

	git clone https://github.com/ythuang0522/homopolish.git
	cd homopolish
	conda env create -f environment.yml
	conda activate homopolish

Homopolish is available in bioconda but not guaranteed to be the latest version. An executable file, homopolish, will be found after conda installation. 

	conda create -n homopolish -c conda-forge -c bioconda homopolish

## Download virus, bacteria, or fungi sketches
Homopolish retrieves homologous sequences by scanning microbial genomes compressed in (Mash) sketches. Three sketches of virus (74Mb), bacteria (719Mb), and fungi (74Mb) can be downloaded from the following addresses using wget or curl.

	Virus: http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/virus.msh.gz
	Bacteria: http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/bacteria.msh.gz
	Fungi: http://bioinfo.cs.ccu.edu.tw/bioinfo/mash_sketches/fungi.msh.gz

Then unzip the downloaded skeches.

```
gunzip bacteria.msh.gz
```
    
## Quick usage for polishing homopolymer errors

Homopolish should be run with a pre-trained model (R9.4.pkl/R10.3.pkl for Nanopore and pb.pkl for PacBio CLR) and one sketch (virus, bacteria, or fungi). For Nanopore sequencing, Homopolish should be run after the Racon-Medaka pipeline as it only removes indel errors. For PacBio CLR sequencing, it can be invoked directly after Flye assembly. For instance, if your Medaka-polished genome (yourgenome.fasta) is bacteria and sequenced by R9.4 flowcell, please type
```
python3 homopolish.py polish -a yourgenome.fasta -s bacteria.msh -m R9.4.pkl -o youroutput
```

If installed from bioconda, you should directly invoke the script.
```
homopolish polish -a yourgenome.fasta -s bacteria.msh -m R9.4.pkl -o youroutput
```

You can specify particular genus and species via ```-g``` (Taxanomic names in [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt)) without mash searching.
```
python3 homopolish.py polish -a yourgenome.fasta -g genusname_speciesname -m R9.4.pkl -o youroutput
```

If you wanna use private local genomes instead of NCBI, specify the path to your own (fasta) database via ```-l```.
```
python3 homopolish.py polish -a yourgenome.fasta -l path_to_your_genomes.fasta -m R9.4.pkl -o youroutput
```

### Other Options and usage

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

### Output File
The program's output files will be saved in the folder named youroutput (```-o youroutput``` ), and the polished genoem is named ```yourgenome_homopolished.fasta```.

## Quick usage for polishing modification errors
Homopolish impleneted a submodule called modpolish for correcting modification-mediated errors. Given a draft genome with ONT reads (fastq), modpolish correcting the modification-mediated errors using reads, quality and homologs.
```
python3 homopolish.py modpolish -a yourgenome.fasta -q fastq_path -m R9.4.pkl -s bacteria.msh -o youroutput
```
You can supply the bam (i.e., reads to draft genome) insetad of reads for skipping the time-consuming alignment.
```
python3 homopolish.py modpolish -a yourgenome.fasta -b bam_path -m R9.4.pkl -s bacteria.msh -o youroutput
```

Run ```python3 homopolish.py modpolish --help``` to view all the options:
```
usage: python3 homopolish.py modpolish -a dragt genome file path -q fastq file path -s SKETCH_PATH
                            (-o output modified_position_csv)                          
                        
optional arguments:
  -h, --help            show this help message and exit

  -a Draft Genome File Path --Draft genome file path
                        [REQUIRED] Path to a assembly genome.
  -s SKETCH_PATH, --sketch_path SKETCH_PATH
                        Path to a mash sketch file.
                        Genus name
                        Path to your local DB (ex: cat closely-related_genomes1.fasta closely-related_genomes2.fasta> DB.fasta)
  -o                    Boolean,output mod_polish fix position csv
  -q                    fastq File
```
  
### Output File of Modpolish
A genome named ```yourgenome_modpolished.fasta``` will be outputted.

## Debugging mode

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

## Accuracy

Comparison of genome accuracy polished by Racon, Medaka, MarginPolish, HELEN, and Homopolish on Nanopore R9.4. Median Q scores were computed by [fastmer](https://github.com/jts/assembly_accuracy/blob/master/fastmer.py). We note that these are based on early ONT basecaller (Guppy 3.2). After Guppy 3.4, we have seen consistent >Q50 genomes.
![Accuracy of Homopolish](https://www.biorxiv.org/content/biorxiv/early/2020/09/20/2020.09.19.304949/F1.large.jpg)
![Accuracy of Homopolish](https://www.biorxiv.org/content/biorxiv/early/2020/09/20/2020.09.19.304949/F2.large.jpg)

Since v0.3, we found [FastANI](https://github.com/ParBLiSS/FastANI) is more accurate than [mash](https://github.com/marbl/Mash) for selecting highly-similar genomes for polishing. The FastANI version improved the accuracy of R9.4 and R10.3 (see below Guppy 3.2 results). We note that short contigs/plasmids (<~6kb) won't benefit from FastANI as it can't output distance estimatino.
![R9.4 of Homopolish](https://github.com/ythuang0522/homopolish/blob/master/images/R9.4.png)
![R10.3 of Homopolish](https://github.com/ythuang0522/homopolish/blob/master/images/R10.3.png)

With the addition of FastANI, genomes sequenced by PacBio CLR (and assembled by Flye) can now be also siginificantly improved by Homopolish.
![PacBioCLR of Homopolish](https://github.com/ythuang0522/homopolish/blob/master/images/PacBio%20CLR.jpg)

## Citation
If you use homopolish, please cite

Huang, Y.-T., Liu, P.-Y., and Shih, P.-W. [Homopolish: a method for the revmoal of systematic errors in nanopore sequencing by homologous polishing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02282-6), Genome Biology, 2021.

## Contact
ythuang at cs.ccu.edu.tw
