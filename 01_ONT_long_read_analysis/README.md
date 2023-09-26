# Long-read RNA seq analysis using StringTie2

This is a `snakemake` pipeline that takes Oxford Nanopore Sequencing (ONT) data (fastq) as input, generates fastq stats using `nanostat`, performs fastq processing and filtering using `pychopper`, map the reads to the genome using `minimap2` and uses `StringTie2` to assemble and quantify transcripts. 

# Getting Started

## Input

- ONT fastq reads
- Reference genome assembly in fasta format
- GTF: [Gencode GTF](https://www.gencodegenes.org/human/)

### Data
We used data from [Aguzzoli Heberle et al., 2023](https://www.biorxiv.org/content/10.1101/2023.08.06.552162v1.full). The data was downloaded from synapse: https://www.synapse.org/#!Synapse:syn52047893/files/. This is long-read RNA seq of frontal cortex from 6 AD samples and 6 controls. The sequencing was done on ONT promethION using one flow cell per sample.

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- The rest of the dependencies (including snakemake) are installed via conda through the `environment.yml` file


## Installation

Clone the directory:

```bash
git clone --recursive https://github.com/egustavsson/ONT_TALON.git
```

Create conda environment for the pipeline which will install all the dependencies:

```bash
cd ONT_TALON
conda env create -f environment.yml
```

## Usage

Edit `config.yml` to set up the working directory and input files/directories. `snakemake` command should be issued from within the pipeline directory. Please note that before you run any of the `snakemake` commands, make sure to first activate the conda environment using the command `conda activate ont_talon`.

```bash
cd ONT_TALON
conda activate ont_talon
snakemake --use-conda -j <num_cores> all
```
It is a good idea to do a dry run (using -n parameter) to view what would be done by the pipeline before executing the pipeline.

```bash
snakemake --use-conda -n all
```

You can visualise the processes to be executed in a DAG:

```bash
snakemake --dag | dot -Tpng > dag.png
```

To exit a running `snakemake` pipeline, hit `ctrl+c` on the terminal. If the pipeline is running in the background, you can send a `TERM` signal which will stop the scheduling of new jobs and wait for all running jobs to be finished.

```bash
killall -TERM snakemake
```

To deactivate the conda environment:
```bash
conda deactivate
```

## Output
```
working directory  
|--- config.yml           # a copy of the parameters used in the pipeline  
|--- Nanostat/  
     |-- # output of nanostat - fastq stats  
|--- Pychopper/  
     |-- # output of pychopper - filtered fastq  
|--- Mapping/  
     |-- # output of minimap2 - aligned reads  
|--- Talon/  
     |-- # output of Talon  
     |-- _talon.gtf                       # assembled transcripts  
     |-- _talon_abundance_filtered.tsv    # transcript abundance  
     
```
