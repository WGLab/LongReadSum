# LongReadSum: A fast and flexible QC tool for long read sequencing data

![image](https://github.com/user-attachments/assets/ae3d2b6a-46da-476a-8fe5-7e0a34751585)

[![build tests](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml/badge.svg)](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml)

LongReadSum supports FASTA, FASTQ, BAM, FAST5, and sequencing_summary.txt file formats for quick generation of QC data in HTML and text format.

# README Contents
- [Installation using Anaconda](#installation-using-anaconda)
- [Installation using Docker](#installation-using-anaconda)
- [Building from source](#building-from-source)
- General usage for common filetypes:
  - [Common parameters](#common-parameters)
  - [WGS BAM](#wgs-bam)
  - [RRMS BAM](#rrms-bam)
  - [RNA-Seq BAM (TIN values)](#rna-seq-bam)
  - [ONT POD5](#ont-pod5)
  - [ONT FAST5](#ont-fast5)
  - [Basecall summary (ONT sequencing_summary.txt)](#basecall-summary)
  - [FASTQ](#fastq)
  - [FASTA](#fasta)
- [Revision history](#revision-history)
- [Getting help](#getting-help)
- [Citing LongReadSum](#citing-longreadsum)

# Installation using Anaconda
First, install [Anaconda](https://www.anaconda.com/).

Next, create a new environment. This installation has been tested with Python 3.10, Linux 64-bit.

```
conda create -n longreadsum python=3.10
conda activate longreadsum
```

LongReadSum and its dependencies can then be installed using the following command:

```
conda install -c bioconda -c wglab longreadsum=1.3.1
```

# Installation using Docker
First, install [Docker](https://docs.docker.com/engine/install/).
Pull the latest image from Docker hub, which contains the latest longreadsum release and its dependencies.

```
docker pull genomicslab/longreadsum
```

## Running

On Unix/Linux:
```
docker run -v C:/Users/.../DataDirectory:/mnt/ -it genomicslab/longreadsum bam -i /mnt/input.bam -o /mnt/output
```
Note that the `-v` command is required for Docker to find the input file. Use a directory under `C:/Users/` to ensure volume files are mounted correctly. In the above example, the local directory `C:/Users/.../DataDirectory` containing the input file `input.bam` is mapped to a directory `/mnt/` in the Docker container. Thus, the input file and output directory arguments are relative to the `/mnt/` directory, but the output files will also be saved locally in `C:/Users/.../DataDirectory` under the specified subdirectory `output`.


# Building from source
To get the latest updates in longreadsum, you can build from source.
First install [Anaconda](https://www.anaconda.com/). Then follow the instructions below to install LongReadSum and its dependencies:

```
# Pull the latest updates
git clone https://github.com/WGLab/LongReadSum
cd LongReadSum

# Create the longreadsum environment, install dependencies, and activate
conda env create -f environment.yml
conda activate longreadsum

# Build the program
make
```

If you are using FAST5 files with VBZ compression, you will need to download and install the VBZ plugin corresponding to your architecture:
https://github.com/nanoporetech/vbz_compression/releases

```
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar -xf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
```

Finally, add the plugin to your path:
```
export HDF5_PLUGIN_PATH=/full/path/to/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
```


## Running
Activate the conda environment and then run with arguments:
```
conda activate longreadsum
python longreadsum [arguments]
```

# General Usage
## Common parameters

To see all parameters for a filetype, run:

```longreadsum [filetype] --help```

This section describes parameters common to all filetypes:

| Parameter	| Description | Default |
| --- | --- | --- |
| -i, --input | A single input filepath
| -I, --inputs | Multiple comma-separated input filepaths
| -P, --pattern | Use pattern matching (*) to specify multiple input files. Enclose the pattern in double quotes.
| -g, --log | Log file path | log_output.log
| -G, --log-level |Logging level (1: DEBUG, 2: INFO, 3: WARNING, 4: ERROR, 5: CRITICAL) | 2
| -o, --outputfolder | Output directory | output_longreadsum
| -t, --threads | The number of threads used | 1
| -Q, --outprefix | Output file prefix |  QC_
| --fontsize | Font size for plots | 14
| --markersize | Marker size for plots | 10

## WGS BAM
## RRMS BAM
## RNA-Seq BAM
## ONT POD5
## ONT FAST5
## Basecall summary
## FASTQ
## FASTA

Specifying input files:

```
usage: longreadsum [-h] {fa,fq,f5,f5s,seqtxt,bam,rrms} ...

Fast and comprehensive QC for long read sequencing data.

positional arguments:
  {fa,fq,f5,seqtxt,bam}
    fa                  FASTA file input
    fq                  FASTQ file input
    f5                  FAST5 file input
    f5s                 FAST5 file input with signal statistics output    
    seqtxt              sequencing_summary.txt input
    bam                 BAM file input
    rrms                RRMS BAM file input

optional arguments:
  -h, --help            show this help message and exit

Example with single inputs:
	longreadsum bam -i input.bam -o output_directory -t 12

Example with multiple inputs:
	longreadsum bam -I input1.bam, input2.bam -o output_directory
	longreadsum bam -P *.bam -o output_directory

RRMS example:
  longreadsum rrms --csv rrms_results.csv --input input.bam --output output_directory --threads 12

FAST5 signal mode example:
  longreadsum f5s --input input.fast5 --output output_directory
```


# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadSum/releases). 

# Getting help
Please refer to the [LongReadSum issue pages](https://github.com/WGLab/LongReadSum/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing LongReadSum
### Please cite the preprint below if you use our tool:

Perdomo, J. E., Ahsan, M. U., Liu, Q., Fang, L. & Wang, K. LongReadSum: A fast and flexible quality control and signal summarization tool for long-read sequencing data. bioRxiv, 2024.2008.2005.606643, [doi:10.1101/2024.08.05.606643](https://doi.org/doi:10.1101/2024.08.05.606643) (2024).
