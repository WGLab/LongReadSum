# LongReadSum: A fast and flexible QC tool for long read sequencing data

![image](https://user-images.githubusercontent.com/14855676/180858677-bba1dda7-15a2-4ba0-8ff5-6c954d00ba85.png)

[![build tests](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml/badge.svg)](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml)

LongReadSum supports FASTA, FASTQ, BAM, FAST5, and sequencing_summary.txt file formats for quick generation of QC data in HTML and text format.

## System requirements

## Software requirements
Please refer to `environment.yml` for detail. For your quick reference, LongReadSum needs
```
  - python=3.9
  - hdf5
  - htslib
  - swig
  - matplotlib
```

# Installation using Docker (recommended)
First, install [Docker](https://docs.docker.com/engine/install/).
Pull the latest image from Docker hub:

```
docker pull genomicslab/longreadsum
```

## Running

On Unix/Linux:
```
docker run -v local/inputs/:/inputs/ -it genomicslab/longreadsum bam -i /inputs/input.bam
```
Note that the `-v` command is required for Docker to find the input file. In the above BAM file example, the local directory `local/inputs/` containing the input file is mapped to a directory `/inputs/` in the Docker container. Thus, the input file is specified as `/inputs/input.bam`

# Installation using Anaconda
First install [Anaconda](https://www.anaconda.com/). Then follow the instructions below to install LongReadSum and its dependencies:

```
git clone https://github.com/WGLab/LongReadSum
cd LongReadSum
conda env create -f environment.yml

export PATH=~/miniconda3/envs/lrst_py39/bin:$PATH
conda activate lrst_py39
make

```
## Running
Activate the conda environment:

`conda activate lrst_py39`

To test that you are using the correct Python interpreter, run:

`which python`

This should point to the environment's Python interpreter path:

`~/miniconda3/envs/lrst_py39/bin/python`

If the path is incorrect, export its location to `PATH`:

`export PATH=~/miniconda3/envs/lrst_py39/bin:$PATH`

Then you can run LongReadSum using the following command:

`python /path/to/LongReadSum [arguments]`

# General Usage

Specifying input files:

```
usage: LongReadSum [-h] {fa,fq,f5,seqtxt,bam} ...

QC tools for long-read sequencing data

positional arguments:
  {fa,fq,f5,seqtxt,bam}
    fa                  FASTA file input
    fq                  FASTQ file input
    f5                  FAST5 file input
    f5s                 FAST5 file input with signal statistics output    
    seqtxt              sequencing_summary.txt input
    bam                 BAM file input

optional arguments:
  -h, --help            show this help message and exit

Example with single inputs:
	python LongReadSum bam -i path/to/input.bam -o /output_directory/

Example with multiple inputs:
	python LongReadSum bam -I "path/to/input1.bam, path/to/input2.bam" -o /output_directory/
	python LongReadSum bam -P "path/to/*.bam" -o /output_directory/
```

# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadSum/releases). 

# Getting help
Please refer to the [LongReadSum issue pages](https://github.com/WGLab/LongReadSum/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing LongReadSum
***Please cite the presentation below if you use our tool***

Perdomo, J. E., M. U. Ahsan, Q. Liu, L. Fang, K. Wang. *LongReadSum: A fast and flexible quality control tool for long-read sequencing data*. Poster presented at: American Society of Human Genetics (ASHG) Annual Meeting; 2022 October 25-29; Los Angeles Convention Center, Los Angeles, CA.
