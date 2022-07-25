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

# Installation
First install [Anaconda](https://www.anaconda.com/). Then follow the instructions below to install LongReadSum and its dependencies:

```
git clone https://github.com/WGLab/LongReadSum
conda env create -f LongReadSum/environment.yml

export PATH=~/miniconda3/envs/lrst_py39/bin:$PATH
conda activate lrst_py39
make

```


# General Usage

First, make sure you have activated the `lrst_py39` conda environment.

To test that you are using the correct Python interpreter, run:

`which python`

This should point to the environment's Python interpreter path:

`~/miniconda3/envs/lrst_py39/bin/python`

If the path is incorrect, export its location to `PATH` as described above.

Then you can run LongReadSum using the following command:

`python /path/to/LongReadSum [arguments]`

Specifying input files:

```
usage: LongReadSum [-h] {fa,fq,f5,seqtxt,bam} ...

QC tools for long-read sequencing data

positional arguments:
  {fa,fq,f5,seqtxt,bam}
    fa                  FASTA file input
    fq                  FASTQ file input
    f5                  FAST5 file input
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
***Please cite the publication below if you use our tool***
