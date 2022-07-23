# LongReadSum: A fast and flexible QC tool for long read sequencing data

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

First, make sure you are in the `lrst_py39` conda environment, and that you have exported its location to `PATH` as described above.

To test that you are using the correct Python interpreter, run:

`which python`

This should point to the environment's Python interpreter path:

`~/miniconda3/bin/python`

Then you can run LongReadSum using the following command:

`python /path/to/LongReadSum`

```
usage: LongReadSum [-h] {fq,fa,bam,f5} ...

Data analysis tools for long-read sequencing data

positional arguments:
  {fq,fa,bam,f5}
    fq            Show data analysis for fq files
    fa            Show data analysis for fa files
    bam           Show data analysis for bam files
    fast5         Show data analysis for FAST5 files
    seqtxt        Show data analysis for sequencing_summary.txt files

optional arguments:
  -h, --help      show this help message and exit

For example,
                                                LongReadSum fq: with fq or fastq input
                                                LongReadSum fa: with fa or fasta input
                                                LongReadSum bam: with bam input
                                                LongReadSum f5: with fast5 input
```


# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadSum/releases). 

# Getting help
Please refer to the [LongReadSum issue pages](https://github.com/WGLab/LongReadSum/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing LongReadSum
***Please cite the publication below if you use our tool***
