# LongReadSum: A fast and flexible QC tool for long read sequencing data

![image](https://user-images.githubusercontent.com/14855676/180858677-bba1dda7-15a2-4ba0-8ff5-6c954d00ba85.png)

[![build tests](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml/badge.svg)](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml)

LongReadSum supports FASTA, FASTQ, BAM, FAST5, and sequencing_summary.txt file formats for quick generation of QC data in HTML and text format.

## Software requirements
Please refer to `environment.yml` for detail. For your quick reference, LongReadSum needs
```
  - python=3.9
  - hdf5
  - htslib
  - swig
  - matplotlib
```
# Installation using Anaconda (Linux 64-bit)
First, install [Anaconda](https://www.anaconda.com/).
LongReadSum can be installed using the following command:

```
conda install -c bioconda -c wglab longreadsum=1.3.0
```

# Installation using Docker
First, install [Docker](https://docs.docker.com/engine/install/).
Pull the latest image from Docker hub:

```
docker pull genomicslab/longreadsum
```

## Running

On Unix/Linux:
```
docker run -v C:/Users/.../DataDirectory:/mnt/ -it genomicslab/longreadsum bam -i /mnt/input.bam -o /mnt/output
```
Note that the `-v` command is required for Docker to find the input file. Use a directory under `C:/Users/` to ensure volume files are mounted correctly. In the above example, the local directory `C:/Users/.../DataDirectory` containing the input file `input.bam` is mapped to a directory `/mnt/` in the Docker container. Thus, the input file and output directory arguments are relative to the `/mnt/` directory, but the output files will also be saved locally in `C:/Users/.../DataDirectory` under the specified subdirectory `output`.


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

If you are using FAST5 files with VBZ compression, you will need to download and install the VBZ plugin corresponding to your architecture:
https://github.com/nanoporetech/vbz_compression/releases

For example:

```
wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
tar -xf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
```

Finally, add the plugin to your path:
```
export HDF5_PLUGIN_PATH=/full/path/to/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
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
