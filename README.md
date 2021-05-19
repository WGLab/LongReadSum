# LongReadDS: a fast and flexible tools for data statistics of long-read sequencing data 

LongReadDS supports fq, fa, bam and fast5 inputs for quick generation of long-read data in a HTML format.

# System requirements
## Hardware requirements
There is specific hardware requirements to use LongReadDS if you can successfully install all dependent packages.

## Software requirements
Please refer to `environment.yml` for detail. For your quick reference, LongReadDS needs
```
  - python=3.9
  - hdf5
  - htslib
  - swig
  - matplotlib
```

# Installation
It is easy to install the dependent packages of LongReadDS using `annoconda`. Thus, please install `annoconda` first, and then follow the commands below to install LongReadDS.

```
git clone https://github.com/WGLab/LongReadDS
cd LongReadDS
conda env create -f environment.yml
source activate lrst_py39   #if you change conda env name, please replace `lrst_py39`
make
```

Then, you can run `python LongReadDS.py`


# General Usage
After installation, simply type `python LongReadDS/bin/LongReadDS.py` will tell you the options.
```
usage: LongReadDS.py [-h] {fq,fa,bam,f5} ...

Data analysis tools for long-read sequencing data

positional arguments:
  {fq,fa,bam,f5}
    fq            Show data analysis for fq files
    fa            Show data analysis for fa files
    bam           Show data analysis for bam files
    f5            Show data analysis for f5 files

optional arguments:
  -h, --help      show this help message and exit

For example,
                                                python LongReadDS.py fq: with fq or fastq input
                                                python LongReadDS.py fa: with fa or fasta input
                                                python LongReadDS.py bam: with bam input
                                                python LongReadDS.py f5: with fast5 input
```


# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadDS/releases). 

# Getting help
Please refer to the [LongReadDS issue pages](https://github.com/WGLab/LongReadDS/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing LongReadDS
***Please cite the publication below if you use our tool***



