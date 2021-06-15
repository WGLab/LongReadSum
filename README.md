# LongReadSum: a fast and flexible tools for data statistics of long-read sequencing data 

LongReadSum supports fq, fa, bam and fast5 inputs for quick generation of long-read data in a HTML format.

# System requirements
## Hardware requirements
There is specific hardware requirements to use LongReadSum if you can successfully install all dependent packages.

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
It is easy to install the dependent packages of LongReadSum using `annoconda`. Thus, please install `annoconda` first, and then follow the commands below to install LongReadSum.

```
git clone https://github.com/WGLab/LongReadSum
cd LongReadSum
conda env create -f environment.yml
source activate lrst_py39   #if you change conda env name, please replace `lrst_py39`
make
```

Then, you can run `python LongReadSum.py`


# General Usage
After installation, simply type `python LongReadSum/LongReadSum.py` will tell you the options.
```
usage: LongReadSum.py [-h] {fq,fa,bam,f5} ...

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
                                                python LongReadSum.py fq: with fq or fastq input
                                                python LongReadSum.py fa: with fa or fasta input
                                                python LongReadSum.py bam: with bam input
                                                python LongReadSum.py f5: with fast5 input
```


# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadSum/releases). 

# Getting help
Please refer to the [LongReadSum issue pages](https://github.com/WGLab/LongReadSum/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing LongReadSum
***Please cite the publication below if you use our tool***



