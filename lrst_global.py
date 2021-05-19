
import os

prg_name = "LongReadDS"
originalError = '!!!Error: !!!!!! \n'

LOG_ALL = 0;
LOG_DEBUG = 1;
LOG_INFO = 2;
LOG_WARN = 3;
LOG_ERROR = 4;
LOG_FATAL = 5;
LOG_OFF = 6;

default_image_path = './img/'
if not os.path.isdir( default_image_path ):
   os.makedirs( default_image_path )
default_image_suf = '.png'

plot_filenames= {\
# for fq/fa
"read_length_distr": [default_image_path+"read_length_distr"+default_image_suf, "Read Length", "The distribution of read length"],\
# for bam
"map_st": [default_image_path+"map_st"+default_image_suf, "Map Information", "The statistics of mapped and unmapped reads"],\

# for fast5
"nanopore_st":[default_image_path+"nanopore_st"+default_image_suf, "ONT-Basecall", "The statistics of passed and failed reads"]\
}


