"""
lrst_global.py:
Format image filenames based on the provided output folder parameter.
Set the log error code variables.
"""


import os

prg_name = "LongReadSum"
originalError = '!!!Error: !!!!!! \n'

LOG_ALL = 0;
LOG_DEBUG = 1;
LOG_INFO = 2;
LOG_WARN = 3;
LOG_ERROR = 4;
LOG_FATAL = 5;
LOG_OFF = 6;

default_image_path = 'img/'
default_image_suf = '.png'

plot_filenames= {\
# for fq/fa
"read_length_distr": {'file':default_image_path+"read_length_distr"+default_image_suf, 'title':"Read Length", 'description':"The distribution of read length"},\
# for bam
"map_st": {'file':default_image_path+"map_st"+default_image_suf, 'title':"Map Information", 'description':"The statistics of mapped and unmapped reads"},\
                                  
"err_st": {'file':default_image_path+"err_st"+default_image_suf, 'title':"Base Alignment and Error Statistics", 'description':"Alignment statistics of mapped bases"},\

"read_length_st": {'file':default_image_path+"read_length_st"+default_image_suf, 'title':"Read Length Statistics", 'description':"Statistics of Read Lengths"},\

"base_st": {'file':default_image_path+"base_st"+default_image_suf, 'title':"Base Count Statistics", 'description':"Statistics of Base Counts", 'summary':""},\

"basic_info": {'file':default_image_path+"basic_info"+default_image_suf, 'title':"Basic Statistics", 'description':"Basic Statistics", 'summary':""},\

"read_length_hist":{'file':default_image_path+"read_length_hist"+default_image_suf, 'title':"Read Length Histogram", 'description':"Read Length Histogram", 'summary':""},

"base_quality":{'file':default_image_path+"base_quality"+default_image_suf, 'title':"Base Quality", 'description':"Base Quality"},
                 
"read_avg_base_quality":{'file':default_image_path+"read_avg_base_quality"+default_image_suf, 'title':"Average Read Base Quality", 'description':"Average Read Base Quality"},
                 
"pos_quality":{'file':default_image_path+"pos_quality"+default_image_suf, 'title':"Base Position Quality", 'description':"Base Position Quality"},
                                  
# for fast5
"nanopore_st": {'file':default_image_path+"nanopore_st"+default_image_suf, 'title':"ONT-Basecall", 'description':"The statistics of passed and failed reads"}\
}
