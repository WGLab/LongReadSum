
prg_name = "LongReadDS"
originalError = '!!!Error: !!!!!! \n'

LOG_ALL = 0;
LOG_DEBUG = 1;
LOG_INFO = 2;
LOG_WARN = 3;
LOG_ERROR = 4;
LOG_FATAL = 5;
LOG_OFF = 6;


plot_file_name= {\
# for fq/fa
"read_length_distr": ["read_length_distr", "The distribution of read length"],\
# for bam
"map_st": ["map_st", "The statistics of mapped and unmapped reads"],\

# for fast5
"nanopore_st":["nanopore_st", "The statistics of passed and failed reads"]\
}


