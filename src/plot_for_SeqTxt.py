"""
plot_for_SeqTxt.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

from src import lrst_global
from src.utils import *
from lib import lrst

def generate_bs( output_statistics, para_dict ):
   lrst_global.plot_filenames["basic_st"] = {};
   lrst_global.plot_filenames["basic_st"]['file'] = ""
   lrst_global.plot_filenames["basic_st"]['title'] = "Table of Statistics Summary"
   lrst_global.plot_filenames["basic_st"]['description'] = "Basic Statistics for F5 summary"
   
   table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Passed</th><th>Failed</th><th>All</th></tr>\n</thead>"
   table_str += "\n<tbody>"
   int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td></tr>"
   double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
   table_str += int_str_for_format.format("#Total Reads", \
                 output_statistics.passed_long_read_info.long_read_info.total_num_reads, \
                 output_statistics.failed_long_read_info.long_read_info.total_num_reads, \
                 output_statistics.all_long_read_info.long_read_info.total_num_reads);
   table_str += int_str_for_format.format("#Total Bases", \
                 output_statistics.passed_long_read_info.long_read_info.total_num_bases, \
                 output_statistics.failed_long_read_info.long_read_info.total_num_bases, \
                 output_statistics.all_long_read_info.long_read_info.total_num_bases);
   table_str += int_str_for_format.format("Longest Read Length", \
                 output_statistics.passed_long_read_info.long_read_info.longest_read_length, \
                 output_statistics.failed_long_read_info.long_read_info.longest_read_length, \
                 output_statistics.all_long_read_info.long_read_info.longest_read_length);
   table_str += int_str_for_format.format("N50", \
                 output_statistics.passed_long_read_info.long_read_info.n50_read_length, \
                 output_statistics.failed_long_read_info.long_read_info.n50_read_length, \
                 output_statistics.all_long_read_info.long_read_info.n50_read_length);
   table_str += double_str_for_format.format("Mean Read Length", \
                 output_statistics.passed_long_read_info.long_read_info.mean_read_length, \
                 output_statistics.failed_long_read_info.long_read_info.mean_read_length, \
                 output_statistics.all_long_read_info.long_read_info.mean_read_length);
   table_str += int_str_for_format.format("Median Read Length", \
                 output_statistics.passed_long_read_info.long_read_info.median_read_length, \
                 output_statistics.failed_long_read_info.long_read_info.median_read_length, \
                 output_statistics.all_long_read_info.long_read_info.median_read_length);
   table_str += "\n</tbody>\n</table>"
   
   lrst_global.plot_filenames["basic_st"]['detail'] = table_str;


def plot( output_statistics, para_dict ):
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
 
    generate_bs( output_statistics, para_dict );
 
    plot_read_length_stats([output_statistics.all_long_read_info.long_read_info, output_statistics.passed_long_read_info.long_read_info, output_statistics.failed_long_read_info.long_read_info], get_image_path('read_length_st'), subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    plot_base_counts([output_statistics.all_long_read_info.long_read_info, output_statistics.passed_long_read_info.long_read_info, output_statistics.failed_long_read_info.long_read_info], get_image_path('base_st'), subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    plot_basic_info([output_statistics.all_long_read_info.long_read_info, output_statistics.passed_long_read_info.long_read_info, output_statistics.failed_long_read_info.long_read_info], get_image_path('basic_info'), categories=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    histogram(output_statistics.all_long_read_info.long_read_info, get_image_path('read_length_hist'))
