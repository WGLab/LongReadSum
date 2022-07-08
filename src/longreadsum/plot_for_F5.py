"""
plot_for_F5.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

import lrst_global
import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

from lib import lrst
from utils import *
    
def generate_bs( f5_output, para_dict ):
   lrst_global.plot_filenames["basic_st"] = {};
   lrst_global.plot_filenames["basic_st"]['file'] = ""
   lrst_global.plot_filenames["basic_st"]['title'] = "Table of Statistics Summary"
   lrst_global.plot_filenames["basic_st"]['description'] = "Basic Statistics for F5 summary"
   
   table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Passed</th><th>Failed</th><th>All</th></tr>\n</thead>"
   table_str += "\n<tbody>"
   int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td></tr>"
   double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
   table_str += int_str_for_format.format("#Total Reads", \
                 f5_output.f5_passed_long_read_info.long_read_info.total_num_reads, \
                 f5_output.f5_failed_long_read_info.long_read_info.total_num_reads, \
                 f5_output.f5_long_read_info.long_read_info.total_num_reads);
   table_str += int_str_for_format.format("#Total Bases", \
                 f5_output.f5_passed_long_read_info.long_read_info.total_num_bases, \
                 f5_output.f5_failed_long_read_info.long_read_info.total_num_bases, \
                 f5_output.f5_long_read_info.long_read_info.total_num_bases);
   table_str += int_str_for_format.format("Longest Read Length", \
                 f5_output.f5_passed_long_read_info.long_read_info.longest_read_length, \
                 f5_output.f5_failed_long_read_info.long_read_info.longest_read_length, \
                 f5_output.f5_long_read_info.long_read_info.longest_read_length);
   table_str += int_str_for_format.format("N50", \
                 f5_output.f5_passed_long_read_info.long_read_info.n50_read_length, \
                 f5_output.f5_failed_long_read_info.long_read_info.n50_read_length, \
                 f5_output.f5_long_read_info.long_read_info.n50_read_length);
   table_str += double_str_for_format.format("Mean Read Length", \
                 f5_output.f5_passed_long_read_info.long_read_info.mean_read_length, \
                 f5_output.f5_failed_long_read_info.long_read_info.mean_read_length, \
                 f5_output.f5_long_read_info.long_read_info.mean_read_length);
   table_str += int_str_for_format.format("Median Read Length", \
                 f5_output.f5_passed_long_read_info.long_read_info.median_read_length, \
                 f5_output.f5_failed_long_read_info.long_read_info.median_read_length, \
                 f5_output.f5_long_read_info.long_read_info.median_read_length);
   table_str += "\n</tbody>\n</table>"
   
   lrst_global.plot_filenames["basic_st"]['detail'] = table_str;


def f5_plot( f5_output, para_dict ):
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
 
    generate_bs( f5_output, para_dict );
 
    plot_read_length_stats([f5_output.f5_long_read_info.long_read_info, f5_output.f5_passed_long_read_info.long_read_info, f5_output.f5_failed_long_read_info.long_read_info], get_image_path('read_length_st'), subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    plot_base_counts([f5_output.f5_long_read_info.long_read_info, f5_output.f5_passed_long_read_info.long_read_info, f5_output.f5_failed_long_read_info.long_read_info], get_image_path('base_st'), subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    plot_basic_info([f5_output.f5_long_read_info.long_read_info, f5_output.f5_passed_long_read_info.long_read_info, f5_output.f5_failed_long_read_info.long_read_info], get_image_path('basic_info'), categories=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    histogram(f5_output.f5_long_read_info.long_read_info, get_image_path('read_length_hist'))
