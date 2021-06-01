import lrst_global
import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

import lrst;
from utils import *


def generate_bs( fq_output, para_dict ):
   lrst_global.plot_filenames["basic_st"] = {};
   lrst_global.plot_filenames["basic_st"]['file'] = ""
   lrst_global.plot_filenames["basic_st"]['title'] = "Basic Statistics"
   lrst_global.plot_filenames["basic_st"]['description'] = "Basic Statistics for FQ summary"

   table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
   table_str += "\n<tbody>"
   int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
   double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
   table_str += int_str_for_format.format("#Total Reads", \
                 fq_output.long_read_info.total_num_reads);
   table_str += int_str_for_format.format("#Total Bases", \
                 fq_output.long_read_info.total_num_bases);
   table_str += int_str_for_format.format("Longest Read Length", \
                 fq_output.long_read_info.longest_read_length);
   table_str += int_str_for_format.format("N50", \
                 fq_output.long_read_info.n50_read_length);
   table_str += double_str_for_format.format("GC Content(%)", \
                 fq_output.long_read_info.gc_cnt*100);
   table_str += double_str_for_format.format("Mean Read Length", \
                 fq_output.long_read_info.mean_read_length);
   table_str += int_str_for_format.format("Median Read Length", \
                 fq_output.long_read_info.median_read_length);
   table_str += "\n</tbody>\n</table>"

   lrst_global.plot_filenames["basic_st"]['detail'] = table_str;
    
def fq_plot( fq_output, para_dict ):
    
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
    
    #print("num_primary_alignment: {}".format(fq_output.num_primary_alignment))
    generate_bs( fq_output, para_dict)

    plot_read_length_stats([fq_output.long_read_info], get_image_path('read_length_st'), subtitles=['Long Reads'])
    plot_base_counts([fq_output.long_read_info], get_image_path('base_st'), subtitles=['Long Reads'])
    plot_basic_info([fq_output.long_read_info], get_image_path('basic_info'), categories=['Long Reads'])
    
    