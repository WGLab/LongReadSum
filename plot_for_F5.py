import lrst_global
import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

import lrst;
from utils import *
    
def f5_plot( f5_output, para_dict ):
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
  
    plot_read_length_stats([f5_output.f5_long_read_info.long_read_info, f5_output.f5_passed_long_read_info.long_read_info, f5_output.f5_failed_long_read_info.long_read_info], get_image_path('read_length_st'), subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    plot_base_counts([f5_output.f5_long_read_info.long_read_info, f5_output.f5_passed_long_read_info.long_read_info, f5_output.f5_failed_long_read_info.long_read_info], get_image_path('base_st'), subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    
    plot_basic_info([f5_output.f5_long_read_info.long_read_info, f5_output.f5_passed_long_read_info.long_read_info, f5_output.f5_failed_long_read_info.long_read_info], get_image_path('basic_info'), categories=['All Reads', 'Passed Reads', 'Failed Reads'])
