import lrst_global
import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

import lrst;
from utils import *

def plot_alignment_numbers(data, path):
    fig, axes = plt.subplots(figsize =(8, 6))
    
    numbers_list=[[data.num_primary_alignment, data.num_supplementary_alignment, data.num_secondary_alignment,
                   data.num_reads_with_supplementary_alignment, data.num_reads_with_secondary_alignment,
                   data.num_reads_with_both_secondary_supplementary_alignment, data.forward_alignment, data.reverse_alignment]]
    
    category=['Primary Alignments', 'Supplementary Alignments', 'Secondary Alignments', 'Reads with Supplementary Alignments','Reads with Secondary Alignments','Reads with Secondary and Supplementary Alignments', 'Forward Alignments', 'Reverse Alignments']
    category=[wrap(x) for x in category]
    
    category_list=itertools.cycle([category])
    xlabel_list=itertools.cycle(['Counts'])
    ylabel_list=itertools.cycle(['Read/Alignment Type'])
    subtitle_list=[None]
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='h')
    
def plot_errors(bam_output, path):
    fig, axes = plt.subplots(1,1, figsize =(8, 6))

    numbers_list=[[bam_output.num_matched_bases, bam_output.num_mismatched_bases, bam_output.num_ins_bases, bam_output.num_del_bases, bam_output.num_clip_bases]]
    
    category=['Matched Bases', 'Mismatched Bases', 'Inserted Bases', 'Deleted Bases', 'Clipped Bases']
    category=[wrap(x) for x in category]
    
    category_list=itertools.cycle([category])
    xlabel_list=itertools.cycle(['Counts'])
    ylabel_list=itertools.cycle([None])
    subtitle_list=[None]
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='h')

    
def bam_plot( bam_output, para_dict ):
    
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
    
    print("num_primary_alignment: {}".format(bam_output.num_primary_alignment))
    
    plot_alignment_numbers(bam_output, get_image_path('map_st'))
    
    plot_errors(bam_output, get_image_path('err_st'))
    plot_read_length_stats(bam_output, get_image_path('read_length_st'))
    plot_base_counts(bam_output, get_image_path('base_st'))
    plot_basic_info(bam_output.mapped_long_read_info,bam_output.unmapped_long_read_info,['Mapped Reads', 'Unmapped Reads'], get_image_path('basic_info'))
    
    