import lrst_global
import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

import lrst;


def fmt(x):
    if abs(x)>=1e9:
        return f'{x:.3E}'
    
    elif abs(x)>=1e3 or float(x).is_integer():
        return f'{x:,.0f}'
    
    elif abs(x)<0.01:
        return f'{x:.3E}'
    
    else:
        return f'{x:.4}'

def wrap(s):
    l=s.split(' ')
    split=list(zip(*[iter(l)]*3))
    if len(l)%3:
        split.append(tuple(l[-(len(l)%3):]))
    return '\n'.join([' '.join(x) for x in split])

    
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
    
def plot_read_length_stats(data, path):
    fig, axes = plt.subplots(2, sharex=True, figsize =(8, 6))
    
    numbers_list=[[data.mapped_long_read_info.n50_read_length,
                 data.mapped_long_read_info.mean_read_length,
                 data.mapped_long_read_info.median_read_length],
    
                 [data.unmapped_long_read_info.n50_read_length,
                 data.unmapped_long_read_info.mean_read_length,
                 data.unmapped_long_read_info.median_read_length]]
    
    category=['N50', 'Mean', 'Median']
    category=[wrap(x) for x in category]
    
    category_list=itertools.cycle([category])
    ylabel_list=itertools.cycle(['Length (bp)'])
    xlabel_list=itertools.cycle([None])
    subtitle_list=["Mapped Reads","Unmapped Reads"]
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
        
    
def plot_base_counts(data, path):
    fig, axes = plt.subplots(2, sharey=True, figsize =(8, 6))

    numbers_list=[[data.mapped_long_read_info.total_a_cnt, 
             data.mapped_long_read_info.total_c_cnt,
             data.mapped_long_read_info.total_g_cnt,
             data.mapped_long_read_info.total_tu_cnt, 
             data.mapped_long_read_info.total_n_cnt],
                  
             [data.unmapped_long_read_info.total_a_cnt, 
             data.unmapped_long_read_info.total_c_cnt,
             data.unmapped_long_read_info.total_g_cnt,
             data.unmapped_long_read_info.total_tu_cnt, 
             data.unmapped_long_read_info.total_n_cnt]]

    category_list=itertools.cycle([['A', 'C', 'G', 'T/U', 'N']])
    ylabel_list=itertools.cycle(['Base'])
    xlabel_list=itertools.cycle(['Counts'])
    subtitle_list=["Mapped Reads","Unmapped Reads"]
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
    #lrst_global.plot_filenames['base_st']['summary']='GC Content: {:.2%}'.format(bam_output.mapped_long_read_info.gc_cnt)

def plot_basic_info(data, path):
    fig, axes= plt.subplots(2,2, figsize =(8, 6))
    

    numbers_list=[( data.mapped_long_read_info.total_num_reads,  data.unmapped_long_read_info.total_num_reads),
             (data.mapped_long_read_info.total_num_bases,  data.unmapped_long_read_info.total_num_bases),
             (data.mapped_long_read_info.longest_read_length,  data.unmapped_long_read_info.longest_read_length),
             (data.mapped_long_read_info.gc_cnt,  data.unmapped_long_read_info.gc_cnt)]
    
    category_list=itertools.cycle([['Mapped Reads', 'Unmapped Reads']])
    subtitle_list=['Number of Reads', 'Number of Bases', 'Longest Read', 'GC Content']
    ylabel_list=['Count', 'Count', 'Length (bp)', '%']
    xlabel_list=itertools.cycle([None])
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
    

    
def read_distr(data,path):
    print()
    
    
def bar_plot(fig, numbers_list,category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='v', print_value =True):
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    for ax, numbers, category, xlabel, ylabel, subtitle in zip(fig.axes, numbers_list,category_list, xlabel_list, ylabel_list, subtitle_list):
        ax.set(ylabel=ylabel, xlabel=xlabel)
        ax.set_title(subtitle, pad=10)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelbottom=True)
        
        if orientation=='h':
            ax.barh(category, numbers)
            for index, value in enumerate(numbers):
                ax.text(value, index, ' %s' %fmt(value) if print_value else '')
            plt.tight_layout()
            
        elif orientation=='v':
            ax.bar(category, numbers)
            for index, value in enumerate(numbers):
                ax.text(index, value+max(numbers)*0.02, '%s' %fmt(value) if print_value else '', ha='center')

    
    plt.savefig(path) 
    
def bam_plot( bam_output, para_dict ):
    
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
    
    print("num_primary_alignment: {}".format(bam_output.num_primary_alignment))
    
    plot_alignment_numbers(bam_output, get_image_path('map_st'))
    
    plot_errors(bam_output, get_image_path('err_st'))
    plot_read_length_stats(bam_output, get_image_path('read_length_st'))
    plot_base_counts(bam_output, get_image_path('base_st'))
    plot_basic_info(bam_output, get_image_path('basic_info'))
    
    x=0
    p=0
    l=np.array(bam_output.mapped_long_read_info.read_length_count)
    with open('ff','w') as file:
        for x in l:
            file.write('%d\n' %x)