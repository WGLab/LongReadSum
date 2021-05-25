import lrst_global
import os
import matplotlib.pyplot as plt
from textwrap import wrap

import lrst;

def wrap(s):
    l=s.split(' ')
    split=list(zip(*[iter(l)]*3))
    if len(l)%3:
        split.append(tuple(l[-(len(l)%3):]))
    return '\n'.join([' '.join(x) for x in split])

def plot_alignment_numbers(bam_output, path):
    fig, ax = plt.subplots(figsize =(8, 6))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    numbers=[bam_output.num_primary_alignment, bam_output.num_supplementary_alignment,bam_output.num_secondary_alignment , bam_output.num_reads_with_supplementary_alignment , bam_output.num_reads_with_secondary_alignment ,bam_output.num_reads_with_both_secondary_supplementary_alignment ,bam_output.forward_alignment ,bam_output.reverse_alignment]
    
    labels=['Primary Alignments', 'Supplementary Alignments', 'Secondary Alignments', 'Reads with Supplementary Alignments','Reads with Secondary Alignments','Reads with Secondary and Supplementary Alignments', 'Forward Alignments', 'Reverse Alignments']
    labels=[wrap(x) for x in labels]
    
    plt.barh(labels,numbers)
    
    for index, value in enumerate(numbers):
        plt.text(value+max(numbers)*0.01, index-0.08, f'{value:,}')
    
   
    plt.ylabel('Read/Alignment Type')
    plt.xlabel('Counts')
    plt.tight_layout()

    plt.savefig(path)

def plot_errors(bam_output, path):
    fig, ax = plt.subplots(figsize =(8, 6))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    numbers=[bam_output.num_matched_bases, bam_output.num_mismatched_bases, bam_output.num_ins_bases, bam_output.num_del_bases, bam_output.num_clip_bases]

    labels=['Matched Bases', 'Mismatched Bases', 'Inserted Bases', 'Deleted Bases', 'Clipped Bases']
    labels=[wrap(x) for x in labels]

    plt.barh(labels,numbers)
    
    for index, value in enumerate(numbers):
        plt.text(value+max(numbers)*0.01, index-0.08, f'{value:,}')
    
    
    plt.ylabel('Base Allignment Type')
    plt.xlabel('Counts')
    plt.tight_layout()

    plt.savefig(path)
    
    
def plot_read_length_stats(bam_output, path):
    fig, (ax_map, ax_unmap) = plt.subplots(2, sharex=True, figsize =(8, 6))
    
    numbers_mapped=[bam_output.mapped_long_read_info.n50_read_length,
                 bam_output.mapped_long_read_info.mean_read_length,
                 bam_output.mapped_long_read_info.median_read_length]
    
    numbers_unmapped=[bam_output.unmapped_long_read_info.n50_read_length,
                 bam_output.unmapped_long_read_info.mean_read_length,
                 bam_output.unmapped_long_read_info.median_read_length]
    
    labels=['N50', 'Mean', 'Median']
    labels=[wrap(x) for x in labels]

    for ax, numbers in [(ax_map, numbers_mapped), (ax_unmap, numbers_unmapped)]:
        ax.set(xlabel='Length (bp)')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.barh(labels, numbers)
        ax.tick_params(labelbottom=True)

        for index, value in enumerate(numbers):
            ax.text(value+max(numbers)*0.01, index-0.08, f'{value:,}')

    #plt.ylabel('Base Allignment Type')

    plt.title('Read Length Statistics')

    ax_map.set_title("Mapped Reads")
    ax_unmap.set_title("Unmapped Reads")

    plt.tight_layout()
    
    plt.savefig(path)
    
def plot_base_counts(bam_output, path):
    fig, (ax_map, ax_unmap) = plt.subplots(2, sharex=True, figsize =(8, 6))

    numbers_mapped=[bam_output.mapped_long_read_info.total_a_cnt, 
             bam_output.mapped_long_read_info.total_c_cnt,
             bam_output.mapped_long_read_info.total_g_cnt,
             bam_output.mapped_long_read_info.total_tu_cnt, 
             bam_output.mapped_long_read_info.total_n_cnt]
    
    numbers_unmapped=[bam_output.unmapped_long_read_info.total_a_cnt, 
             bam_output.unmapped_long_read_info.total_c_cnt,
             bam_output.unmapped_long_read_info.total_g_cnt,
             bam_output.unmapped_long_read_info.total_tu_cnt, 
             bam_output.unmapped_long_read_info.total_n_cnt]

    labels=['A', 'C', 'G', 'T/U', 'N']
 
    for ax, numbers in [(ax_map, numbers_mapped), (ax_unmap, numbers_unmapped)]:
        ax.set(xlabel='Base', ylabel='Counts')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.barh(labels, numbers)
        ax.tick_params(labelbottom=True)

        for index, value in enumerate(numbers):
            ax.text(value+max(numbers)*0.01, index-0.08, f'{value:,}')
            
    plt.title('Base Counts')
    ax_map.set_title("Mapped Reads")
    ax_unmap.set_title("Unmapped Reads")
    plt.tight_layout()
    
    plt.savefig(path)    
    
    #lrst_global.plot_filenames['base_st']['summary']='GC Content: {:.2%}'.format(bam_output.mapped_long_read_info.gc_cnt)

def plot_basic_info(bam_output, path):
    fig, axes= plt.subplots(2,2, figsize =(8, 6))
    
    numbers_list=[( bam_output.mapped_long_read_info.total_num_reads,  bam_output.unmapped_long_read_info.total_num_reads),
             (bam_output.mapped_long_read_info.total_num_bases,  bam_output.unmapped_long_read_info.total_num_bases),
             (bam_output.mapped_long_read_info.longest_read_length,  bam_output.unmapped_long_read_info.longest_read_length),
             (bam_output.mapped_long_read_info.gc_cnt,  bam_output.unmapped_long_read_info.gc_cnt)]
    
    xlabels=['Mapped Reads', 'Unmapped Reads']
    subtitles=['Number of Reads', 'Number of Bases', 'Longest Read', 'GC Content']
    ylabels=['Count', 'Count', 'Length (bp)', '%']

    for ax, numbers, ylabel, subtitle in zip(axes.flat, numbers_list, ylabels, subtitles):
        ax.set(ylabel=ylabel)
        ax.set_title(subtitle)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.bar(xlabels, numbers)
        ax.tick_params(labelbottom=True)

        for index, value in enumerate(numbers):
            ax.text(index-0.08, value+max(numbers)*0.01, f'{value:,}')
    plt.tight_layout()
    
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
