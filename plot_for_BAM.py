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
    plt.title('Counts of different categories of reads and alignments')
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
    plt.title('Counts of different categories of base alignments')
    plt.tight_layout()

    plt.savefig(path)
    
    
def plot_read_length_stats(bam_output, path):
    fig, ax = plt.subplots(figsize =(8, 6))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    numbers=[5000,2000,1000]

    labels=['N50', 'Mean', 'Median']
    labels=[wrap(x) for x in labels]

    plt.barh(labels,numbers)
    
    for index, value in enumerate(numbers):
        plt.text(value+max(numbers)*0.01, index-0.08, f'{value:,}')
    
    
    #plt.ylabel('Base Allignment Type')
    plt.xlabel('Length')
    plt.title('Read Length Statistics')
    plt.tight_layout()

    plt.savefig(path)
    
def bam_plot( bam_output, para_dict ):
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
    
    print("num_primary_alignment: {}".format(bam_output.num_primary_alignment))
    
    plot_alignment_numbers(bam_output, get_image_path('map_st'))
    
    plot_errors(bam_output, get_image_path('err_st'))
    plot_read_length_stats(bam_output, get_image_path('read_length_st'))