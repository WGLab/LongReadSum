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

def plot_read_length_stats(f5_output, path):
    fig, (ax_passed, ax_failed) = plt.subplots(2, sharex=True, figsize =(8, 6))
    
    numbers_passed=[f5_output.f5_passed_long_read_info.long_read_info.n50_read_length,
                 f5_output.f5_passed_long_read_info.long_read_info.mean_read_length,
                 f5_output.f5_passed_long_read_info.long_read_info.median_read_length]
    
    numbers_failed=[f5_output.f5_failed_long_read_info.long_read_info.n50_read_length,
                 f5_output.f5_failed_long_read_info.long_read_info.mean_read_length,
                 f5_output.f5_failed_long_read_info.long_read_info.median_read_length]
    
    labels=['N50', 'Mean', 'Median']
    labels=[wrap(x) for x in labels]

    for ax, numbers in [(ax_passed, numbers_passed), (ax_failed, numbers_failed)]:
        ax.set(xlabel='Length (bp)')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.barh(labels, numbers)
        ax.tick_params(labelbottom=True)

        for index, value in enumerate(numbers):
            ax.text(value+max(numbers)*0.01, index-0.08, f'{value:,}')

    #plt.ylabel('Base Allignment Type')

    plt.title('Read Length Statistics')

    ax_passed.set_title("Passed Reads")
    ax_failed.set_title("Failed Reads")

    plt.tight_layout()
    
    plt.savefig(path)
    
def plot_basic_info(f5_output, path):
    fig, axes= plt.subplots(2,2, figsize =(8, 6))
    
    numbers_list=[( f5_output.f5_passed_long_read_info.long_read_info.total_num_reads,  f5_output.f5_failed_long_read_info.long_read_info.total_num_reads),
             (f5_output.f5_passed_long_read_info.long_read_info.total_num_bases,  f5_output.f5_failed_long_read_info.long_read_info.total_num_bases),
             (f5_output.f5_passed_long_read_info.long_read_info.longest_read_length,  f5_output.f5_failed_long_read_info.long_read_info.longest_read_length),
             (f5_output.f5_passed_long_read_info.long_read_info.gc_cnt,  f5_output.f5_failed_long_read_info.long_read_info.gc_cnt)]
    
    xlabels=['Passed Reads', 'Failed Reads']
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
    
def f5_plot( f5_output, para_dict ):
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])
  
    print("test array {}".format( f5_output.f5_long_read_info.read_quality_distribution[4]))
 
    plot_read_length_stats(f5_output, get_image_path('read_length_st'))
    plot_basic_info(f5_output, get_image_path('basic_info'))
