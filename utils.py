import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np

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

def plot_read_length_stats(data, path, subtitles=None, categories=None):
    fig, axes = plt.subplots(len(data), sharex=True, figsize =(8, 6))
    
    numbers_list=[[x.n50_read_length, x.mean_read_length, x.median_read_length] for x in data]

    category=['N50', 'Mean', 'Median']    
    category_list=itertools.cycle([category])
    ylabel_list=itertools.cycle(['Length (bp)'])
    xlabel_list=itertools.cycle([None])
    subtitle_list=subtitles
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
        
    
def plot_base_counts(data, path, subtitles=None, categories=None):
    fig, axes = plt.subplots(len(data), sharey=True, figsize =(8, 6))

    numbers_list=[[x.total_a_cnt, x.total_c_cnt, x.total_g_cnt, x.total_tu_cnt, x.total_n_cnt] for x in data]
    
    category_list=itertools.cycle([['A', 'C', 'G', 'T/U', 'N']])
    xlabel_list=itertools.cycle([None])
    ylabel_list=itertools.cycle(['Counts'])
    subtitle_list=subtitles
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
    #lrst_global.plot_filenames['base_st']['summary']='GC Content: {:.2%}'.format(bam_output.mapped_long_read_info.gc_cnt)

def plot_basic_info(data, path, subtitles=None, categories=None):
    fig, axes= plt.subplots(2,2, figsize =(8, 6))
    

    numbers_list=[[x.total_num_reads, x.total_num_bases, x.longest_read_length, x.gc_cnt] for x in data]
    
    numbers_list=zip(*numbers_list)
    
    category_list=itertools.cycle([categories])
    subtitle_list=['Number of Reads', 'Number of Bases', 'Longest Read', 'GC Content']
    xlabel_list=['Count', 'Count', 'Length (bp)', '%']
    ylabel_list=itertools.cycle([None])
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='h')
    
    
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
    
def histogram(data, path):
    plt.subplots_adjust(hspace=0.5,wspace=0.5)
    fig, axes= plt.subplots(2,1, figsize=(8,6))
    colors = {'Mean':'g', 'Median':'r', 'N50':'k'}
    
    mat=data.read_length_count
    stats = {'Mean':int(data.mean_read_length), 'Median':data.median_read_length, 'N50':data.n50_read_length}
    
    
    mat=np.array(mat)[:,np.newaxis]
    mat=mat[:np.max(np.nonzero(mat))+1]
    bin_size=500
    mat_full=np.vstack([mat, np.zeros(bin_size-len(mat)%bin_size)[:,np.newaxis]])
    mat_full=mat_full.ravel()

    lengths=np.arange(0,len(mat_full))
    
    ax=axes[0]
    ax.set_title('Read Length Histogram')
    ax.set(ylabel='Frequency', xlabel='Length(bp)')
    ax.hist(lengths, weights=mat_full/np.sum(mat_full), bins=np.arange(0,len(lengths), bin_size))
    for x in stats:
        ax.axvline(x=stats[x], color=colors[x], linestyle='--',linewidth=1, label='{} = {}'.format(x,stats[x]))

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.legend()
    
    ax=axes[1]
    
    logbins = np.logspace(0,np.log10(lengths[-1]),len(lengths)//500)
    
    ax.set_title('Log Read Length Histogram')
    ax.set(ylabel='Frequency', xlabel='Length(bp)')

    ax.hist(lengths,weights=mat_full/np.sum(mat_full), bins=logbins)

    ax.set_xscale('log')

    for x in stats:
        ax.axvline(x=stats[x], color=colors[x], linestyle='--',linewidth=1)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(xmin=100)
    
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=2.0)
    plt.savefig(path) 