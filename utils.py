import os, itertools
import matplotlib.pyplot as plt
from textwrap import wrap
import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots


def fmt(x):
    if abs(x)>=1e7:
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
    fig, axes = plt.subplots(len(data), sharey=True, figsize =(8, 6))
    
    numbers_list=[[x.n50_read_length, x.mean_read_length, x.median_read_length] for x in data]

    category=['N50', 'Mean', 'Median']    
    category_list=itertools.cycle([category])
    ylabel_list=itertools.cycle(['Length (bp)'])
    xlabel_list=itertools.cycle([None])
    subtitle_list=subtitles
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
        
    
def plot_base_counts(data, path, subtitles=None, categories=None):
    fig, axes = plt.subplots(len(data), figsize =(8, 6))

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
    #plt.ticklabel_format(axis='both',style='sci', scilimits=(0,0))
        
    for ax, numbers, category, xlabel, ylabel, subtitle in zip(fig.axes, numbers_list,category_list, xlabel_list, ylabel_list, subtitle_list):
        #ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set(ylabel=ylabel, xlabel=xlabel)
        ax.set_title(subtitle, pad=10)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelbottom=True)
        ax.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')

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
    mat=data.read_length_count
    mean, median, n50=int(data.mean_read_length), data.median_read_length, data.n50_read_length
    
    
    mat=np.array(mat)[:,np.newaxis]
    mat=mat[:np.max(np.nonzero(mat))+1]
    bin_size=1000
    mat_full=np.vstack([mat, np.zeros(bin_size-len(mat)%bin_size)[:,np.newaxis]])
    mat_full=mat_full.ravel()

    lengths=np.arange(0,len(mat_full))
    
    binsize=1000
    hist, bins = np.histogram(lengths, weights=mat_full, bins=np.arange(0,len(lengths)+1, bin_size))
    log_hist, log_bins = np.histogram(np.log10(lengths+1), weights=mat_full, bins=len(lengths)//binsize)
    
    fig = make_subplots(
    rows=2, cols=1,
    subplot_titles=("Read Length Histogram", "Log Read Length Histogram"), vertical_spacing = 0.18,)
    fig.update_layout(showlegend=False, autosize=False,
        width=800,
        height=600)
    
    fig.update_annotations(font_size=16)

    xd=bins[1:]
    customdata=np.dstack((bins[:-1],bins[1:], hist))[0,:,:]
    yd=hist
    y_max=np.max(hist)
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=customdata, hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata[2]:.0f}<extra></extra>', marker_color='#36a5c7'), row=1, col=1)

    fig.add_vline(mean, line_width=1, line_dash="dash",annotation_text='Mean', annotation_bgcolor="white", annotation_textangle=90, annotation_font=dict(size=10), row=1, col=1)
    fig.add_vline(median, line_width=1, line_dash="dash",annotation_text='Median', annotation_bgcolor="white", annotation_textangle=90, annotation_font=dict(size=10), row=1, col=1)
    fig.add_vline(n50, line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="white", annotation_textangle=90, annotation_font=dict(size=10), row=1, col=1)


    xd=log_bins[1:]
    customdata=np.dstack((np.power(10,log_bins)[:-1],np.power(10,log_bins)[1:],log_hist))[0,:,:]
    yd=log_hist
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=customdata, hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata[2]:.0f}<extra></extra>', marker_color='#36a5c7') , row=2, col=1)

    fig.add_vline(np.log10(mean), line_width=1, line_dash="dash",annotation_text='Mean', annotation_bgcolor="white", annotation_textangle=90, annotation_font=dict(size=10), row=2, col=1)
    fig.add_vline(np.log10(median), line_width=1, line_dash="dash",annotation_text='Median', annotation_bgcolor="white", annotation_textangle=90, annotation_font=dict(size=10), row=2, col=1)
    fig.add_vline(np.log10(n50), line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="white", annotation_textangle=90, annotation_font=dict(size=10), row=2, col=1)

    fig.update_xaxes(
        tickmode = 'array',
        tickvals = list(range(0, 12)),
        ticktext = ['0']+['{:,}'.format(10**x) for x in range(1, 12)],
        ticks="outside", row=2, col=1)

    fig.update_xaxes(ticks="outside", title_text='Read Length', title_standoff= 0)
    fig.update_yaxes(ticks="outside", title_text='Counts', title_standoff= 0)    

    fig.write_image(path,engine="kaleido")
    
    return fig.to_html(full_html=False)


def base_quality(data, path):
    
    xd=np.arange(256)
    yd=np.array(data.base_quality_distribution)
    fig = go.Figure()
    
    fig.add_trace(go.Bar(x=xd, y=yd, marker_color='#36a5c7'))
    
    fig.update_xaxes(ticks="outside", dtick=10, title_text='Base Quality', title_standoff= 0)
    fig.update_yaxes(ticks="outside", title_text='Number of bases', title_standoff= 0)
    
    fig.write_image(path,engine="kaleido")
    
    return fig.to_html(full_html=False)


def read_avg_base_quality(data, path):
    
    xd=np.arange(256)
    yd=np.array(data.read_average_base_quality_distribution)
    fig = go.Figure()
    
    fig.add_trace(go.Bar(x=xd, y=yd, marker_color='#36a5c7'))
    
    fig.update_xaxes(ticks="outside", dtick=10, title_text='Average Base Quality', title_standoff= 0)
    fig.update_yaxes(ticks="outside", title_text='Number of Reads', title_standoff= 0)
    
    fig.write_image(path,engine="kaleido")
    
    return fig.to_html(full_html=False)
