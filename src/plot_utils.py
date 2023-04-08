import os
import logging
import numpy as np
import itertools

import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots

if __package__ == 'src':
    from src import lrst_global
else:
    import lrst_global

logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)


def setDefaultFontSize(font_size):
    """Set the default font size for matplotlib plots."""
    plt.rcParams.update({'font.size': font_size})


def fmt(x):
    """Format numbers for plots."""
    format_x = np.format_float_scientific(x, exp_digits=4)
    return format_x


def wrap(s):
    l = s.split(' ')
    split = list(zip(*[iter(l)] * 3))
    if len(l) % 3:
        split.append(tuple(l[-(len(l) % 3):]))
    return '\n'.join([' '.join(x) for x in split])


def plot_read_length_stats(data, path, subtitles=None, categories=None):
    fig, axes = plt.subplots(len(data), sharey=True, figsize=(8, 6))

    numbers_list = [[x.n50_read_length, x.mean_read_length, x.median_read_length] for x in data]

    category = ['N50', 'Mean', 'Median']
    category_list = itertools.cycle([category])
    ylabel_list = itertools.cycle(['Length (bp)'])
    xlabel_list = itertools.cycle([None])
    subtitle_list = subtitles
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)


def plot_base_counts(data, path, subtitles=None, categories=None):
    fig, axes = plt.subplots(len(data), figsize=(8, 6))

    numbers_list = [[x.total_a_cnt, x.total_c_cnt, x.total_g_cnt, x.total_tu_cnt, x.total_n_cnt] for x in data]

    category_list = itertools.cycle([['A', 'C', 'G', 'T/U', 'N']])
    xlabel_list = itertools.cycle([None])
    ylabel_list = itertools.cycle(['Counts'])
    subtitle_list = subtitles
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path)
    # lrst_global.plot_filenames['base_st']['summary']='GC Content: {:.2%}'.format(bam_output.mapped_long_read_info.gc_cnt)


def plot_basic_info(data, path, subtitles=None, categories=None):
    fig, axes = plt.subplots(2, 2, figsize=(8, 6))
    numbers_list = [[x.total_num_reads, x.total_num_bases, x.longest_read_length, x.gc_cnt] for x in data]
    numbers_list = zip(*numbers_list)

    category_list = itertools.cycle([categories])
    subtitle_list = ['Number of Reads', 'Number of Bases', 'Longest Read', 'GC Content']
    xlabel_list = ['Count', 'Count', 'Length (bp)', '%']
    ylabel_list = itertools.cycle([None])
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='h')


def bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='v',
             print_value=True):
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    # plt.ticklabel_format(axis='both',style='sci', scilimits=(0,0))

    for ax, numbers, category, xlabel, ylabel, subtitle in zip(fig.axes, numbers_list, category_list, xlabel_list,
                                                               ylabel_list, subtitle_list):
        # ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set(ylabel=ylabel, xlabel=xlabel)
        ax.set_title(subtitle, pad=10)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.tick_params(labelbottom=True)
        ax.ticklabel_format(style='sci', scilimits=(-3, 4), axis='both')

        if orientation == 'h':
            ax.barh(category, numbers)
            for index, value in enumerate(numbers):
                ax.text(value, index, ' %s' % fmt(value) if print_value else '')
            plt.tight_layout()

        elif orientation == 'v':
            ax.bar(category, numbers)
            for index, value in enumerate(numbers):
                ax.text(index, value + max(numbers) * 0.02, '%s' % fmt(value) if print_value else '', ha='center')

    plt.savefig(path)


def histogram(data, path, font_size):
    """Plot a histogram."""
    annotation_size = 10  # Annotation font size
    mat = data.read_length_count
    mean, median, n50 = int(data.mean_read_length), data.median_read_length, data.n50_read_length

    mat = np.array(mat)[:, np.newaxis]
    mat = mat[:data.longest_read_length + 1]
    bin_size = 1000
    mat_full = np.vstack([mat, np.zeros(bin_size - len(mat) % bin_size)[:, np.newaxis]])
    mat_full = mat_full.ravel()

    lengths = np.arange(0, len(mat_full))

    binsize = 1000
    hist, bins = np.histogram(lengths, weights=mat_full, bins=np.arange(0, len(lengths) + 1, bin_size))
    log_hist, log_bins = np.histogram(np.log10(lengths + 1), weights=mat_full, bins=len(lengths) // binsize)

    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Read Length Histogram", "Log Read Length Histogram"), vertical_spacing=0.3)
    fig.update_layout(showlegend=False, autosize=True)

    xd = bins[1:]
    customdata = np.dstack((bins[:-1], bins[1:], hist))[0, :, :]
    yd = hist
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=customdata,
                         hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata[2]:.0f}<extra></extra>',
                         marker_color='#36a5c7'), row=1, col=1)

    fig.add_vline(mean, line_width=1, line_dash="dash", annotation_text='Mean', annotation_bgcolor="black",
                  annotation_textangle=90, row=1, col=1)
    fig.add_vline(median, line_width=1, line_dash="dash", annotation_text='Median', annotation_bgcolor="blue",
                  annotation_textangle=90, row=1, col=1)
    fig.add_vline(n50, line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="green",
                  annotation_textangle=90, row=1, col=1)

    xd = log_bins[1:]
    customdata = np.dstack((np.power(10, log_bins)[:-1], np.power(10, log_bins)[1:], log_hist))[0, :, :]
    yd = log_hist
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=customdata,
                         hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata[2]:.0f}<extra></extra>',
                         marker_color='#36a5c7'), row=2, col=1)

    fig.add_vline(np.log10(mean), line_width=1, line_dash="dash", annotation_text='Mean', annotation_bgcolor="black",
                  annotation_textangle=90, row=2, col=1)
    fig.add_vline(np.log10(median), line_width=1, line_dash="dash", annotation_text='Median', annotation_bgcolor="blue",
                  annotation_textangle=90, row=2, col=1)
    fig.add_vline(np.log10(n50), line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="green",
                  annotation_textangle=90, row=2, col=1)
    fig.update_annotations(font=dict(color="white"))

    fig.update_xaxes(
        tickmode='array',
        tickvals=list(range(0, 12)),
        ticktext=['0'] + ['{:,}'.format(10 ** x) for x in range(1, 12)],
        ticks="outside", row=2, col=1)

    fig.update_xaxes(ticks="outside", title_text='Read Length', title_standoff=0, row=1, col=1)
    fig.update_xaxes(ticks="outside", title_text='Read Length (Log scale)', title_standoff=0, row=2, col=1)
    fig.update_yaxes(ticks="outside", title_text='Counts', title_standoff=0)

    # Set font sizes
    fig.update_layout(showlegend=False, autosize=True,
                      font=dict(size=font_size))

    fig.update_annotations(font_size=annotation_size)
    html_obj = fig.to_html(full_html=False)
    fig.write_image(path, engine="kaleido")

    return html_obj


def base_quality(data, path, font_size):
    """
    Save the 'Base quality' plot image.
    """
    xd = np.arange(256)
    yd = np.array(data.base_quality_distribution)
    fig = go.Figure()

    customdata = np.dstack((xd, yd))[0, :, :]
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=customdata,
                         hovertemplate='Base Quality: %{customdata[0]:.0f}<br>Base Counts:%{customdata[1]:.0f}<extra></extra>',
                         marker_color='#36a5c7'))

    fig.update_xaxes(ticks="outside", dtick=10, title_text='Base Quality', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of bases', title_standoff=0)
    fig.update_layout(font=dict(size=font_size))  # Set font size
    fig.write_image(path, engine="kaleido")

    return fig.to_html(full_html=False)


def read_avg_base_quality(data, path, font_size):
    """
    Save the 'Average base quality' plot image.
    """
    xd = np.arange(256)
    yd = np.array(data.read_average_base_quality_distribution)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=xd, y=yd, marker_color='#36a5c7'))

    fig.update_xaxes(ticks="outside", dtick=10, title_text='Average Base Quality', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of Reads', title_standoff=0)
    fig.update_layout(font=dict(size=font_size))  # Set font size

    fig.write_image(path, engine="kaleido")

    return fig.to_html(full_html=False)


def create_statistics_table(module_output, table_title="Basic Statistics"):
    lrst_global.plot_filenames["basic_st"] = {}
    lrst_global.plot_filenames["basic_st"]['file'] = ""
    lrst_global.plot_filenames["basic_st"]['title'] = "Summary Table"
    lrst_global.plot_filenames["basic_st"]['description'] = table_title

    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
    table_str += int_str_for_format.format("#Total Reads",
                                           module_output.long_read_info.total_num_reads)
    table_str += int_str_for_format.format("#Total Bases",
                                           module_output.long_read_info.total_num_bases)
    table_str += int_str_for_format.format("Longest Read Length",
                                           module_output.long_read_info.longest_read_length)
    table_str += int_str_for_format.format("N50",
                                           module_output.long_read_info.n50_read_length)
    table_str += double_str_for_format.format("GC Content(%)",
                                              module_output.long_read_info.gc_cnt * 100)
    table_str += double_str_for_format.format("Mean Read Length",
                                              module_output.long_read_info.mean_read_length)
    table_str += int_str_for_format.format("Median Read Length",
                                           module_output.long_read_info.median_read_length)
    table_str += "\n</tbody>\n</table>"

    lrst_global.plot_filenames["basic_st"]['detail'] = table_str


def create_base_quality_plots(module_output, para_dict, table_title):
    """
    Generate HTML plots for base and base quality (BAM, FASTQ, FAST5).
    """
    out_path = para_dict["output_folder"]
    get_image_path = lambda x: os.path.join(out_path, lrst_global.plot_filenames[x]['file'])

    # Set the default matplotlib font size
    setDefaultFontSize(12)

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create the statistics table
    create_statistics_table(module_output, table_title)

    # Create basic plots
    basic_data = module_output.long_read_info
    plot_read_length_stats([basic_data], get_image_path('read_length_st'), subtitles=[''])
    plot_base_counts([basic_data], get_image_path('base_st'), subtitles=[''])
    plot_basic_info([basic_data], get_image_path('basic_info'), categories=['Read'])

    # Read length histogram
    length_hist_path = get_image_path('read_length_hist')
    lrst_global.plot_filenames['read_length_hist']['dynamic'] = histogram(basic_data, length_hist_path, font_size)

    # Base quality histogram
    quality_data = module_output.seq_quality_info
    quality_hist_path = get_image_path('base_quality')
    lrst_global.plot_filenames['base_quality']['dynamic'] = base_quality(quality_data, quality_hist_path, font_size)

    # Read quality histogram
    read_quality_dynamic = read_avg_base_quality(quality_data, get_image_path('read_avg_base_quality'), font_size)
    lrst_global.plot_filenames['read_avg_base_quality']['dynamic'] = read_quality_dynamic
