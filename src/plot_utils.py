import os
import logging
import numpy as np
import itertools

import plotly.graph_objs as go
from plotly.subplots import make_subplots


# Return the default image path
def getDefaultImageFolder():
    return 'img/'


# Return the default image suffix
def getDefaultImageSuffix():
    return '.png'


# Return a dictionary of default plot filenames
def getDefaultPlotFilenames():
    default_image_path = getDefaultImageFolder()
    default_image_suf = getDefaultImageSuffix()

    plot_filenames = {  # for fq/fa
        "read_length_distr": {'file': default_image_path + "read_length_distr" + default_image_suf,
                              'title': "Read Length", 'description': "Read Length Distribution"},  # for bam
        "read_alignments_bar": {'file': default_image_path + "map_st" + default_image_suf, 'title': "Map Information",
                   'description': "Read Mapping Statistics"},
        "base_alignments_bar": {'file': default_image_path + "err_st" + default_image_suf,
                   'title': "Base Alignment and Error Statistics",
                   'description': "Base Alignment and Error Statistics"},
        "read_length_bar": {'file': default_image_path + "read_length_st" + default_image_suf,
                           'title': "Read Length Statistics", 'description': "Read Length Statistics"},
        "base_counts": {'file': default_image_path + "base_st" + default_image_suf, 'title': "Base Count Statistics",
                    'description': "Base Count Statistics", 'summary': ""},
        "basic_info": {'file': default_image_path + "basic_info" + default_image_suf, 'title': "Basic Statistics",
                       'description': "Basic Statistics", 'summary': ""},
        "read_length_hist": {'file': default_image_path + "read_length_hist" + default_image_suf,
                             'title': "Read Length Histogram", 'description': "Read Length Histogram", 'summary': ""},

        "base_quality": {'file': default_image_path + "base_quality" + default_image_suf,
                         'title': "Base Quality Histogram", 'description': "Base Quality Histogram"},

        "read_avg_base_quality": {'file': default_image_path + "read_avg_base_quality" + default_image_suf,
                                  'title': "Read Base Quality Histogram", 'description': "Read Base Quality Histogram"},

        "pos_quality": {'file': default_image_path + "pos_quality" + default_image_suf,
                        'title': "Base Position Quality", 'description': "Base Position Quality"},
    }

    return plot_filenames


def fmt(x):
    """Format numbers for plots."""
    format_x = "{:,}".format(round(x))

    return format_x


def wrap(s):
    l = s.split(' ')
    split = list(zip(*[iter(l)] * 3))
    if len(l) % 3:
        split.append(tuple(l[-(len(l) % 3):]))
    return '\n'.join([' '.join(x) for x in split])


def plot_read_length_stats(output_data):
    # Define the three categories
    category = ['N50', 'Mean', 'Median']

    # Create a bar trace for each type of read length statistic
    bar_titles = ['All Reads', 'Mapped Reads', 'Unmapped Reads']
    data_objects = [output_data.long_read_info, output_data.mapped_long_read_info, output_data.unmapped_long_read_info]
    all_traces = []
    for i in range(3):
        plot_title = bar_titles[i]
        data = data_objects[i]
        values = [data.n50_read_length, data.mean_read_length, data.median_read_length]
        trace = go.Bar(x=category, y=values, name=plot_title)
        all_traces.append(trace)

    # Create the layout
    layout = go.Layout(title='Read Length Statistics', xaxis=dict(title='Statistics'), yaxis=dict(title='Length (bp)'), barmode='group')

    # Create the figure and add the traces
    fig = go.Figure(data=all_traces, layout=layout)

    # Generate the HTML
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj

def plot_base_counts(output_data, filetype):
    # Define the five categories
    category = ['A', 'C', 'G', 'T/U', 'N']

    # Create a bar trace for each type of data
    all_traces = []
    if filetype == 'BAM':
        bar_titles = ['All Reads', 'Mapped Reads', 'Unmapped Reads']
        data_objects = [output_data.long_read_info, output_data.mapped_long_read_info, output_data.unmapped_long_read_info]
        for i in range(3):
            plot_title = bar_titles[i]
            data = data_objects[i]
            values = [data.total_a_cnt, data.total_c_cnt, data.total_g_cnt, data.total_tu_cnt, data.total_n_cnt]
            trace = go.Bar(x=category, y=values, name=plot_title)
            all_traces.append(trace)
    else:
        bar_titles = ['All Reads']
        data_objects = [output_data.long_read_info]
        for i in range(1):
            plot_title = bar_titles[i]
            data = data_objects[i]
            values = [data.total_a_cnt, data.total_c_cnt, data.total_g_cnt, data.total_tu_cnt, data.total_n_cnt]
            trace = go.Bar(x=category, y=values, name=plot_title)
            all_traces.append(trace)

    # Create the layout
    layout = go.Layout(title='Base Counts', xaxis=dict(title='Base'), yaxis=dict(title='Counts'), barmode='group')

    # Create the figure and add the traces
    fig = go.Figure(data=all_traces, layout=layout)

    # Generate the HTML
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj


def plot_basic_info(output_data, file_type):
    # Define the four categories
    category = ['Number of Reads', 'Number of Bases', 'Longest Read', 'GC Content']

    if file_type == 'BAM':

        # Create a bar trace for each type of data
        bar_titles = ['All Reads', 'Mapped Reads', 'Unmapped Reads']
        data_objects = [output_data.long_read_info, output_data.mapped_long_read_info, output_data.unmapped_long_read_info]
        all_traces = []

        # Create subplots for each category
        fig = make_subplots(rows=2, cols=2, subplot_titles=("Number of Reads", "Number of Bases", "Longest Read", "GC Content"), horizontal_spacing=0.3, vertical_spacing=0.2)

        # Add traces for each category
        key_list = ['total_num_reads', 'total_num_bases', 'longest_read_length', 'gc_cnt']
        for i in range(4):
            # Get the data for this category
            key_name = key_list[i]

            # Add the traces for each type of data
            data = [getattr(data_objects[0], key_name), getattr(data_objects[1], key_name), getattr(data_objects[2], key_name)]

            # Create the trace
            trace = go.Bar(x=data, y=bar_titles, orientation='h')

            # Add the trace to the figure
            fig.add_trace(trace, row=(i // 2) + 1, col=(i % 2) + 1)

    else:
        bar_titles = ['All Reads']
        data_objects = [output_data.long_read_info]
        all_traces = []

        # Create subplots for each category
        fig = make_subplots(rows=2, cols=2, subplot_titles=("Number of Reads", "Number of Bases", "Longest Read", "GC Content"), horizontal_spacing=0.3, vertical_spacing=0.2)

        # Add traces for each category
        key_list = ['total_num_reads', 'total_num_bases', 'longest_read_length', 'gc_cnt']
        for i in range(4):


    # Create the layout with different y-axis for each category
    layout = go.Layout(title='', xaxis=dict(title='Statistic'), yaxis=dict(title='Value'), barmode='group', showlegend=False)

    fig.update_layout(showlegend=False)

    # Customize the y-axes titles
    # fig.update_yaxes(title_text='Y-axis 1', row=1, col=1)
    # fig.update_yaxes(title_text='Y-axis 2', row=1, col=2)

    # Generate the HTML
    html_obj = fig.to_html(full_html=False, default_height=800, default_width=1200)

    return html_obj


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
    fig.write_image(path, engine="auto", default_height=500, default_width=700)

    return html_obj


def read_lengths_histogram(data, path, font_size):
    """Plot the read length histograms."""
    annotation_size = 10  # Annotation font size
    mean, median, n50 = data.mean_read_length, data.median_read_length, data.n50_read_length

    # Read the read lengths array in float64 format
    read_lengths = np.array(data.read_lengths, dtype=np.float64)

    # Calculate a histogram of read lengths
    hist, edges = np.histogram(read_lengths, bins=10)

    # Create a figure with two subplots
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=("Read Length Histogram", "Log Read Length Histogram"), vertical_spacing=0.3)
    fig.update_layout(showlegend=False, autosize=False)

    customdata = np.dstack((edges[:-1], edges[1:], hist))[0, :, :]
    fig.add_trace(go.Bar(x=edges, y=hist, customdata=customdata,
                         hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata['
                                       '2]:.0f}<extra></extra>',
                         marker_color='#36a5c7'), row=1, col=1)

    fig.add_vline(mean, line_width=1, line_dash="dash", annotation_text='Mean', annotation_bgcolor="black",
                  annotation_textangle=90, row=1, col=1)
    fig.add_vline(median, line_width=1, line_dash="dash", annotation_text='Median', annotation_bgcolor="blue",
                  annotation_textangle=90, row=1, col=1)
    fig.add_vline(n50, line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="green",
                  annotation_textangle=90, row=1, col=1)

    # Log histogram
    # Get the log10 histogram of read lengths
    read_lengths_log = np.log10(read_lengths, out=np.zeros_like(read_lengths), where=(read_lengths != 0))
    log_hist, log_edges = np.histogram(read_lengths_log, bins=len(edges))

    xd = log_edges
    customdata = np.dstack((np.power(10, log_edges)[:-1], np.power(10, log_edges)[1:], log_hist))[0, :, :]
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

    # Set tick value range for the log scale
    tick_vals = list(range(0, 5))
    fig.update_xaxes(
        range=[0, 5],
        tickmode='array',
        tickvals=tick_vals,
        ticktext=['{:,}'.format(10 ** x) for x in tick_vals],
        ticks="outside", title_text='Read Length (Log Scale)', title_standoff=0, row=2, col=1)

    fig.update_xaxes(ticks="outside", title_text='Read Length', title_standoff=0, row=1, col=1)
    fig.update_yaxes(ticks="outside", title_text='Counts', title_standoff=0)

    # Set font sizes
    fig.update_layout(font=dict(size=font_size), autosize=True)

    fig.update_annotations(font_size=annotation_size)
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)
    fig.write_image(path, engine="auto")

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
                         hovertemplate='Base Quality: %{customdata[0]:.0f}<br>Base Counts:%{customdata['
                                       '1]:.0f}<extra></extra>',
                         marker_color='#36a5c7'))

    fig.update_xaxes(ticks="outside", dtick=10, title_text='Base Quality', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of bases', title_standoff=0)
    fig.update_layout(font=dict(size=font_size))  # Set font size
    fig.write_image(path, engine="auto")

    return fig.to_html(full_html=False, default_height=500, default_width=700)


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

    fig.write_image(path, engine="auto")

    return fig.to_html(full_html=False, default_height=500, default_width=700)


def create_statistics_table(output_data, plot_filepaths, table_title="Basic Statistics"):
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Summary Table"
    plot_filepaths["basic_st"]['description'] = table_title
    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
    table_str += int_str_for_format.format("#Total Reads",
                                           output_data.long_read_info.total_num_reads)
    table_str += int_str_for_format.format("#Total Bases",
                                           output_data.long_read_info.total_num_bases)
    table_str += int_str_for_format.format("Longest Read Length",
                                           output_data.long_read_info.longest_read_length)
    table_str += int_str_for_format.format("N50",
                                           output_data.long_read_info.n50_read_length)
    table_str += double_str_for_format.format("GC Content(%)",
                                              output_data.long_read_info.gc_cnt * 100)
    table_str += double_str_for_format.format("Mean Read Length",
                                              output_data.long_read_info.mean_read_length)
    table_str += int_str_for_format.format("Median Read Length",
                                           output_data.long_read_info.median_read_length)
    table_str += "\n</tbody>\n</table>"

    plot_filepaths["basic_st"]['detail'] = table_str

    return plot_filepaths


def create_base_quality_plots(output_data, para_dict, table_title):
    """
    Generate HTML plots for base and base quality (BAM, FASTQ, FAST5).
    """
    out_path = para_dict["output_folder"]
    plot_filepaths = getDefaultPlotFilenames()
    get_image_path = lambda x: os.path.join(out_path, plot_filepaths[x]['file'])

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create the statistics table
    plot_filepaths = create_statistics_table(output_data, plot_filepaths, table_title)

    # Create basic plots
    basic_data = output_data.long_read_info
    plot_read_length_stats([basic_data], get_image_path('read_length_st'), subtitles=[''])
    plot_base_counts([basic_data], get_image_path('base_st'), subtitles=[''])
    plot_basic_info([basic_data], get_image_path('basic_info'), categories=['Read'])

    # Read length histogram
    length_hist_path = get_image_path('read_length_hist')
    plot_filepaths['read_length_hist']['dynamic'] = histogram(basic_data, length_hist_path, font_size)

    # Base quality histogram
    quality_data = output_data.seq_quality_info
    quality_hist_path = get_image_path('base_quality')
    plot_filepaths['base_quality']['dynamic'] = base_quality(quality_data, quality_hist_path, font_size)

    # Read quality histogram
    read_quality_dynamic = read_avg_base_quality(quality_data, get_image_path('read_avg_base_quality'), font_size)
    plot_filepaths['read_avg_base_quality']['dynamic'] = read_quality_dynamic

    return plot_filepaths


def plot_fastq(output_data, para_dict):
    """
    Generate HTML plots for FASTQ.
    """
    plot_filepaths = plot(output_data, para_dict, 'FASTQ')

    return plot_filepaths


def plot(output_data, para_dict, file_type):
    out_path = para_dict["output_folder"]
    plot_filepaths = getDefaultPlotFilenames()
    get_image_path = lambda x: os.path.join(out_path, plot_filepaths[x]['file'])

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create the summary table
    create_summary_table(output_data, plot_filepaths, file_type)

    # Generate plots
    plot_filepaths['base_counts']['dynamic'] = plot_base_counts(output_data, file_type)
    plot_filepaths['basic_info']['dynamic'] = plot_basic_info(output_data, file_type)
    plot_filepaths['read_length_hist']['dynamic'] = read_lengths_histogram(output_data.long_read_info,
                                                                           get_image_path('read_length_hist'),
                                                                           font_size)

    plot_filepaths['base_quality']['dynamic'] = base_quality(output_data.seq_quality_info,
                                                             get_image_path('base_quality'), font_size)

    plot_filepaths['read_length_bar']['dynamic'] = plot_read_length_stats(output_data)

    if file_type == 'BAM':
        plot_filepaths['read_alignments_bar']['dynamic'] = plot_alignment_numbers(output_data,
                                                                                  get_image_path('read_alignments_bar'))
        plot_filepaths['base_alignments_bar']['dynamic'] = plot_errors(output_data, get_image_path('base_alignments_bar'))

    return plot_filepaths


def create_summary_table(output_data, plot_filepaths, file_type):
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Summary Table"
    plot_filepaths["basic_st"]['description'] = "Basic statistics of the {} file.".format(file_type)

    if file_type == 'BAM':
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Mapped</th><th>Unmapped</th><th>All</th></tr>\n" \
                    "</thead> "
        table_str += "\n<tbody>"
        int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:," \
                             "d}</td><td style=\"text-align:right\">{:,d}</td></tr> "
        double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td " \
                                "style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr> "
        table_str += int_str_for_format.format("#Total Reads", output_data.mapped_long_read_info.total_num_reads,
                                               output_data.unmapped_long_read_info.total_num_reads,
                                               output_data.long_read_info.total_num_reads)
        table_str += int_str_for_format.format("#Total Bases",
                                               output_data.mapped_long_read_info.total_num_bases,
                                               output_data.unmapped_long_read_info.total_num_bases,
                                               output_data.long_read_info.total_num_bases)
        table_str += int_str_for_format.format("Longest Read Length",
                                               output_data.mapped_long_read_info.longest_read_length,
                                               output_data.unmapped_long_read_info.longest_read_length,
                                               output_data.long_read_info.longest_read_length)
        table_str += int_str_for_format.format("N50",
                                               output_data.mapped_long_read_info.n50_read_length,
                                               output_data.unmapped_long_read_info.n50_read_length,
                                               output_data.long_read_info.n50_read_length)
        table_str += double_str_for_format.format("GC Content(%)",
                                                  output_data.mapped_long_read_info.gc_cnt * 100,
                                                  output_data.unmapped_long_read_info.gc_cnt * 100,
                                                  output_data.long_read_info.gc_cnt * 100)
        table_str += double_str_for_format.format("Mean Read Length",
                                                  output_data.mapped_long_read_info.mean_read_length,
                                                  output_data.unmapped_long_read_info.mean_read_length,
                                                  output_data.long_read_info.mean_read_length)
        table_str += int_str_for_format.format("Median Read Length",
                                               output_data.mapped_long_read_info.median_read_length,
                                               output_data.unmapped_long_read_info.median_read_length,
                                               output_data.long_read_info.median_read_length)
    else:
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
        table_str += "\n<tbody>"
        int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
        double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
        table_str += int_str_for_format.format("#Total Reads",
                                               output_data.long_read_info.total_num_reads)
        table_str += int_str_for_format.format("#Total Bases",
                                               output_data.long_read_info.total_num_bases)
        table_str += int_str_for_format.format("Longest Read Length",
                                               output_data.long_read_info.longest_read_length)
        table_str += int_str_for_format.format("N50",
                                               output_data.long_read_info.n50_read_length)
        table_str += double_str_for_format.format("GC Content(%)",
                                                  output_data.long_read_info.gc_cnt * 100)
        table_str += double_str_for_format.format("Mean Read Length",
                                                  output_data.long_read_info.mean_read_length)
        table_str += int_str_for_format.format("Median Read Length",
                                               output_data.long_read_info.median_read_length)
        
    table_str += "\n</tbody>\n</table>"
    plot_filepaths["basic_st"]['detail'] = table_str


def plot_alignment_numbers(data, path):
    category = ['Primary Alignments', 'Supplementary Alignments', 'Secondary Alignments',
                'Reads with Supplementary Alignments', 'Reads with Secondary Alignments',
                'Reads with Secondary and Supplementary Alignments', 'Forward Alignments', 'Reverse Alignments']
    category = [wrap(x) for x in category]

    # Create a horizontally aligned bar plot trace from the data using plotly
    trace = go.Bar(x=[data.num_primary_alignment, data.num_supplementary_alignment, data.num_secondary_alignment,
                      data.num_reads_with_supplementary_alignment, data.num_reads_with_secondary_alignment,
                      data.num_reads_with_both_secondary_supplementary_alignment, data.forward_alignment,
                      data.reverse_alignment], y=category, orientation='h')

    # Create the layout for the plot
    layout = go.Layout(title=go.layout.Title(text=""),
                       xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Counts")),
                       yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Generate the HTML object for the plot
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=1000)

    return html_obj


def plot_errors(output_data, path):
    category = ['Matched Bases', 'Mismatched Bases', 'Inserted Bases', 'Deleted Bases', 'Clipped Bases']
    category = [wrap(x) for x in category]

    # Create a horizontally aligned bar plot trace from the data using plotly
    trace = go.Bar(x=[output_data.num_matched_bases, output_data.num_mismatched_bases, output_data.num_ins_bases,
                      output_data.num_del_bases, output_data.num_clip_bases], y=category, orientation='h')

    # Create the layout for the plot
    layout = go.Layout(title=go.layout.Title(text=""),
                       xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Counts")),
                       yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Generate the HTML object for the plot
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj

