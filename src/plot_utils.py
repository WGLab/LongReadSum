import os
import numpy as np
import csv
from random import sample

import plotly.graph_objs as go
from plotly.subplots import make_subplots

if __name__ == 'src.plot_utils':
    from lib import lrst  # For debugging
else:
    import lrst

# Set up logging
import logging
logging.basicConfig(level=logging.INFO)

# Constants
MAX_BASE_QUALITY = 100
MAX_READ_QUALITY = 100
PLOT_FONT_SIZE = 16


def getDefaultPlotFilenames():
    """Create a default HTML plot data structure."""
    plot_filenames = {  # for fq/fa
        "read_length_distr": {'title': "Read Length", 'description': "Read Length Distribution"},  # for bam
        "read_alignments_bar": {'title': "Read Alignments",
                   'description': "Read Alignments"},
        "base_alignments_bar": {'title': "Base Alignment and Error",
                   'description': "Base Alignment and Error"},
        "read_length_bar": {'title': "Read Length Statistics", 'description': "Read Length Statistics"},
        "base_counts": {'title': "Base Counts",
                    'description': "Base Counts", 'summary': ""},
        "read_length_hist": {'title': "Read Length Histogram", 'description': "Read Length Histogram", 'summary': ""},
        
        "gc_content_hist": {'title': "GC Content Histogram", 'description': "GC Content Histogram", 'summary': ""},

        "read_length_mod_rates": {'title': "Read Length vs. Modification Rates", 'description': "Read Length vs. Modification Rates", 'summary': ""},

        "base_quality": {'title': "Base Quality Histogram", 'description': "Base Quality Histogram"},

        "read_avg_base_quality": {'title': "Read Base Quality Histogram", 'description': "Read Base Quality Histogram"},

        "pos_quality": {'title': "Base Position Quality", 'description': "Base Position Quality"},
        "ont_signal": {'title': "ONT Signal", 'description': "ONT Signal"},
    }

    return plot_filenames


def wrap(label):
    """Wrap the label text."""
    # First split the string into a list of words
    words = label.split(' ')

    # Then join the words back together with <br> tags if the total length is greater than 30
    new_label = ''
    current_length = 0
    max_length = 30
    for word in words:
        if current_length > max_length:
            new_label += '<br>'
            current_length = 0

        new_label += word + ' '
        current_length = len(new_label)

    # Remove the last space
    new_label = new_label[:-1]

    return new_label


def plot_read_length_stats(output_data, file_type, plot_filepaths):
    """Plot the read length statistics."""

    # Define the three categories
    category = ['N50', 'Mean', 'Median']
    all_traces = []
    error_flag = False

    if file_type == 'BAM':
        # Create a bar trace for each type of read length statistic
        bar_titles = ['All Reads', 'Mapped Reads', 'Unmapped Reads']
        data_objects = [output_data.long_read_info, output_data.mapped_long_read_info, output_data.unmapped_long_read_info]
        for i in range(3):
            plot_title = bar_titles[i]
            data = data_objects[i]
            values = [data.n50_read_length, data.mean_read_length, data.median_read_length]
            trace = go.Bar(x=category, y=values, name=plot_title)
            all_traces.append(trace)

            # Set the error flag if any of the values are zero (except for unmapped reads)
            if i != 2 and (values[0] == 0 or values[1] == 0 or values[2] == 0):
                error_flag = True

    elif file_type == 'SeqTxt':
        # Create a bar trace for each type of read length statistic
        bar_titles = ['All Reads', 'Passed Reads', 'Failed Reads']
        data_objects = [output_data.all_long_read_info.long_read_info, output_data.passed_long_read_info.long_read_info, output_data.failed_long_read_info.long_read_info]
        for i in range(3):
            plot_title = bar_titles[i]
            data = data_objects[i]
            values = [data.n50_read_length, data.mean_read_length, data.median_read_length]
            trace = go.Bar(x=category, y=values, name=plot_title)
            all_traces.append(trace)

            # Set the error flag if any of the values are zero (except for failed reads)
            if i != 2 and (values[0] == 0 or values[1] == 0 or values[2] == 0):
                error_flag = True

    else:
        # Get the data for all reads
        key_list = ['n50_read_length', 'mean_read_length', 'median_read_length']

        # Create a bar trace
        bar_title = 'All Reads'
        data = output_data.long_read_info
        values = [getattr(data, key_name) for key_name in key_list]
        trace = go.Bar(x=category, y=values, name=bar_title)
        all_traces.append(trace)

        # Set the error flag if any of the values are zero
        if values[0] == 0 or values[1] == 0 or values[2] == 0:
            error_flag = True


    # Create the layout
    layout = go.Layout(title='', xaxis=dict(title='Statistics'), yaxis=dict(title='Length (bp)'), barmode='group', font=dict(size=PLOT_FONT_SIZE))

    # Create the figure and add the traces
    fig = go.Figure(data=all_traces, layout=layout)

    # Generate the HTML
    # html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)
    plot_filepaths['read_length_bar']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=700)

    # Set the error flag
    plot_filepaths['read_length_bar']['error_flag'] = error_flag


def plot_base_counts(output_data, filetype, plot_filepaths):
    """Plot overall base counts for the reads."""

    # Create a bar trace for each base
    error_flag = False
    category = ['A', 'C', 'G', 'T/U', 'N']
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

            # Set the error flag if the N count is greater than 10% or the A, C,
            # G, or T/U counts are zero (except for unmapped reads)
            if i != 2:
                if data.total_num_bases == 0:
                    error_flag = True
                elif data.total_n_cnt / data.total_num_bases > 0.1:
                    error_flag = True
                elif data.total_a_cnt == 0 or data.total_c_cnt == 0 or data.total_g_cnt == 0 or data.total_tu_cnt == 0:
                    error_flag = True

    elif filetype == 'SeqTxt':
        bar_titles = ['All Reads', 'Passed Reads', 'Failed Reads']
        data_objects = [output_data.all_long_read_info.long_read_info, output_data.passed_long_read_info.long_read_info, output_data.failed_long_read_info.long_read_info]
        for i in range(3):
            plot_title = bar_titles[i]
            data = data_objects[i]
            values = [data.total_a_cnt, data.total_c_cnt, data.total_g_cnt, data.total_tu_cnt, data.total_n_cnt]
            trace = go.Bar(x=category, y=values, name=plot_title)
            all_traces.append(trace)

            # Set the error flag if the N count is greater than 10% or the A, C,
            # G, or T/U counts are zero
            if data.total_num_bases == 0:
                error_flag = True
            elif data.total_n_cnt / data.total_num_bases > 0.1:
                error_flag = True
            elif data.total_a_cnt == 0 or data.total_c_cnt == 0 or data.total_g_cnt == 0 or data.total_tu_cnt == 0:
                error_flag = True

    else:
        plot_title = 'All Reads'
        data = output_data.long_read_info
        values = [data.total_a_cnt, data.total_c_cnt, data.total_g_cnt, data.total_tu_cnt, data.total_n_cnt]
        trace = go.Bar(x=category, y=values, name=plot_title)
        all_traces.append(trace)

        # Set the error flag if the N count is greater than 10% or the A, C,
        # G, or T/U counts are zero
        if data.total_num_bases == 0:
            error_flag = True
        elif data.total_n_cnt / data.total_num_bases > 0.1:
            error_flag = True
        elif data.total_a_cnt == 0 or data.total_c_cnt == 0 or data.total_g_cnt == 0 or data.total_tu_cnt == 0:
            error_flag = True

    # Create the figure and add the traces
    layout = go.Layout(title='', xaxis=dict(title='Base'), yaxis=dict(title='Counts'), barmode='group', font=dict(size=PLOT_FONT_SIZE))
    fig = go.Figure(data=all_traces, layout=layout)

    # Generate the HTML
    plot_filepaths['base_counts']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=700)
    plot_filepaths['base_counts']['error_flag'] = error_flag


def read_lengths_histogram(data, font_size, plot_filepaths):
    """Plot the read length histogram."""
    linear_bin_count = 10
    log_bin_count = 10

    annotation_size = 10  # Annotation font size
    mean, median, n50 = data.mean_read_length, data.median_read_length, data.n50_read_length

    # Read the read lengths array in float64 format
    read_lengths = np.array(data.read_lengths, dtype=np.float64)

    # If there are no read lengths, throw an error
    if len(read_lengths) == 0:
        raise ValueError("No read lengths found")

    # Calculate a histogram of read lengths, but don't center the bins
    # edges = np.linspace(1, np.max(read_lengths), num=bin_count + 1)
    edges = np.linspace(np.min(read_lengths), np.max(read_lengths), num=linear_bin_count + 1)
    hist, _ = np.histogram(read_lengths, bins=edges)

    # Create a figure with two subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Read Length Histogram", "Log Read Length Histogram"), vertical_spacing=0.0)
    linear_col=1
    log_col=2

    linear_bindata = np.dstack((edges[:-1], edges[1:], hist))[0, :, :]
    fig.add_trace(go.Bar(x=edges, y=hist, customdata=linear_bindata,
                         hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata[2]:.0f}<extra></extra>',
                         marker_color='#36a5c7'), row=1, col=linear_col)

    fig.add_vline(mean, line_width=1, line_dash="dash", annotation_text='Mean', annotation_bgcolor="black",
                  annotation_textangle=90, row=1, col=linear_col)
    fig.add_vline(median, line_width=1, line_dash="dash", annotation_text='Median', annotation_bgcolor="blue",
                  annotation_textangle=90, row=1, col=linear_col)
    fig.add_vline(n50, line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="green",
                  annotation_textangle=90, row=1, col=linear_col)

    # Log scale histogram
    read_lengths_log = np.log10(read_lengths, out=np.zeros_like(read_lengths), where=(read_lengths != 0))
    log_edges = np.linspace(0, np.max(read_lengths_log), num=log_bin_count + 1)
    log_hist, _ = np.histogram(read_lengths_log, bins=log_edges)

    xd = log_edges
    log_bindata = np.dstack((np.power(10, log_edges)[:-1], np.power(10, log_edges)[1:], log_hist))[0, :, :]
    yd = log_hist
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=log_bindata,
                         hovertemplate='Length: %{customdata[0]:.0f}-%{customdata[1]:.0f}bp<br>Counts:%{customdata[2]:.0f}<extra></extra>',
                         marker_color='#36a5c7'), row=1, col=log_col)

    fig.add_vline(np.log10(mean), line_width=1, line_dash="dash", annotation_text='Mean', annotation_bgcolor="black",
                  annotation_textangle=90, row=1, col=log_col)
    fig.add_vline(np.log10(median), line_width=1, line_dash="dash", annotation_text='Median', annotation_bgcolor="blue",
                  annotation_textangle=90, row=1, col=log_col)
    fig.add_vline(np.log10(n50), line_width=1, line_dash="dash", annotation_text='N50', annotation_bgcolor="green",
                  annotation_textangle=90, row=1, col=log_col)
    fig.update_annotations(font=dict(color="white"))

    # Set tick value range for the log scale
    tick_vals = log_edges
    tick_labels = []
    for i in range(len(log_bindata)):
        # Format the tick labels to be in kilobases (kb) if the value is greater
        # than 1000, in megabases (Mb) if the value is greater than 1,000,000,
        # and in bases (b) if less than 1000
        left_val = log_bindata[i][0]
        left_val_str = '{:,}Mb'.format(int(left_val / 1000000)) if left_val > 1000000 else '{:,}kb'.format(int(left_val / 1000)) if left_val > 1000 else '{:,}bp'.format(int(left_val))

        right_val = log_bindata[i][1]
        right_val_str = '{:,}Mb'.format(int(right_val / 1000000)) if right_val > 1000000 else '{:,}kb'.format(int(right_val / 1000)) if right_val > 1000 else '{:,}bp'.format(int(right_val))

        tick_labels.append('{}-{}'.format(left_val_str, right_val_str))

    fig.update_xaxes(ticks="outside", title_text='Read Length (Log Scale)', title_standoff=0, row=1, col=log_col, tickvals=tick_vals, ticktext=tick_labels, tickangle=45)
    tick_vals = edges
    
    # Format the tick labels to be the range of the bin centers
    tick_labels = []
    for i in range(len(linear_bindata)):
        # Format the tick labels to be in kilobases (kb) if the value is greater
        # than 1000, in megabases (Mb) if the value is greater than 1,000,000,
        # and in bases (b) if less than 1000
        left_val = linear_bindata[i][0]
        left_val_str = '{:,}Mb'.format(int(left_val / 1000000)) if left_val > 1000000 else '{:,}kb'.format(int(left_val / 1000)) if left_val > 1000 else '{:,}bp'.format(int(left_val))

        right_val = linear_bindata[i][1]
        right_val_str = '{:,}Mb'.format(int(right_val / 1000000)) if right_val > 1000000 else '{:,}kb'.format(int(right_val / 1000)) if right_val > 1000 else '{:,}bp'.format(int(right_val))

        tick_labels.append('{}-{}'.format(left_val_str, right_val_str))
        
    linear_col=1
    fig.update_xaxes(ticks="outside", title_text='Read Length', title_standoff=0, row=1, col=linear_col, tickvals=tick_vals, ticktext=tick_labels, tickangle=45)
    fig.update_yaxes(ticks="outside", title_text='Counts', title_standoff=0)

    # Update the layout
    fig.update_layout(showlegend=False, autosize=True, font=dict(size=PLOT_FONT_SIZE))
    fig.update_annotations(font_size=annotation_size)

    # Generate the HTML
    plot_filepaths['read_length_hist']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=1200)
                           

def read_gc_content_histogram(data, font_size, plot_filepaths):
    """Plot the per-read GC content histogram."""
    bin_size = 1
    gc_content = np.array(data.read_gc_content_count)

    # Calculate the percentage of reads with a GC content of <30%
    gc_content_below_30 = np.sum(gc_content[:30])
    logging.info("[TEST] Percentage of reads with GC content <30%: {}".format(gc_content_below_30 / np.sum(gc_content)))

    # Calculate the percentage of reads with a GC content of >70%
    gc_content_above_70 = np.sum(gc_content[70:])
    logging.info("[TEST] Percentage of reads with GC content >70%: {}".format(gc_content_above_70 / np.sum(gc_content)))

    # Calculate the percentage of reads with a GC content of <20%
    gc_content_below_20 = np.sum(gc_content[:20])
    logging.info("[TEST] Percentage of reads with GC content <20%: {}".format(gc_content_below_20 / np.sum(gc_content)))

    # Calculate the percentage of reads with a GC content of >60%
    gc_content_above_60 = np.sum(gc_content[60:])
    logging.info("[TEST] Percentage of reads with GC content >60%: {}".format(gc_content_above_60 / np.sum(gc_content)))

    # Set the error flag if the GC content is below 20% for more than 10% of the
    # reads
    error_flag = False
    if np.sum(gc_content) == 0:
        error_flag = True
    elif np.sum(gc_content[:20]) / np.sum(gc_content) > 0.1:
        error_flag = True

    # Bin the GC content if the bin size is greater than 1
    if bin_size > 1:
        gc_content = np.array([np.sum(gc_content[i:i + bin_size]) for i in range(0, 101, bin_size)])

    gc_content_bins = [i for i in range(0, 101, bin_size)]

    # Generate hover text for each bin
    hover_text = []
    if bin_size > 1:
        for i in range(len(gc_content_bins)):
            hover_text.append('GC content: {}-{}%<br>Counts: {}'.format(gc_content_bins[i], gc_content_bins[i] + bin_size, gc_content[i]))
    else:
        for i in range(len(gc_content_bins)):
            hover_text.append('GC content: {}%<br>Counts: {}'.format(gc_content_bins[i], gc_content[i]))

    # Set the X values to be the center of the bins
    if bin_size > 1:
        x_values = [gc_content_bins[i] + bin_size / 2 for i in range(len(gc_content_bins))]
    else:
        x_values = gc_content_bins

    # Create the figure
    fig = go.Figure()
    fig.add_trace(go.Bar(x=x_values, y=gc_content, marker_color='#36a5c7', hovertext=hover_text, hoverinfo='text'))

    # Update the layout
    fig.update_xaxes(ticks="outside", dtick=10, title_text='GC Content (%)', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of Reads', title_standoff=0)
    fig.update_layout(font=dict(size=PLOT_FONT_SIZE))  # Set font size

    plot_filepaths['gc_content_hist']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=700)
    plot_filepaths['gc_content_hist']['error_flag'] = error_flag


def base_quality(data, font_size, plot_filepaths):
    """Plot the base quality distribution."""
    xd = np.arange(MAX_BASE_QUALITY)
    yd = np.array(data.base_quality_distribution)
    fig = go.Figure()

    customdata = np.dstack((xd, yd))[0, :, :]
    fig.add_trace(go.Bar(x=xd, y=yd, customdata=customdata,
                         hovertemplate='Base Quality: %{customdata[0]:.0f}<br>Base Counts:%{customdata['
                                       '1]:.0f}<extra></extra>',
                         marker_color='#36a5c7'))

    fig.update_xaxes(ticks="outside", dtick=10, title_text='Base Quality', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of bases', title_standoff=0)
    fig.update_layout(font=dict(size=PLOT_FONT_SIZE))  # Set font size

    # return fig.to_html(full_html=False, default_height=500, default_width=700)
    plot_filepaths['base_quality']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=700)

    # Set the error flag if the base quality is below 20 for more than 10% of
    # the bases
    error_flag = False
    if np.sum(yd) == 0:
        error_flag = True
    elif np.sum(yd[:20]) / np.sum(yd) > 0.1:
        error_flag = True

    plot_filepaths['base_quality']['error_flag'] = error_flag


def read_avg_base_quality(data, font_size, plot_filepaths):
    """Plot the read average base quality distribution."""
    xd = np.arange(MAX_READ_QUALITY)
    yd = np.array(data.read_average_base_quality_distribution)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=xd, y=yd, marker_color='#36a5c7'))

    fig.update_xaxes(ticks="outside", dtick=10, title_text='Average Base Quality', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of Reads', title_standoff=0)
    fig.update_layout(font=dict(size=PLOT_FONT_SIZE))  # Set font size

    # return fig.to_html(full_html=False, default_height=500, default_width=700)
    plot_filepaths['read_avg_base_quality']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=700)

    # Set the error flag if the average base quality is below 20 for more than
    # 10% of the reads
    error_flag = False
    if np.sum(yd) == 0:
        error_flag = True
    elif np.sum(yd[:20]) / np.sum(yd) > 0.1:
        error_flag = True

    plot_filepaths['read_avg_base_quality']['error_flag'] = error_flag


def plot_base_modifications(base_modifications):
    """Plot the base modifications per location."""
    # Get the modification types
    modification_types = list(base_modifications.keys())

    # Create the figure
    fig = go.Figure()

    # Add a trace for each modification type
    for mod_type in modification_types:
        # Get the modification data
        mod_data = base_modifications[mod_type]

        # Create the trace
        trace = go.Scattergl(x=mod_data['positions'], y=mod_data['counts'], mode='markers', name=mod_type)

        # Add the trace to the figure
        fig.add_trace(trace)

    # Update the layout
    fig.update_layout(title='Base Modifications', xaxis_title='Position', yaxis_title='Counts', showlegend=True, font=dict(size=PLOT_FONT_SIZE))

    # Generate the HTML
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj


def plot(output_data, para_dict, file_type):
    """Generate the plots for the output data."""
    logging.info("Generating plots for file type: {}".format(file_type))
    plot_filepaths = getDefaultPlotFilenames()
    font_size = 14  # Font size for the plots
    create_summary_table(output_data, plot_filepaths, file_type)  # Create the summary table

    # Modified base table and plots
    try:
        para_dict["mod"]
    except KeyError:
        para_dict["mod"] = False

    if file_type == 'BAM' and para_dict["mod"]:
        # Output file for the read length vs. modification rates plot
        output_folder = para_dict["output_folder"]
        read_length_mod_rate_file = os.path.join(output_folder, 'read_length_hist.png')
        plot_filepaths['read_length_mod_rates']['file'] = read_length_mod_rate_file

        # Generate the modified base table and read length vs. modification rates plot
        base_modification_threshold = para_dict["modprob"]
        create_modified_base_table(output_data, plot_filepaths, base_modification_threshold)
        if 'base_mods' not in plot_filepaths:
            logging.warning("WARNING: Modified base table not created")

    # Create the TIN table if available
    try:
        para_dict["genebed"]
    except KeyError:
        para_dict["genebed"] = ""
        
    if file_type == 'BAM' and para_dict["genebed"] != "":
        input_files = para_dict["input_files"]
        create_tin_table(output_data, input_files, plot_filepaths)

        # Check if the TIN table is available
        if 'tin' in plot_filepaths:
            logging.info("SUCCESS: TIN table created")
        else:
            logging.warning("WARNING: TIN table not created")

    plot_base_counts(output_data, file_type, plot_filepaths)

    # Read length histogram
    if file_type == 'SeqTxt':
        long_read_data = output_data.all_long_read_info.long_read_info
    else:
        long_read_data = output_data.long_read_info

    if file_type != 'FAST5s':
        read_lengths_histogram(long_read_data, font_size, plot_filepaths)
        plot_read_length_stats(output_data, file_type, plot_filepaths)

    # GC content histogram
    if file_type == 'BAM':
        read_gc_content_histogram(output_data.mapped_long_read_info, font_size, plot_filepaths)
    elif file_type == 'SeqTxt':
        read_gc_content_histogram(output_data.passed_long_read_info.long_read_info, font_size, plot_filepaths)
    elif file_type == 'FASTQ' or file_type == 'FASTA':
        read_gc_content_histogram(output_data.long_read_info, font_size, plot_filepaths)

    # Base quality histogram
    if file_type != 'FASTA' and file_type != 'FAST5s' and file_type != 'SeqTxt':
        seq_quality_info = output_data.seq_quality_info

        # Base quality histogram
        base_quality(seq_quality_info, font_size, plot_filepaths)
        
    # Read average base quality histogram
    if file_type == 'FASTQ':
        read_avg_base_quality(seq_quality_info, font_size, plot_filepaths)

    if file_type == 'BAM':
        # Plot read alignment QC
        plot_alignment_numbers(output_data, plot_filepaths)
        
        # Plot base alignment and error QC
        plot_errors(output_data, plot_filepaths)
        
    elif file_type == 'FAST5s':
        plot_filepaths['ont_signal']['dynamic'] = plot_signal(output_data, para_dict)

    return plot_filepaths

def plot_pod5(pod5_output, para_dict, bam_output=None):
    """Plot the ONT POD5 signal data for a random sample of reads."""
    out_path = para_dict["output_folder"]
    plot_filepaths = getDefaultPlotFilenames()

    # Create the summary table
    create_pod5_table(pod5_output, plot_filepaths)

    # Generate the signal plots
    marker_size = 10
    read_count_max = para_dict["read_count"]

    read_count = len(pod5_output.keys())
    logging.info("Plotting signal data for {} reads".format(read_count))

    # Randomly sample a small set of reads if it is a large dataset
    read_sample_size = min(read_count_max, read_count)
    unsampled_indices = list(range(0, read_sample_size))
    read_indices = sample(unsampled_indices, read_sample_size)

    if read_sample_size < read_count:
        logging.info("Randomly sampling {} reads from the total of {} reads".format(read_sample_size, read_count))
    else:
        logging.info("Plotting signal data for all {} reads".format(read_count))


    # Plot the reads
    output_html_plots = {}
    for read_index in read_indices:
        # Create the figure
        fig = go.Figure()

        # Get the read data
        nth_read_name = list(pod5_output.keys())[read_index]
        nth_read_data = pod5_output[nth_read_name]['signal']
        signal_length = len(nth_read_data)
        logging.info("Signal data count for read {}: {}".format(nth_read_name, signal_length))
        nth_read_mean = pod5_output[nth_read_name]['mean']
        nth_read_std = pod5_output[nth_read_name]['std']
        nth_read_median = pod5_output[nth_read_name]['median']
        nth_read_skewness = pod5_output[nth_read_name]['skewness']
        nth_read_kurtosis = pod5_output[nth_read_name]['kurtosis']

        # Set up the output CSV
        csv_qc_filepath = os.path.join(out_path, nth_read_name + '_QC.csv')
        qc_file = open(csv_qc_filepath, 'w', encoding='utf-8')
        qc_writer = csv.writer(qc_file)
        qc_writer.writerow(["Raw_Signal", "Length", "Mean", "Median", "StdDev", "PearsonSkewnessCoeff", "Kurtosis"])
        
        # Update CSV
        raw_row = [nth_read_data, signal_length, nth_read_mean, nth_read_median, nth_read_std, nth_read_skewness, nth_read_kurtosis]
        qc_writer.writerow(raw_row)

        # Close CSV
        qc_file.close()

        # Plot the base sequence if available
        if bam_output:
            move_table = bam_output.getReadMoveTable(nth_read_name)
            read_sequence = bam_output.getReadSequence(nth_read_name)
            start_index = bam_output.getReadSequenceStart(nth_read_name)
            end_index = bam_output.getReadSequenceEnd(nth_read_name)

            # Print the first couple of indices from the table.
            # Each index in the move table represents a k-mer move. Thus, for
            # each base, the signal is between two indices in the move table, starting
            # from the first index.
            logging.info("Move table for read {}: {}".format(nth_read_name, move_table[:5]))
            logging.info("Move table range: {}-{}".format(min(move_table), max(move_table)))
            logging.info("Read sequence for read {}: {}".format(nth_read_name, read_sequence[:5]))
            logging.info("Read sequence length for read {}: {}".format(nth_read_name, len(read_sequence)))
            logging.info("Signal data length for read {}: {}".format(nth_read_name, len(move_table)))
            logging.info("Signal interval for read {}: {}-{}".format(nth_read_name, start_index, end_index))

            # Filter the signal data. Use the last index of the move table + 20
            # as the end index, since the signal data can be much longer than the
            # read sequence.
            end_index = max(move_table) + 20
            nth_read_data = nth_read_data[start_index:end_index]
            signal_length = len(nth_read_data)

            # Set up the X tick values
            base_tick_values = move_table

            # Set up the X tick labels
            x_tick_labels = list(read_sequence)

            # Update the plot style
            fig.update_xaxes(title="Base",
                                tickangle=0,
                                tickmode='array',
                                tickvals=base_tick_values,
                                ticktext=x_tick_labels)
        else:
            fig.update_xaxes(title="Index")

        # Plot the signal data
        x = np.arange(signal_length)
        fig.add_trace(go.Scattergl(
            x=x, y=nth_read_data,
            mode='markers',
            marker=dict(color='LightSkyBlue',
                        size=5,
                        line=dict(color='MediumPurple', width=2)),
            opacity=0.5))

        # Update the plot style (using 0-100 to improve performance)
        fig.update_layout(
            title=nth_read_name,
            yaxis_title="Signal",
            showlegend=False,
            font=dict(size=PLOT_FONT_SIZE),
            xaxis=dict(range=[0, 100])
        )
        fig.update_traces(marker={'size': marker_size})
        # fig.update_xaxes(title="Index")

        # Append the dynamic HTML object to the output structure
        dynamic_html = fig.to_html(full_html=False)
        output_html_plots.update({nth_read_name: dynamic_html})

    # Update the plot filepaths
    plot_filepaths['ont_signal']['dynamic'] = output_html_plots

    return plot_filepaths


# Plot the ONT FAST5 signal data
def plot_signal(output_data, para_dict):
    """Plot the ONT FAST5 signal data for a random sample of reads."""
    
    # Get input parameters
    output_dir = para_dict["output_folder"]
    marker_size = 10
    read_count_max = para_dict["read_count"]
    
    # Get read and base counts
    read_count = output_data.getReadCount()
    if read_count == 0:
        raise ValueError("No reads found in the dataset")

    # Randomly sample a small set of reads if it is a large dataset
    read_sample_size = min(read_count_max, read_count)
    unsampled_indices = list(range(0, read_sample_size))
    read_indices = sample(unsampled_indices, read_sample_size)

    if read_sample_size < read_count:
        logging.info("Randomly sampling {} reads from the total of {} reads".format(read_sample_size, read_count))
    else:
        logging.info("Plotting signal data for all {} reads".format(read_count))

    # Plot the reads
    output_html_plots = {}
    for read_index in read_indices:
        # Create the figure
        fig = go.Figure()

        # Get the read data
        nth_read_name = output_data.getNthReadName(read_index)
        nth_read_data = output_data.getNthReadBaseSignals(read_index)
        logging.info("Signal data count for read {}: {}".format(nth_read_name, len(nth_read_data)))
        nth_read_means = output_data.getNthReadBaseMeans(read_index)
        nth_read_stds = output_data.getNthReadBaseStds(read_index)
        nth_read_medians = output_data.getNthReadBaseMedians(read_index)
        nth_read_skewness = output_data.getNthReadPearsonSkewnessCoeff(read_index)
        nth_read_kurtosis = output_data.getNthReadKurtosis(read_index)
        nth_read_sequence = output_data.getNthReadSequence(read_index)
        sequence_length = len(nth_read_data)

        # Check if sequence data is available
        sequence_available = True if nth_read_sequence else False

        # Set up the output CSVs
        csv_qc_filepath = os.path.join(output_dir, nth_read_name + '_QC.csv')
        qc_file = open(csv_qc_filepath, 'w', encoding='utf-8')
        qc_writer = csv.writer(qc_file)
        qc_writer.writerow(["Base", "Raw_Signal", "Length", "Mean", "Median", "StdDev", "PearsonSkewnessCoeff", "Kurtosis"])

        # Loop through the data
        first_index = 0
        last_index = sequence_length
        start_index = 0
        sequence_list = list(nth_read_sequence)
        base_tick_values = []  # Append the last indices of the base signal to use for tick values
        for i in range(first_index, last_index):
            base_signals = nth_read_data[i]  # Get the base's signal
            signal_length = len(base_signals)
            end_index = start_index + signal_length
            base_tick_values.append(end_index)

            # Plot
            x = np.arange(start_index, end_index, 1)
            fig.add_trace(go.Scattergl(
                x=x, y=base_signals,
                mode='markers',
                marker=dict(color='LightSkyBlue',
                            size=5,
                            line=dict(color='MediumPurple', width=2)),
                opacity=0.5))

            # Update CSVs
            base_value = sequence_list[i] if sequence_available else ''
            signal_mean = nth_read_means[i]
            signal_median = nth_read_medians[i]
            signal_stds = nth_read_stds[i]
            signal_skewness = nth_read_skewness[i]
            signal_kurtosis = nth_read_kurtosis[i]
            raw_row = \
                [base_value, base_signals, signal_length,
                 signal_mean, signal_median, signal_stds,
                 signal_skewness, signal_kurtosis]

            qc_writer.writerow(raw_row)

            # Update the index
            start_index = end_index

        # Close CSVs
        qc_file.close()

        # Update the plot style
        fig.update_layout(
            title=nth_read_name,
            yaxis_title="Signal",
            showlegend=False,
            font=dict(size=PLOT_FONT_SIZE),
            xaxis=dict(range=[0, 100])
        )
        fig.update_traces(marker={'size': marker_size})

        if sequence_available:
            # Set up X tick labels
            x_tick_labels = sequence_list[first_index:last_index]
            fig.update_xaxes(title="Base",
                             tickangle=0,
                             tickmode='array',
                             tickvals=base_tick_values,
                             ticktext=x_tick_labels)
        else:
            fig.update_xaxes(title="Index")

        # Append the dynamic HTML object to the output structure
        dynamic_html = fig.to_html(full_html=False)
        output_html_plots.update({nth_read_name: dynamic_html})

    return output_html_plots

def format_cell(value, type_str='int', error_flag=False):
    """Format the cell value for the summary table."""
    style = "background-color: #F88379;" if error_flag else ""
    if type_str == 'int':
        return "<td style=\"text-align:right;{}\">{:,d}</td>".format(style, value)
    elif type_str == 'float':
        return "<td style=\"text-align:right;{}\">{:.1f}</td>".format(style, value)
    else:
        logging.error("ERROR: Invalid type for formatting cell value")

def format_row(row_name, values, type_str='int', col_ignore=None):
    """Format the row for the summary table. Skip flagging null values in specific columns."""
    cell_str = []
    row_flag = False
    for i, value in enumerate(values):
        # Set the error flag if the value is 0 except for unmapped reads
        error_flag = value == 0 and i != col_ignore
        row_flag = row_flag or error_flag  # Flag for the entire row
        cell_str.append(format_cell(value, type_str, error_flag))

    return "<tr><td>{}</td>{}</tr>".format(row_name, "".join(cell_str)), row_flag


def create_summary_table(output_data, plot_filepaths, file_type):
    """Create the summary table for the basic statistics."""
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Summary Table"

    # Decide the file type label
    file_type_label = file_type
    if file_type == 'FAST5s':
        file_type_label = 'FAST5'
    elif file_type == 'SeqTxt':
        file_type_label = 'Basecall Summary'
        
    plot_filepaths["basic_st"]['description'] = "{} Basic Statistics".format(file_type_label)
    table_error_flag = False

    if file_type == 'BAM':

        # Add alignment statistics to the summary table
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Mapped</th><th>Unmapped</th><th>All</th></tr>\n" \
                    "</thead> "
        table_str += "\n<tbody>"

        # Total reads
        row_str, row_flag = format_row("Total Reads", \
                                        [output_data.mapped_long_read_info.total_num_reads, \
                                            output_data.unmapped_long_read_info.total_num_reads, \
                                            output_data.long_read_info.total_num_reads], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag
        
        # Total bases
        row_str, row_flag = format_row("Total Bases", \
                                        [output_data.mapped_long_read_info.total_num_bases, \
                                         output_data.unmapped_long_read_info.total_num_bases, \
                                         output_data.long_read_info.total_num_bases], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Longest read length
        row_str, row_flag = format_row("Longest Read Length", \
                                        [output_data.mapped_long_read_info.longest_read_length, \
                                         output_data.unmapped_long_read_info.longest_read_length, \
                                         output_data.long_read_info.longest_read_length], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # N50
        row_str, row_flag = format_row("N50", \
                                        [output_data.mapped_long_read_info.n50_read_length, \
                                            output_data.unmapped_long_read_info.n50_read_length, \
                                            output_data.long_read_info.n50_read_length], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # GC content
        row_str, row_flag = format_row("GC Content(%)", \
                                        [output_data.mapped_long_read_info.gc_cnt * 100, \
                                            output_data.unmapped_long_read_info.gc_cnt * 100, \
                                            output_data.long_read_info.gc_cnt * 100], \
                                        'float', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Mean read length
        row_str, row_flag = format_row("Mean Read Length", \
                                        [output_data.mapped_long_read_info.mean_read_length, \
                                            output_data.unmapped_long_read_info.mean_read_length, \
                                            output_data.long_read_info.mean_read_length], \
                                        'float', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Median read length
        row_str, row_flag = format_row("Median Read Length", \
                                        [output_data.mapped_long_read_info.median_read_length, \
                                            output_data.unmapped_long_read_info.median_read_length, \
                                            output_data.long_read_info.median_read_length], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag
        
    elif file_type == 'SeqTxt':
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Passed</th><th>Failed</th><th>All</th></tr>\n</thead>"
        table_str += "\n<tbody>"
        
        # Total reads
        row_str, row_flag = format_row("Total Reads", \
                                        [output_data.passed_long_read_info.long_read_info.total_num_reads, \
                                            output_data.failed_long_read_info.long_read_info.total_num_reads, \
                                            output_data.all_long_read_info.long_read_info.total_num_reads], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Total bases
        row_str, row_flag = format_row("Total Bases", \
                                        [output_data.passed_long_read_info.long_read_info.total_num_bases, \
                                            output_data.failed_long_read_info.long_read_info.total_num_bases, \
                                            output_data.all_long_read_info.long_read_info.total_num_bases], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Longest read length
        row_str, row_flag = format_row("Longest Read Length", \
                                        [output_data.passed_long_read_info.long_read_info.longest_read_length, \
                                            output_data.failed_long_read_info.long_read_info.longest_read_length, \
                                            output_data.all_long_read_info.long_read_info.longest_read_length], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # N50
        row_str, row_flag = format_row("N50", \
                                        [output_data.passed_long_read_info.long_read_info.n50_read_length, \
                                            output_data.failed_long_read_info.long_read_info.n50_read_length, \
                                            output_data.all_long_read_info.long_read_info.n50_read_length], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Mean read length
        row_str, row_flag = format_row("Mean Read Length", \
                                        [output_data.passed_long_read_info.long_read_info.mean_read_length, \
                                            output_data.failed_long_read_info.long_read_info.mean_read_length, \
                                            output_data.all_long_read_info.long_read_info.mean_read_length], \
                                        'float', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Median read length
        row_str, row_flag = format_row("Median Read Length", \
                                        [output_data.passed_long_read_info.long_read_info.median_read_length, \
                                            output_data.failed_long_read_info.long_read_info.median_read_length, \
                                            output_data.all_long_read_info.long_read_info.median_read_length], \
                                        'int', 1)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

    elif file_type == 'FAST5s':
        # Get values
        read_count = output_data.getReadCount()
        total_base_count = output_data.getTotalBaseCount()

        # Set up the HTML table
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
        table_str += "\n<tbody>"

        # Total reads
        row_str, row_flag = format_row("Total Reads", [read_count], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Total bases
        row_str, row_flag = format_row("Total Bases", [total_base_count], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

    else:
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
        table_str += "\n<tbody>"
        # Total reads
        row_str, row_flag = format_row("Total Reads", [output_data.long_read_info.total_num_reads], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Total bases
        row_str, row_flag = format_row("Total Bases", [output_data.long_read_info.total_num_bases], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Longest read length
        row_str, row_flag = format_row("Longest Read Length", [output_data.long_read_info.longest_read_length], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # N50
        row_str, row_flag = format_row("N50", [output_data.long_read_info.n50_read_length], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # GC content
        row_str, row_flag = format_row("GC Content(%)", [output_data.long_read_info.gc_cnt * 100], 'float', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Mean read length
        row_str, row_flag = format_row("Mean Read Length", [output_data.long_read_info.mean_read_length], 'float', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        # Median read length
        row_str, row_flag = format_row("Median Read Length", [output_data.long_read_info.median_read_length], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag
        
    table_str += "\n</tbody>\n</table>"
    # table_str += """
    #     <div class="help-icon">
    #         💡
    #         <div class="tooltip">This is your help text explaining the feature!</div>
    #     </div>
    #     """
    plot_filepaths["basic_st"]['detail'] = table_str
    plot_filepaths["basic_st"]['error_flag'] = table_error_flag


def get_axis_name(row, axis_type='x'):
    """Get the axis name for the plot."""
    axis_number = row + 1
    return f"{axis_type}axis{axis_number}" if axis_number > 1 else f"{axis_type}axis"


def create_modified_base_table(output_data, plot_filepaths, base_modification_threshold):
    """Create a summary table for the base modifications."""
    help_text = "Total unfiltered predictions are all predictions prior to applying the base modification probability threshold.\n" \
                "This threshold is set by the user (default: 0.5) and is used to filter out low-confidence base modifications.\n" \
                "Total modification counts are the number of base modifications that pass the threshold.\n" \
                "These counts are also separated by forward and reverse strand predictions.\n" \
                "CpG modification counts are the total CpG modifications that pass the threshold.\n" \
                "These are total counts and not site-specific counts." \
                
    plot_filepaths["base_mods"] = {}
    plot_filepaths["base_mods"]['file'] = ""
    plot_filepaths["base_mods"]['title'] = "Base Modifications"
    plot_filepaths["base_mods"]['description'] = "Base modification statistics"
    table_error_flag = False

    # Print the types of modifications
    logging.info("Getting base modification types")
    base_mod_types = output_data.getBaseModTypes()
    logging.info("[TEST] Modification types: ")
    if base_mod_types:
        logging.info("Modification types: ")
        for mod_type in base_mod_types:
            logging.info(mod_type)

        logging.info("Getting base modification statistics")

        # Get the read length vs. base modification rate data for each
        # # modification type
        # logging.info("Getting mod data size")
        # read_mod_data_size = output_data.getReadModDataSize()
        # logging.info("Mod data size: {}".format(read_mod_data_size))

        # # Choose a maximum of 10,000 reads to randomly sample for the plot
        # max_reads = min(read_mod_data_size, 10000)        
        # # read_indices = set(sample(range(read_mod_data_size), max_reads))
        # read_indices = np.random.choice(read_mod_data_size, max_reads, replace=False)
        # read_length_mod_rates = {}

        # Get the read length (%) vs. base modification probability data for
        # each sampled read
        sample_count = 10000
        read_len_pct = []
        mod_prob = []
        for mod_type in base_mod_types:
            for i in range(sample_count):
                try:
                    pct = output_data.getNthReadLenPct(i, mod_type)
                    prob = output_data.getNthReadModProb(i, mod_type)
                    read_len_pct.append(pct)
                    mod_prob.append(prob)
                except Exception as e:
                    logging.error(f"Error getting read length vs. base modification probability data: {e}")

        # Convert the lists to numpy arrays
        read_len_pct = np.array(read_len_pct) * 100  # Convert to percentage
        mod_prob = np.array(mod_prob)
        
        # Dictionary of modification character to full name
        mod_char_to_name = {'m': '5mC', 'h': '5hmC', 'f': '5fC', 'c': '5caC', \
                            'g': '5hmU', 'e': '5fu', 'b': '5caU', \
                            'a': '6mA', 'o': '8oxoG', 'n': 'Xao', \
                            'C': 'Amb. C', 'A': 'Amb. A', 'T': 'Amb. T', 'G': 'Amb. G',\
                            'N': 'Amb. N', \
                            'v': 'pseU'}

        # Create a plot of pct read length vs. base modification probability for
        # each modification type, as well as a histogram of the average base
        # modification probability for 100 bins of the read length

        # Make a subplot of two columns for the read length vs. base
        # modification probability and the histogram of the average base
        # modification probability for each modification type
        fig = make_subplots(rows=len(base_mod_types), cols=2, shared_xaxes=False, shared_yaxes=False, vertical_spacing=0.1, subplot_titles=[f"{mod_char_to_name[mod_type]} Modification Probability" for mod_type in base_mod_types])

        for i, mod_type in enumerate(base_mod_types):
            logging.info(f"Creating trace for modification type: {mod_type} at row: {i + 1}")

            # Add the trace for the read length vs. base modification
            # probability scatter plot
            fig.add_trace(go.Scatter
                (x=read_len_pct, y=mod_prob, mode='markers', name=mod_char_to_name[mod_type], marker=dict(size=5), showlegend=False),
                row=i + 1, col=1)
            
            # Print the first 50 pairs sorted by read length for debugging
            # read_len_pct, mod_prob = zip(*sorted(zip(read_len_pct, mod_prob)))
            # if i == 0:
            #     for j in range(50):
            #         logging.info(f"Read length: {read_len_pct[j]}, Modification probability: {mod_prob[j]}")
            
            # Create a histogram of the base modification probabilities
            base_mod_prob_hist = go.Histogram(x=mod_prob, name=mod_char_to_name[mod_type], showlegend=False, nbinsx=20)
            fig.add_trace(base_mod_prob_hist, row=i + 1, col=2)
            
            # Add a bar plot of the average base modification probability for
            # 100 bins of the read length
            # bins = np.linspace(0, 100, 11)  # 10 bins (0-10%, 10-20%, ..., 90-100%)
            # bin_centers = (bins[:-1] + bins[1:]) / 2  # Bin centers for plotting

            # # Get the average probability per bin
            # avg_prob_per_bin = np.zeros(10)
            # bin_indices = np.digitize(read_len_pct, bins) - 1
            # for j in range(10):  # Loop over bins
            #     bin_mask = (bin_indices == j)
            #     if np.any(bin_mask):
            #         avg_prob_per_bin[j] = np.mean(mod_prob[bin_mask])
            #         logging.info(f"Bin {j}: {avg_prob_per_bin[j]}")

            # # Create the bar plot

            # # Print the bins and read length percentages for the first 10 reads
            # # for debugging
            # if i == 0:
            #     logging.info("Bins: {}".format(bins))
            #     logging.info("Bin indices: {}".format(bin_indices[:10]))
            #     logging.info("Read length percentages: {}".format(read_len_pct[:10]))

            # # Create the bar plot
            # fig.add_trace(go.Bar(x=bin_centers, y=avg_prob_per_bin, name=mod_char_to_name[mod_type], showlegend=False), row=i + 1, col=2)

            # Update the plot style
            fig.update_xaxes(title="Read Length (%)", row=i + 1, col=1)
            fig.update_yaxes(title="Modification Probability", row=i + 1, col=1)
            fig.update_xaxes(title="Modification Probability", row=i + 1, col=2)
            fig.update_yaxes(title="Frequency", row=i + 1, col=2)
            # fig.update_xaxes(title="Read Length (%)", row=i + 1, col=2)
            # fig.update_yaxes(title="Average Modification Probability", row=i + 1, col=2)

            # Set the range of the y-axis to 0-1
            fig.update_yaxes(range=[0, 1], row=i + 1, col=1)
            # fig.update_yaxes(range=[0, 1], row=i + 1, col=2)

        # Update the plot layout
        fig.update_layout(title="Read Length vs. Base Modification Probability", font=dict(size=PLOT_FONT_SIZE))
            
        # Generate the HTML
        if len(base_mod_types) > 0:
            plot_height = 500 * len(base_mod_types)
            plot_width = 700 * 2
            logging.info("Saving the read length vs. modification rates plot")
            plot_filepaths["read_length_mod_rates"]['dynamic'] = fig.to_html(full_html=False, default_height=plot_height, default_width=plot_width)
    else:
        logging.warning("WARNING: No modification types found")

    # Create the base modification statistics table'
    logging.info("Creating the base modification statistics table")
    table_str = "<table>\n<tbody>"
    row_str, row_flag = format_row("Total Unfiltered Predictions", [output_data.modified_prediction_count], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    row_str, row_flag = format_row("Probability Threshold", [base_modification_threshold], 'float', 0)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    row_str, row_flag = format_row("Total Modification Counts", [output_data.sample_modified_base_count], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    row_str, row_flag = format_row("Total Modification Counts (Forward Strand Only)", [output_data.sample_modified_base_count_forward], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    row_str, row_flag = format_row("Total Modification Counts (Reverse Strand Only)", [output_data.sample_modified_base_count_reverse], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    row_str, row_flag = format_row("Total CpG Modification Counts (Forward Strand Only)", [output_data.sample_cpg_forward_count], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    row_str, row_flag = format_row("Total CpG Modification Counts (Reverse Strand Only)", [output_data.sample_cpg_reverse_count], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    # Add the modification type data
    for mod_type in base_mod_types:
        # mod_name = mod_char_to_name[mod_type]
        try:
            mod_name = mod_char_to_name[mod_type]
        except KeyError:
            logging.warning("WARNING: Unknown modification type: {}".format(mod_type))
            mod_name = mod_type

        mod_count = output_data.getModTypeCount(mod_type)
        mod_count_fwd = output_data.getModTypeCount(mod_type, 0)
        mod_count_rev = output_data.getModTypeCount(mod_type, 1)

        row_str, row_flag = format_row("Total {} Counts in the Sample".format(mod_name), [mod_count], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        row_str, row_flag = format_row("Total {} Counts in the Sample (Forward Strand)".format(mod_name), [mod_count_fwd], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

        row_str, row_flag = format_row("Total {} Counts in the Sample (Reverse Strand)".format(mod_name), [mod_count_rev], 'int', None)
        table_str += row_str
        table_error_flag = table_error_flag or row_flag

    # Finish the table
    table_str += "\n</tbody>\n</table>"

    # Add the help text
    table_str += """
        <div class="help-icon">
            💡
            <div class="tooltip">{}</div>
        </div>
        """.format(help_text)
    
    # Add text below the table suggesting the user to use Modkit for more
    # detailed analysis on per-site modification rates
    table_str += "<p><i>For per-site modification rates, please use \
        <a href=\"https://github.com/nanoporetech/modkit\">Modkit</a> by Oxford Nanopore Technologies..</i></p>"


    plot_filepaths["base_mods"]['detail'] = table_str
    plot_filepaths["base_mods"]['error_flag'] = table_error_flag

def create_tin_table(output_data, input_files, plot_filepaths):
    """Create a summary table for the RNA-Seq TIN values."""
    plot_filepaths["tin"] = {}
    plot_filepaths["tin"]['file'] = ""
    plot_filepaths["tin"]['title'] = "Transcript Integrity Number (TIN)"
    plot_filepaths["tin"]['description'] = "RNA-Seq TIN values"

    # Create a table with the first column showing the BAM filepath, and the
    # following columns showing TIN count, mean, median, and standard deviation
    table_str = "<table>\n<thead>\n<tr><th>BAM File</th><th>Median TIN Score</th><th>Number of Transcripts</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    
    # Loop through each BAM file
    error_flag = False
    for bam_file in input_files:
        # Format the filepath as filename only
        bam_filename = os.path.basename(bam_file)

        # Get the file data
        # tin_count = output_data.getTINCount(bam_file)
        # tin_mean = output_data.getTINMean(bam_file)
        tin_median = output_data.getTINMedian(bam_file)
        # tin_std = output_data.getTINStdDev(bam_file)

        # Add the data to the table
        # row_str, row_flag = format_row(bam_filename, [tin_count, tin_mean,
        # tin_median, tin_std], 'float', None)
        row_str, row_flag = format_row(bam_filename, [tin_median, output_data.getTINCount(bam_file)], 'float', None)
        table_str += row_str
        error_flag = error_flag or row_flag

    table_str += "\n</tbody>\n</table>"

    # Add the table to the plot filepaths
    plot_filepaths["tin"]['detail'] = table_str
    plot_filepaths["tin"]['error_flag'] = error_flag


def create_pod5_table(output_dict, plot_filepaths):
    """Create a summary table for the ONT POD5 signal data."""
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Summary Table"
    file_type_label = "POD5"
    plot_filepaths["basic_st"]['description'] = f"{file_type_label} Basic Statistics"
    table_error_flag = False
    
    # Set up the HTML table
    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    # int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    # table_str += int_str_for_format.format("Total Reads", read_count)
    read_count = len(output_dict.keys())
    row_str, row_flag = format_row("Total Reads", [read_count], 'int', None)
    table_str += row_str
    table_error_flag = table_error_flag or row_flag

    table_str += "\n</tbody>\n</table>"
    plot_filepaths["basic_st"]['detail'] = table_str
    plot_filepaths["basic_st"]['error_flag'] = table_error_flag


def plot_alignment_numbers(data, plot_filepaths):
    category = ['Primary Alignments', 'Supplementary Alignments', 'Secondary Alignments',
                'Reads with Supplementary Alignments', 'Reads with Secondary Alignments',
                'Reads with Secondary and Supplementary Alignments', 'Forward Alignments', 'Reverse Alignments']
    category = [wrap(x) for x in category]

    # Set the error flag if primary alignments equal 0
    error_flag = data.num_primary_alignment == 0

    logging.info("[TEST] Number of reverse alignments: {}".format(data.reverse_alignment))
    logging.info("[TEST] Number of forward alignments: {}".format(data.forward_alignment))

    # Create a horizontally aligned bar plot trace from the data using plotly
    trace = go.Bar(x=[data.num_primary_alignment, data.num_supplementary_alignment, data.num_secondary_alignment,
                      data.num_reads_with_supplementary_alignment, data.num_reads_with_secondary_alignment,
                      data.num_reads_with_both_secondary_supplementary_alignment, data.forward_alignment,
                      data.reverse_alignment], y=category, orientation='h')

    # Create the layout for the plot
    layout = go.Layout(title=go.layout.Title(text=""),
                       xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Counts")),
                       yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")),
                       font=dict(size=PLOT_FONT_SIZE))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Update the HTML data for the plot
    plot_filepaths['read_alignments_bar']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=1000)
    plot_filepaths['read_alignments_bar']['error_flag'] = error_flag


def plot_errors(output_data, plot_filepaths):
    """Plot the error statistics for the alignment data."""
    category = \
        ['Matched Bases', 'Mismatched Bases', 'Inserted Bases', 'Deleted Bases', 'Clipped Bases\n(Primary Alignments)']
    category = [wrap(x) for x in category]

    # Create a horizontally aligned bar plot trace from the data using plotly
    trace = go.Bar(x=[output_data.num_matched_bases, output_data.num_mismatched_bases, output_data.num_ins_bases,
                      output_data.num_del_bases, output_data.num_clip_bases], y=category, orientation='h')

    # Create the layout for the plot
    layout = go.Layout(title=go.layout.Title(text=""),
                       xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Counts")),
                       yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")),
                       font=dict(size=PLOT_FONT_SIZE))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Generate the HTML object for the plot
    # html_obj = fig.to_html(full_html=False, default_height=500,
    # default_width=700)
    plot_filepaths['base_alignments_bar']['dynamic'] = fig.to_html(full_html=False, default_height=500, default_width=700)

    # Set the error flag if mismatch or clipped bases > matched bases
    error_flag = output_data.num_mismatched_bases > output_data.num_matched_bases or \
                 output_data.num_clip_bases > output_data.num_matched_bases
    plot_filepaths['base_alignments_bar']['error_flag'] = error_flag

    # return html_obj

