import os
import numpy as np
import csv
from random import sample

import plotly.graph_objs as go
from plotly.subplots import make_subplots

# Set up logging
import logging
logging.basicConfig(level=logging.INFO)

# Constants
MAX_BASE_QUALITY = 100
MAX_READ_QUALITY = 100
PLOT_FONT_SIZE = 16

# Return a dictionary of default plot filenames
def getDefaultPlotFilenames():
    plot_filenames = {  # for fq/fa
        "read_length_distr": {'title': "Read Length", 'description': "Read Length Distribution"},  # for bam
        "read_alignments_bar": {'title': "Read Alignments",
                   'description': "Read Alignments"},
        "base_alignments_bar": {'title': "Base Alignment and Error",
                   'description': "Base Alignment and Error"},
        "read_length_bar": {'title': "Read Length Statistics", 'description': "Read Length Statistics"},
        "base_counts": {'title': "Base Counts",
                    'description': "Base Counts", 'summary': ""},
        "basic_info": {'title': "Basic Statistics",
                       'description': "Basic Statistics", 'summary': ""},
        "read_length_hist": {'title': "Read Length Histogram", 'description': "Read Length Histogram", 'summary': ""},

        "base_quality": {'title': "Base Quality Histogram", 'description': "Base Quality Histogram"},

        "read_avg_base_quality": {'title': "Read Base Quality Histogram", 'description': "Read Base Quality Histogram"},

        "pos_quality": {'title': "Base Position Quality", 'description': "Base Position Quality"},
        "ont_signal": {'title': "ONT Signal", 'description': "ONT Signal"},
    }

    return plot_filenames

# Wrap the text in the table 
def wrap(label):
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

# Plot the read alignment numbers
def plot_read_length_stats(output_data, file_type):

    # Define the three categories
    category = ['N50', 'Mean', 'Median']
    all_traces = []

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

    else:
        # Get the data for all reads
        key_list = ['n50_read_length', 'mean_read_length', 'median_read_length']

        # Create a bar trace
        bar_title = 'All Reads'
        data = output_data.long_read_info
        values = [getattr(data, key_name) for key_name in key_list]
        trace = go.Bar(x=category, y=values, name=bar_title)
        all_traces.append(trace)

    # Create the layout
    layout = go.Layout(title='', xaxis=dict(title='Statistics'), yaxis=dict(title='Length (bp)'), barmode='group', font=dict(size=PLOT_FONT_SIZE))

    # Create the figure and add the traces
    fig = go.Figure(data=all_traces, layout=layout)

    # Generate the HTML
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj

# Plot the base counts
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

    elif filetype == 'SeqTxt':
        bar_titles = ['All Reads', 'Passed Reads', 'Failed Reads']
        data_objects = [output_data.all_long_read_info.long_read_info, output_data.passed_long_read_info.long_read_info, output_data.failed_long_read_info.long_read_info]
        for i in range(3):
            plot_title = bar_titles[i]
            data = data_objects[i]
            values = [data.total_a_cnt, data.total_c_cnt, data.total_g_cnt, data.total_tu_cnt, data.total_n_cnt]
            trace = go.Bar(x=category, y=values, name=plot_title)
            all_traces.append(trace)

    else:
        plot_title = 'All Reads'
        data = output_data.long_read_info
        values = [data.total_a_cnt, data.total_c_cnt, data.total_g_cnt, data.total_tu_cnt, data.total_n_cnt]
        trace = go.Bar(x=category, y=values, name=plot_title)
        all_traces.append(trace)

    # Create the layout
    layout = go.Layout(title='', xaxis=dict(title='Base'), yaxis=dict(title='Counts'), barmode='group', font=dict(size=PLOT_FONT_SIZE))

    # Create the figure and add the traces
    fig = go.Figure(data=all_traces, layout=layout)

    # Generate the HTML
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj

# Plot basic information about the reads in bar chart format
def plot_basic_info(output_data, file_type):
    html_obj = ''
    if file_type == 'BAM':

        # Create a bar trace for each type of data
        bar_titles = ['All Reads', 'Mapped Reads', 'Unmapped Reads']
        data_objects = [output_data.long_read_info, output_data.mapped_long_read_info, output_data.unmapped_long_read_info]

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
            fig.update_layout(showlegend=False)

        # Update the layout
        fig.update_layout(showlegend=False, font=dict(size=PLOT_FONT_SIZE))

        # Generate the HTML
        html_obj = fig.to_html(full_html=False, default_height=800, default_width=1200)

    elif file_type == 'SeqTxt':

        # Create a bar trace for each type of data
        bar_titles = ['All Reads', 'Passed Reads', 'Failed Reads']
        data_objects = [output_data.all_long_read_info.long_read_info, output_data.passed_long_read_info.long_read_info, output_data.failed_long_read_info.long_read_info]

        # Create subplots for each category
        fig = make_subplots(rows=1, cols=3, subplot_titles=("Number of Reads", "Number of Bases", "Longest Read"), horizontal_spacing=0.1)

        # Add traces for each category
        key_list = ['total_num_reads', 'total_num_bases', 'longest_read_length']
        for i in range(3):
            # Get the data for this category
            key_name = key_list[i]

            # Add the traces for each type of data
            data = [getattr(data_objects[0], key_name), getattr(data_objects[1], key_name), getattr(data_objects[2], key_name)]

            # Create the trace
            trace = go.Bar(x=data, y=bar_titles, orientation='h')

            # Add the trace to the figure
            fig.add_trace(trace, row=1, col=i + 1)

        # Update the layout
        fig.update_layout(showlegend=False, font=dict(size=PLOT_FONT_SIZE))

        # Generate the HTML
        html_obj = fig.to_html(full_html=False, default_height=500, default_width=1600)

    return html_obj


# Plot the read length histograms
def read_lengths_histogram(data, font_size):
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

    # Update the layout
    fig.update_layout(showlegend=False, autosize=True, font=dict(size=PLOT_FONT_SIZE))
    # Set font sizes
    # fig.update_layout(showlegend=False, autosize=False)
    # fig.update_layout(font=dict(size=font_size), autosize=True)

    fig.update_annotations(font_size=annotation_size)
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj

# Save the 'Base quality' plot image.
def base_quality(data, font_size):
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

    return fig.to_html(full_html=False, default_height=500, default_width=700)

# Save the 'Average base quality' plot image.
def read_avg_base_quality(data, font_size):
    xd = np.arange(MAX_READ_QUALITY)
    yd = np.array(data.read_average_base_quality_distribution)
    fig = go.Figure()
    fig.add_trace(go.Bar(x=xd, y=yd, marker_color='#36a5c7'))

    fig.update_xaxes(ticks="outside", dtick=10, title_text='Average Base Quality', title_standoff=0)
    fig.update_yaxes(ticks="outside", title_text='Number of Reads', title_standoff=0)
    fig.update_layout(font=dict(size=PLOT_FONT_SIZE))  # Set font size

    return fig.to_html(full_html=False, default_height=500, default_width=700)

# Main plot function
def plot(output_data, para_dict, file_type):
    out_path = para_dict["output_folder"]
    plot_filepaths = getDefaultPlotFilenames()

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create the summary table
    create_summary_table(output_data, plot_filepaths, file_type)

    # Generate plots
    plot_filepaths['base_counts']['dynamic'] = plot_base_counts(output_data, file_type)
    plot_filepaths['basic_info']['dynamic'] = plot_basic_info(output_data, file_type)

    # Read length histogram
    if file_type == 'SeqTxt':
        long_read_data = output_data.all_long_read_info.long_read_info
    else:
        long_read_data = output_data.long_read_info

    if file_type != 'FAST5s':
        plot_filepaths['read_length_hist']['dynamic'] = read_lengths_histogram(long_read_data, font_size)

        plot_filepaths['read_length_bar']['dynamic'] = plot_read_length_stats(output_data, file_type)

    if file_type != 'FASTA' and file_type != 'FAST5s':
        if file_type == 'SeqTxt':
            seq_quality_info = output_data.all_long_read_info.seq_quality_info
        else:
            seq_quality_info = output_data.seq_quality_info

        # Base quality histogram
        plot_filepaths['base_quality']['dynamic'] = base_quality(seq_quality_info, font_size)

        # Read quality histogram
        read_quality_dynamic = read_avg_base_quality(seq_quality_info, font_size)
        plot_filepaths['read_avg_base_quality']['dynamic'] = read_quality_dynamic

    if file_type == 'BAM':
        plot_filepaths['read_alignments_bar']['dynamic'] = plot_alignment_numbers(output_data)
        plot_filepaths['base_alignments_bar']['dynamic'] = plot_errors(output_data)
        
    elif file_type == 'FAST5s':
        plot_filepaths['ont_signal']['dynamic'] = plot_signal(output_data, para_dict)

    return plot_filepaths

# Plot the ONT FAST5 signal data
def plot_signal(output_data, para_dict):
    
    # Get input parameters
    output_dir = para_dict["output_folder"]
    font_size = para_dict["fontsize"]
    marker_size = para_dict["markersize"]
    read_count_max = para_dict["read_count"]
    
    # Get read and base counts
    read_count = output_data.getReadCount()
    logging.info("Plotting signal data for {} reads".format(read_count))

    # Randomly sample a small set of reads if it is a large dataset
    read_sample_size = min(read_count_max, read_count)
    unsampled_indices = list(range(0, read_sample_size))
    read_indices = sample(unsampled_indices, read_sample_size)

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
        qc_file = open(csv_qc_filepath, 'w')
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
            fig.add_trace(go.Scatter(
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
            font=dict(size=PLOT_FONT_SIZE)
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

# Create a summary table for the basic statistics
def create_summary_table(output_data, plot_filepaths, file_type):
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Summary Table"

    # Decide the file type label
    file_type_label = file_type
    if file_type == 'FAST5s':
        file_type_label = 'FAST5'

    plot_filepaths["basic_st"]['description'] = "{} Basic statistics".format(file_type_label)

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
        
    elif file_type == 'SeqTxt':
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Passed</th><th>Failed</th><th>All</th></tr>\n</thead>"
        table_str += "\n<tbody>"
        int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td></tr>"
        double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
        table_str += int_str_for_format.format("#Total Reads",
                                               output_data.passed_long_read_info.long_read_info.total_num_reads,
                                               output_data.failed_long_read_info.long_read_info.total_num_reads,
                                               output_data.all_long_read_info.long_read_info.total_num_reads)
        table_str += int_str_for_format.format("#Total Bases",
                                               output_data.passed_long_read_info.long_read_info.total_num_bases,
                                               output_data.failed_long_read_info.long_read_info.total_num_bases,
                                               output_data.all_long_read_info.long_read_info.total_num_bases)
        table_str += int_str_for_format.format("Longest Read Length",
                                               output_data.passed_long_read_info.long_read_info.longest_read_length,
                                               output_data.failed_long_read_info.long_read_info.longest_read_length,
                                               output_data.all_long_read_info.long_read_info.longest_read_length)
        table_str += int_str_for_format.format("N50",
                                               output_data.passed_long_read_info.long_read_info.n50_read_length,
                                               output_data.failed_long_read_info.long_read_info.n50_read_length,
                                               output_data.all_long_read_info.long_read_info.n50_read_length)
        table_str += double_str_for_format.format("Mean Read Length",
                                                  output_data.passed_long_read_info.long_read_info.mean_read_length,
                                                  output_data.failed_long_read_info.long_read_info.mean_read_length,
                                                  output_data.all_long_read_info.long_read_info.mean_read_length)
        table_str += int_str_for_format.format("Median Read Length",
                                               output_data.passed_long_read_info.long_read_info.median_read_length,
                                               output_data.failed_long_read_info.long_read_info.median_read_length,
                                               output_data.all_long_read_info.long_read_info.median_read_length)

    elif file_type == 'FAST5s':
        # Get values
        read_count = output_data.getReadCount()
        total_base_count = output_data.getTotalBaseCount()

        # Set up the HTML table
        table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
        table_str += "\n<tbody>"
        int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
        table_str += int_str_for_format.format("#Total Reads", read_count)
        table_str += int_str_for_format.format("#Total Bases", total_base_count)

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


def plot_alignment_numbers(data):
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
                       yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")),
                       font=dict(size=PLOT_FONT_SIZE))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Generate the HTML object for the plot
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=1000)

    return html_obj


# Plot base alignment statistics
def plot_errors(output_data):
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
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj

