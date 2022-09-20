"""
plot_for_FAST5s.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

from src import lrst_global
# from src.utils import *

import os
import numpy as np
import plotly.graph_objs as go


def plot(fast5_output, para_dict):
    """
    Update the global variables with HTML strings using the output data.
    """
    out_path = para_dict["output_folder"]
    get_image_path = lambda x: os.path.join(out_path, lrst_global.plot_filenames[x]['file'])

    # Set up the global variable with HTML titles
    lrst_global.plot_filenames["basic_st"] = {}
    lrst_global.plot_filenames["basic_st"]['file'] = ""
    lrst_global.plot_filenames["basic_st"]['title'] = "Basic statistics"
    lrst_global.plot_filenames["basic_st"]['description'] = "FAST5: Basic statistics"

    # Get values
    read_count = fast5_output.getReadCount()
    total_base_count = fast5_output.getTotalBaseCount()

    # Set up the HTML table
    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    table_str += int_str_for_format.format("#Total Reads", read_count)
    table_str += int_str_for_format.format("#Total Bases", total_base_count)
    table_str += "\n</tbody>\n</table>"
    lrst_global.plot_filenames["basic_st"]['detail'] = table_str

    # plot_read_length_stats([fast5_output.long_read_info], get_image_path('read_length_st'), subtitles=['Long Reads'])
    # plot_base_counts([fast5_output.long_read_info], get_image_path('base_st'), subtitles=['Long Reads'])
    # plot_basic_info([fast5_output.long_read_info], get_image_path('basic_info'), categories=['Long Reads'])

    # Plot the reads
    for read_index in range(read_count):
        # Create the figure
        fig = go.Figure()

        # Get the read data
        nth_read_name = fast5_output.getNthReadName(read_index)
        nth_read_data = fast5_output.getNthReadBaseSignals(read_index)
        nth_read_means = fast5_output.getNthReadBaseMeans(read_index)
        nth_read_stds = fast5_output.getNthReadBaseStds(read_index)
        nth_read_medians = fast5_output.getNthReadBaseMedians(read_index)
        nth_read_sequence = fast5_output.getNthReadSequence(read_index)

        sequence_length = len(nth_read_sequence)
        start_index = 0
        for i in range(sequence_length):
            base_signals = nth_read_data[i]
            end_index = start_index + len(base_signals)
            x = np.arange(start_index, end_index, 1)
            fig.add_trace(go.Scatter(x=x, y=base_signals, mode='markers', opacity=0.5, fillcolor='green'))

        start_index = end_index

        test_filepath = os.path.join(out_path, 'fig1.png')
        print("saving plot to ", test_filepath, "...")
        fig.write_image(test_filepath)
        print("Plot generated.")
        0

    # lrst_global.plot_filenames['read_length_hist']['dynamic'] = histogram(fast5_output.long_read_info,
    #                                                                       get_image_path('read_length_hist'))
    # lrst_global.plot_filenames['base_quality']['dynamic'] = base_quality(fast5_output.seq_quality_info,
    #                                                                      get_image_path('base_quality'))
    # lrst_global.plot_filenames['read_avg_base_quality']['dynamic'] = read_avg_base_quality(fast5_output.seq_quality_info,
    #                                                                                        get_image_path(
    #                                                                                            'read_avg_base_quality'))
