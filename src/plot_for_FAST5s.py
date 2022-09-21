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

    # Plot the reads
    output_html_plots = {}
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
        end_index = 0
        window_size = len(nth_read_data[0])
        for i in range(sequence_length):
            base_signals = nth_read_data[i]
            end_index = start_index + window_size
            x = np.arange(start_index, end_index, 1)

            # df = pd.DataFrame({'Base': x, 'Signal': base_signals})
            fig.add_trace(go.Scatter(
                x=x, y=base_signals,
                mode='markers',
                marker=dict(color='LightSkyBlue',
                            size=5,
                            line=dict(color='MediumPurple', width=2)),
                opacity=0.5))
            start_index = end_index

        # Add labels
        fig.update_layout(
            title=nth_read_name,
            xaxis_title="Base",
            yaxis_title="Signal",
            showlegend=False
        )
        tick_sequence = list(nth_read_sequence)
        fig.update_xaxes(tickangle=45,
                         tickmode='array',
                         tickvals=np.arange(0, end_index, window_size),
                         ticktext=tick_sequence)

        # Save
        image_filepath = os.path.join(out_path, nth_read_name + '_BaseSignal.png')
        print("saving plot to ", image_filepath, "...")
        fig.write_image(image_filepath)

        # Append the dynamic HTML object to the output structure
        dynamic_html = fig.to_html(full_html=False)
        #output_html_plots.append(dynamic_html)
        output_html_plots.update({nth_read_name: dynamic_html})
        print("Plot generated.")

    return output_html_plots
