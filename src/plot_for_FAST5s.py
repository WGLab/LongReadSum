"""
plot_for_FAST5s.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

import os
import logging
import csv
import numpy as np
import plotly.graph_objs as go
from random import sample
\
if __package__ == 'src':
    from src import lrst_global
else:
    import lrst_global

def plot(fast5_output, para_dict):
    """
    Update the global variables with HTML strings using the output data.
    """
    out_path = para_dict["output_folder"]

    # Set up the global variable with HTML titles
    lrst_global.plot_filenames["basic_st"] = {}
    lrst_global.plot_filenames["basic_st"]['file'] = ""
    lrst_global.plot_filenames["basic_st"]['title'] = "Summary Table"
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

    # Randomly sample a small set of reads if it is a large dataset
    read_count_max = para_dict["read_count"]
    read_sample_size = min(read_count_max, read_count)
    unsampled_indices = list(range(0, read_sample_size))
    read_indices = sample(unsampled_indices, read_sample_size)

    # Plot the reads
    output_html_plots = {}
    for read_index in read_indices:
        # Create the figure
        fig = go.Figure()

        # Get the read data
        nth_read_name = fast5_output.getNthReadName(read_index)
        nth_read_data = fast5_output.getNthReadBaseSignals(read_index)
        nth_read_means = fast5_output.getNthReadBaseMeans(read_index)
        nth_read_stds = fast5_output.getNthReadBaseStds(read_index)
        nth_read_medians = fast5_output.getNthReadBaseMedians(read_index)
        nth_read_skewness = fast5_output.getNthReadPearsonSkewnessCoeff(read_index)
        nth_read_kurtosis = fast5_output.getNthReadKurtosis(read_index)
        nth_read_sequence = fast5_output.getNthReadSequence(read_index)
        sequence_length = len(nth_read_data)

        # Check if sequence data is available
        sequence_available = True if nth_read_sequence else False

        # Set up the output CSVs
        csv_qc_filepath = os.path.join(out_path, nth_read_name + '_QC.csv')
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
        font_size = para_dict["fontsize"]
        marker_size = para_dict["markersize"]
        fig.update_layout(
            title=nth_read_name,
            yaxis_title="Signal",
            showlegend=False,
            font=dict(size=font_size)
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

        # Save image
        image_filepath = os.path.join(out_path, "img", nth_read_name + '_BaseSignal.png')
        fig.write_image(image_filepath)
        save_msg = "Plot image saved to: " + image_filepath
        logging.info(save_msg)

        # Append the dynamic HTML object to the output structure
        dynamic_html = fig.to_html(full_html=False)
        output_html_plots.update({nth_read_name: dynamic_html})

    return output_html_plots
