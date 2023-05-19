"""
plot_for_BAM.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

if __package__ == 'src':
    from src.plot_utils import *
else:
    from plot_utils import *


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
    layout = go.Layout(title=go.layout.Title(text=""), xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Counts")), yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Generate the HTML object for the plot
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=1000)

    return html_obj


def plot_errors(bam_output, path):
    category = ['Matched Bases', 'Mismatched Bases', 'Inserted Bases', 'Deleted Bases', 'Clipped Bases']
    category = [wrap(x) for x in category]

    # Create a horizontally aligned bar plot trace from the data using plotly
    trace = go.Bar(x=[bam_output.num_matched_bases, bam_output.num_mismatched_bases, bam_output.num_ins_bases,
                        bam_output.num_del_bases, bam_output.num_clip_bases], y=category, orientation='h')

    # Create the layout for the plot
    layout = go.Layout(title=go.layout.Title(text=""), xaxis=go.layout.XAxis(title=go.layout.xaxis.Title(text="Counts")), yaxis=go.layout.YAxis(title=go.layout.yaxis.Title(text="")))

    # Create the figure object
    fig = go.Figure(data=[trace], layout=layout)

    # Generate the HTML object for the plot
    html_obj = fig.to_html(full_html=False, default_height=500, default_width=700)

    return html_obj


def create_summary_table(bam_output, plot_filepaths):
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Basic Statistics"
    plot_filepaths["basic_st"]['description'] = "BAM: Basic Statistics"

    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Mapped</th><th>Unmapped</th><th>All</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:," \
                         "d}</td><td style=\"text-align:right\">{:,d}</td></tr> "
    double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td " \
                            "style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr> "
    table_str += int_str_for_format.format("#Total Reads", bam_output.mapped_long_read_info.total_num_reads,
                                           bam_output.unmapped_long_read_info.total_num_reads,
                                           bam_output.long_read_info.total_num_reads)
    table_str += int_str_for_format.format("#Total Bases",
                                           bam_output.mapped_long_read_info.total_num_bases,
                                           bam_output.unmapped_long_read_info.total_num_bases,
                                           bam_output.long_read_info.total_num_bases)
    table_str += int_str_for_format.format("Longest Read Length",
                                           bam_output.mapped_long_read_info.longest_read_length,
                                           bam_output.unmapped_long_read_info.longest_read_length,
                                           bam_output.long_read_info.longest_read_length)
    table_str += int_str_for_format.format("N50",
                                           bam_output.mapped_long_read_info.n50_read_length,
                                           bam_output.unmapped_long_read_info.n50_read_length,
                                           bam_output.long_read_info.n50_read_length)
    table_str += double_str_for_format.format("GC Content(%)",
                                              bam_output.mapped_long_read_info.gc_cnt * 100,
                                              bam_output.unmapped_long_read_info.gc_cnt * 100,
                                              bam_output.long_read_info.gc_cnt * 100)
    table_str += double_str_for_format.format("Mean Read Length",
                                              bam_output.mapped_long_read_info.mean_read_length,
                                              bam_output.unmapped_long_read_info.mean_read_length,
                                              bam_output.long_read_info.mean_read_length)
    table_str += int_str_for_format.format("Median Read Length",
                                           bam_output.mapped_long_read_info.median_read_length,
                                           bam_output.unmapped_long_read_info.median_read_length,
                                           bam_output.long_read_info.median_read_length)
    table_str += "\n</tbody>\n</table>"

    plot_filepaths["basic_st"]['detail'] = table_str


def plot(bam_output, para_dict):
    out_path = para_dict["output_folder"]
    plot_filepaths = getDefaultPlotFilenames()
    get_image_path = lambda x: os.path.join(out_path, plot_filepaths[x]['file'])

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create the summary table
    create_summary_table(bam_output, plot_filepaths)

    # Generate plots
    plot_filepaths['read_alignments_bar']['dynamic'] = plot_alignment_numbers(bam_output, get_image_path('read_alignments_bar'))
    plot_filepaths['base_alignments_bar']['dynamic'] = plot_errors(bam_output, get_image_path('base_alignments_bar'))
    plot_filepaths['base_counts']['dynamic'] = plot_base_counts(bam_output)
    plot_filepaths['basic_info']['dynamic'] = plot_basic_info(bam_output)
    plot_filepaths['read_length_hist']['dynamic'] = read_lengths_histogram(bam_output.long_read_info,
                                                              get_image_path('read_length_hist'),
                                                              font_size)

    plot_filepaths['base_quality']['dynamic'] = base_quality(bam_output.seq_quality_info,
                                                             get_image_path('base_quality'), font_size)

    plot_filepaths['read_length_bar']['dynamic'] = plot_read_length_stats(bam_output)

    return plot_filepaths
