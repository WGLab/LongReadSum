"""
plot_for_SeqTxt.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

if __package__ == 'src':
    from src.plot_utils import *
else:
    from plot_utils import *


def create_summary_table(seqtxt_output, plot_filepaths):
    plot_filepaths["basic_st"] = {}
    plot_filepaths["basic_st"]['file'] = ""
    plot_filepaths["basic_st"]['title'] = "Basic statistics"
    plot_filepaths["basic_st"]['description'] = "Sequencing Summary: Basic statistics"

    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Passed</th><th>Failed</th><th>All</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
    table_str += int_str_for_format.format("#Total Reads",
                                           seqtxt_output.passed_long_read_info.long_read_info.total_num_reads,
                                           seqtxt_output.failed_long_read_info.long_read_info.total_num_reads,
                                           seqtxt_output.all_long_read_info.long_read_info.total_num_reads)
    table_str += int_str_for_format.format("#Total Bases",
                                           seqtxt_output.passed_long_read_info.long_read_info.total_num_bases,
                                           seqtxt_output.failed_long_read_info.long_read_info.total_num_bases,
                                           seqtxt_output.all_long_read_info.long_read_info.total_num_bases)
    table_str += int_str_for_format.format("Longest Read Length",
                                           seqtxt_output.passed_long_read_info.long_read_info.longest_read_length,
                                           seqtxt_output.failed_long_read_info.long_read_info.longest_read_length,
                                           seqtxt_output.all_long_read_info.long_read_info.longest_read_length)
    table_str += int_str_for_format.format("N50",
                                           seqtxt_output.passed_long_read_info.long_read_info.n50_read_length,
                                           seqtxt_output.failed_long_read_info.long_read_info.n50_read_length,
                                           seqtxt_output.all_long_read_info.long_read_info.n50_read_length)
    table_str += double_str_for_format.format("Mean Read Length",
                                              seqtxt_output.passed_long_read_info.long_read_info.mean_read_length,
                                              seqtxt_output.failed_long_read_info.long_read_info.mean_read_length,
                                              seqtxt_output.all_long_read_info.long_read_info.mean_read_length)
    table_str += int_str_for_format.format("Median Read Length",
                                           seqtxt_output.passed_long_read_info.long_read_info.median_read_length,
                                           seqtxt_output.failed_long_read_info.long_read_info.median_read_length,
                                           seqtxt_output.all_long_read_info.long_read_info.median_read_length)
    table_str += "\n</tbody>\n</table>"

    plot_filepaths["basic_st"]['detail'] = table_str


def plot(seqtxt_output, para_dict):
    out_path = para_dict["output_folder"]
    plot_filepaths = getDefaultPlotFilenames()
    get_image_path = lambda x: os.path.join(out_path, plot_filepaths[x]['file'])

    # Set the default matplotlib font size
    setDefaultFontSize(12)

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create table
    create_summary_table(seqtxt_output, plot_filepaths)

    # Generate plots
    plot_read_length_stats(
        [seqtxt_output.all_long_read_info.long_read_info, seqtxt_output.passed_long_read_info.long_read_info,
         seqtxt_output.failed_long_read_info.long_read_info], get_image_path('read_length_st'),
        subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    plot_base_counts(
        [seqtxt_output.all_long_read_info.long_read_info, seqtxt_output.passed_long_read_info.long_read_info,
         seqtxt_output.failed_long_read_info.long_read_info], get_image_path('base_st'),
        subtitles=['All Reads', 'Passed Reads', 'Failed Reads'])
    plot_basic_info(
        [seqtxt_output.all_long_read_info.long_read_info, seqtxt_output.passed_long_read_info.long_read_info,
         seqtxt_output.failed_long_read_info.long_read_info], get_image_path('basic_info'),
        categories=['All Reads', 'Passed Reads', 'Failed Reads'])
    histogram(seqtxt_output.all_long_read_info.long_read_info, get_image_path('read_length_hist'), font_size)

    return plot_filepaths
