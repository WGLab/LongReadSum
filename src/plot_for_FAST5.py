"""
plot_for_FAST%.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

from src import lrst_global
from src.utils import *


def generate_bs(f5_output, para_dict):
    lrst_global.plot_filenames["basic_st"] = {};
    lrst_global.plot_filenames["basic_st"]['file'] = ""
    lrst_global.plot_filenames["basic_st"]['title'] = "Table of Statistics Summary"
    lrst_global.plot_filenames["basic_st"]['description'] = "Basic Statistics for F5 summary"

    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Statistics</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
    table_str += int_str_for_format.format("#Total Reads", \
                                           f5_output.long_read_info.total_num_reads);
    table_str += int_str_for_format.format("#Total Bases", \
                                           f5_output.long_read_info.total_num_bases);
    table_str += int_str_for_format.format("Longest Read Length", \
                                           f5_output.long_read_info.longest_read_length);
    table_str += int_str_for_format.format("N50", \
                                           f5_output.long_read_info.n50_read_length);
    table_str += double_str_for_format.format("GC Content(%)", \
                                              f5_output.long_read_info.gc_cnt * 100);
    table_str += double_str_for_format.format("Mean Read Length", \
                                              f5_output.long_read_info.mean_read_length);
    table_str += int_str_for_format.format("Median Read Length", \
                                           f5_output.long_read_info.median_read_length);
    table_str += "\n</tbody>\n</table>"

    lrst_global.plot_filenames["basic_st"]['detail'] = table_str;


def plot(f5_output, para_dict):
    out_path = para_dict["output_folder"]
    get_image_path = lambda x: os.path.join(out_path, lrst_global.plot_filenames[x]['file'])

    # print("num_primary_alignment: {}".format(f5_output.num_primary_alignment))
    generate_bs(f5_output, para_dict)

    plot_read_length_stats([f5_output.long_read_info], get_image_path('read_length_st'), subtitles=['Long Reads'])
    plot_base_counts([f5_output.long_read_info], get_image_path('base_st'), subtitles=['Long Reads'])
    plot_basic_info([f5_output.long_read_info], get_image_path('basic_info'), categories=['Long Reads'])

    lrst_global.plot_filenames['read_length_hist']['dynamic'] = histogram(f5_output.long_read_info,
                                                                          get_image_path('read_length_hist'))
    lrst_global.plot_filenames['base_quality']['dynamic'] = base_quality(f5_output.seq_quality_info,
                                                                         get_image_path('base_quality'))
    lrst_global.plot_filenames['read_avg_base_quality']['dynamic'] = read_avg_base_quality(f5_output.seq_quality_info,
                                                                                           get_image_path(
                                                                                               'read_avg_base_quality'))