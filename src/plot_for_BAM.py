"""
plot_for_BAM.py:
Use the formatted statistics from our C++ module output text files to generate summary plots in image format.
"""

from src import lrst_global
from src.utils import *

def plot_alignment_numbers(data, path):
    fig, axes = plt.subplots(figsize =(8, 6))
    
    numbers_list=[[data.num_primary_alignment, data.num_supplementary_alignment, data.num_secondary_alignment,
                   data.num_reads_with_supplementary_alignment, data.num_reads_with_secondary_alignment,
                   data.num_reads_with_both_secondary_supplementary_alignment, data.forward_alignment, data.reverse_alignment]]
    
    category=['Primary Alignments', 'Supplementary Alignments', 'Secondary Alignments', 'Reads with Supplementary Alignments','Reads with Secondary Alignments','Reads with Secondary and Supplementary Alignments', 'Forward Alignments', 'Reverse Alignments']
    category=[wrap(x) for x in category]
    
    category_list=itertools.cycle([category])
    xlabel_list=itertools.cycle(['Counts'])
    ylabel_list=itertools.cycle([''])
    subtitle_list=[None]
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='h')
    
def plot_errors(bam_output, path):
    fig, axes = plt.subplots(1,1, figsize =(8, 6))

    numbers_list=[[bam_output.num_matched_bases, bam_output.num_mismatched_bases, bam_output.num_ins_bases, bam_output.num_del_bases, bam_output.num_clip_bases]]
    
    category=['Matched Bases', 'Mismatched Bases', 'Inserted Bases', 'Deleted Bases', 'Clipped Bases']
    category=[wrap(x) for x in category]
    
    category_list=itertools.cycle([category])
    xlabel_list=itertools.cycle(['Counts'])
    ylabel_list=itertools.cycle([None])
    subtitle_list=[None]
    bar_plot(fig, numbers_list, category_list, xlabel_list, ylabel_list, subtitle_list, path, orientation='h')

def generate_bs( bam_output, para_dict ):
    lrst_global.plot_filenames["basic_st"] = {};
    lrst_global.plot_filenames["basic_st"]['file'] = ""
    lrst_global.plot_filenames["basic_st"]['title'] = "Basic statistics"
    lrst_global.plot_filenames["basic_st"]['description'] = "BAM: Basic statistics"

    table_str = "<table>\n<thead>\n<tr><th>Measurement</th><th>Mapped</th><th>Unmapped</th><th>All</th></tr>\n</thead>"
    table_str += "\n<tbody>"
    int_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td><td style=\"text-align:right\">{:,d}</td></tr>"
    double_str_for_format = "<tr><td>{}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td><td style=\"text-align:right\">{:.1f}</td></tr>"
    table_str += int_str_for_format.format("#Total Reads", \
                 bam_output.mapped_long_read_info.total_num_reads, \
                 bam_output.unmapped_long_read_info.total_num_reads, \
                 bam_output.long_read_info.total_num_reads);
    table_str += int_str_for_format.format("#Total Bases", \
                 bam_output.mapped_long_read_info.total_num_bases, \
                 bam_output.unmapped_long_read_info.total_num_bases, \
                 bam_output.long_read_info.total_num_bases);
    table_str += int_str_for_format.format("Longest Read Length", \
                 bam_output.mapped_long_read_info.longest_read_length, \
                 bam_output.unmapped_long_read_info.longest_read_length, \
                 bam_output.long_read_info.longest_read_length);
    table_str += int_str_for_format.format("N50", \
                 bam_output.mapped_long_read_info.n50_read_length, \
                 bam_output.unmapped_long_read_info.n50_read_length, \
                 bam_output.long_read_info.n50_read_length);
    table_str += double_str_for_format.format("GC Content(%)", \
                 bam_output.mapped_long_read_info.gc_cnt*100, \
                 bam_output.unmapped_long_read_info.gc_cnt*100, \
                 bam_output.long_read_info.gc_cnt*100);
    table_str += double_str_for_format.format("Mean Read Length", \
                 bam_output.mapped_long_read_info.mean_read_length, \
                 bam_output.unmapped_long_read_info.mean_read_length, \
                 bam_output.long_read_info.mean_read_length);
    table_str += int_str_for_format.format("Median Read Length", \
                 bam_output.mapped_long_read_info.median_read_length, \
                 bam_output.unmapped_long_read_info.median_read_length, \
                 bam_output.long_read_info.median_read_length);
    table_str += "\n</tbody>\n</table>"

    lrst_global.plot_filenames["basic_st"]['detail'] = table_str;
    

def bam_plot( bam_output, para_dict ):
    
    out_path=para_dict["output_folder"]
    get_image_path=lambda x: os.path.join(out_path,lrst_global.plot_filenames[x]['file'])

    # Set the default matplotlib font size
    setDefaultFontSize(12)

    # Get the font size for plotly plots
    font_size = para_dict["fontsize"]

    # Create table
    generate_bs(bam_output, para_dict)

    # Generate plots
    plot_alignment_numbers(bam_output, get_image_path('map_st'))
    plot_errors(bam_output, get_image_path('err_st'))
    
    plot_read_length_stats([bam_output.long_read_info, bam_output.mapped_long_read_info, bam_output.unmapped_long_read_info], get_image_path('read_length_st'), subtitles=['All Reads', 'Mapped Reads', 'Unmapped Reads'])
    plot_base_counts([bam_output.long_read_info, bam_output.mapped_long_read_info, bam_output.unmapped_long_read_info], get_image_path('base_st'), subtitles=['All Reads', 'Mapped Reads', 'Unmapped Reads'])
    plot_basic_info([bam_output.long_read_info, bam_output.mapped_long_read_info, bam_output.unmapped_long_read_info], get_image_path('basic_info'), categories=['All Reads', 'Mapped Reads', 'Unmapped Reads'])
    
    lrst_global.plot_filenames['read_length_hist']['dynamic'] = histogram(bam_output.long_read_info, get_image_path('read_length_hist'), font_size)
    lrst_global.plot_filenames['base_quality']['dynamic'] = base_quality(bam_output.seq_quality_info, get_image_path('base_quality'), font_size)
