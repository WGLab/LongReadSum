import os
import logging
import pod5 as p5

# Import the longreadsum module
if __name__ == 'src.pod5_module':
    from lib import lrst  # For debugging
else:
    import lrst

def generate_pod5_qc(input_data: lrst.Input_Para, output_data: lrst.Output_FAST5) -> int:
    """Generate QC for POD5"""
    exit_code = 0

    # Determine which statistics (base or signal) to generate
    signal_mode = (bool)(input_data.other_flags)
    logging.info("POD5 mode: %s QC", 'Signal' if signal_mode else 'Base')

    # # Get the list of read IDs to process (if specified)
    # std::vector<std::string> read_id_list;
    # if (_input_data.read_ids.empty() == false) {
    #     // Parse the comma-separated list of read IDs
    #     std::stringstream ss(_input_data.read_ids);
    #     std::string read_id;
    #     while (std::getline(ss, read_id, ',')) {
    #         read_id_list.push_back(read_id);
    #     }
    # }

    # if (signal_mode == true) {
    #     // Loop through each input file and get the QC data across files
    #     size_t file_count = _input_data.num_input_files;
    #     for (size_t i = 0; i < file_count; i++)
    #     {
    #         input_file = _input_data.input_files[i].c_str();
    #         std::cout << "File name: " << input_file << std::endl;

    #         // Write QC details to the file
    #         exit_code = writeSignalQCDetails(input_file, output_data, read_id_list);
    #     }
    # } else {
    #     // Generate the usual read and base QC output
    #     read_details_file = _input_data.output_folder + "/FAST5_details.txt";
    #     read_summary_file = _input_data.output_folder + "/FAST5_summary.txt";

    #     // Set up the output summary text file
    #     read_details_fp = fopen(read_details_file.c_str(), "w");
    #     if (NULL == read_details_fp)
    #     {
    #         std::cerr << "Failed to set up output file: " << read_details_file << std::endl;
    #         exit_code = 3;
    #     } else {
    #         fprintf(read_details_fp, "#read_name\tlength\tGC\taverage_baseq_quality\n");

    #         // Loop through each input file and get the QC data across files
    #         size_t file_count = _input_data.num_input_files;

    #         // Write QC details to the file
    #         for (size_t i = 0; i < file_count; i++)
    #         {
    #             input_file = _input_data.input_files[i].c_str();
    #             exit_code = writeBaseQCDetails(input_file, output_data, read_details_fp);
    #         }
    #         fclose(read_details_fp);

    #         // Check if the GC content was calculated successfully
    #         if (exit_code == 0) {

    #             // Add the G + C bases
    #             double g_c = output_data.long_read_info.total_g_cnt + output_data.long_read_info.total_c_cnt;

    #             // Add all bases
    #             double a_tu_g_c = g_c + output_data.long_read_info.total_a_cnt + output_data.long_read_info.total_tu_cnt;

    #             // Calculate read length statistics if base counts are not zero
    #             uint64_t total_num_bases = output_data.long_read_info.total_num_bases;
    #             if (total_num_bases == 0) {
    #                 std::cerr << "No bases found in input files." << std::endl;
    #                 exit_code = 3;
    #             } else {
    #                 // Calculate GC-content
    #                 output_data.long_read_info.gc_cnt = g_c / a_tu_g_c;

    #                 // Sort the read lengths in descending order
    #                 std::vector<int> read_lengths = output_data.long_read_info.read_lengths;
    #                 std::sort(read_lengths.begin(), read_lengths.end(), std::greater<int64_t>());

    #                 // Get the max read length
    #                 int64_t max_read_length = read_lengths.at(0);
    #                 output_data.long_read_info.longest_read_length = max_read_length;

    #                 // Get the median read length
    #                 int64_t median_read_length = read_lengths[read_lengths.size() / 2];
    #                 output_data.long_read_info.median_read_length = median_read_length;

    #                 // Get the mean read length
    #                 float mean_read_length = (double)total_num_bases / (double)read_lengths.size();
    #                 output_data.long_read_info.mean_read_length = mean_read_length;

    #                 // Calculate N50 and other N-scores
    #                 for (int percent_value = 1; percent_value <= 100; percent_value++)
    #                 {
    #                     // Get the base percentage threshold for this N-score
    #                     double base_threshold = (double)total_num_bases * (percent_value / 100.0);

    #                     // Calculate the NXX score
    #                     double current_base_count = 0;
    #                     int current_read_index = -1;
    #                     while (current_base_count < base_threshold) {
    #                         current_read_index ++;
    #                         current_base_count += read_lengths.at(current_read_index);
    #                     }
    #                     int nxx_read_length = read_lengths.at(current_read_index);
    #                     output_data.long_read_info.NXX_read_length[percent_value] = nxx_read_length;
    #                 }

    #                 // Set common score variables
    #                 output_data.long_read_info.n50_read_length = output_data.long_read_info.NXX_read_length[50];
    #                 output_data.long_read_info.n95_read_length = output_data.long_read_info.NXX_read_length[95];
    #                 output_data.long_read_info.n05_read_length = output_data.long_read_info.NXX_read_length[5];

    #                 // Create the summary file
    #                 std::cout << "Writing summary file: " << read_summary_file.c_str() << std::endl;
    #                 read_summary_fp = fopen(read_summary_file.c_str(), "w");
    #                 fprintf(read_summary_fp, "total number of reads\t%d\n", output_data.long_read_info.total_num_reads);
    #                 fprintf(read_summary_fp, "total number of bases\t%ld\n", output_data.long_read_info.total_num_bases);
    #                 fprintf(read_summary_fp, "longest read length\t%d\n", output_data.long_read_info.longest_read_length);
    #                 fprintf(read_summary_fp, "N50 read length\t%d\n", output_data.long_read_info.n50_read_length);
    #                 fprintf(read_summary_fp, "mean read length\t%.2f\n", output_data.long_read_info.mean_read_length);
    #                 fprintf(read_summary_fp, "median read length\t%d\n", output_data.long_read_info.median_read_length);
    #                 fprintf(read_summary_fp, "GC%%\t%.2f\n", output_data.long_read_info.gc_cnt * 100);
    #                 fprintf(read_summary_fp, "\n\n");
    #                 for (int percent = 5; percent < 100; percent += 5)
    #                 {
    #                     fprintf(read_summary_fp, "N%02d read length\t%.d\n", percent, output_data.long_read_info.NXX_read_length[percent]);
    #                 }

    #                 fprintf(read_summary_fp, "\n\n");

    #                 fprintf(read_summary_fp, "GC content\tnumber of reads\n");
    #                 for (int gc_ratio = 0; gc_ratio < 100; gc_ratio++)
    #                 {
    #                     fprintf(read_summary_fp, "GC=%d%%\t%d\n", gc_ratio, output_data.long_read_info.read_gc_content_count[gc_ratio]);
    #                 }


    #                 fprintf(read_summary_fp, "\n\n");
    #                 fprintf(read_summary_fp, "base quality\tnumber of bases\n");
    #                 for (int baseq = 0; baseq <= 60; baseq++)
    #                 {
    #                     fprintf(read_summary_fp, "%d\t%ld\n", baseq, output_data.seq_quality_info.base_quality_distribution[baseq]);
    #                 }

    #                 fprintf(read_summary_fp, "\n\n");
    #                 fprintf(read_summary_fp, "read average base quality\tnumber of reads\n");
    #                 for (int baseq = 0; baseq <= 60; baseq++)
    #                 {
    #                     fprintf(read_summary_fp, "%d\t%d\n", baseq, output_data.seq_quality_info.read_average_base_quality_distribution[baseq]);
    #                 }
    #                 fclose(read_summary_fp);
    #             }
    #         }
    #     }
    # }

    return exit_code
    