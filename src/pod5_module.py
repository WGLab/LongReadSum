import os
import numpy as np
import logging
import pod5 as p5

# Import the longreadsum module
if __name__ == 'src.pod5_module':
    from lib import lrst  # For debugging
else:
    import lrst

def generate_pod5_qc(input_data: dict) -> dict:
    """Generate QC for POD5"""
    logging.info("Generating QC for POD5")

    # Get the list of read IDs to process (if specified)
    read_id_list = []
    if (input_data['read_ids']):
        # Parse the comma-separated list of read IDs
        read_id_list = input_data['read_ids'].split(',')
        logging.info("Processing %d read IDs", len(read_id_list))

    # Generate signal QC
    # Loop through each input file and get the QC data across files
    # Dictionary to store read signal QC (key: read_id, [signal, mean_signal, median_signal, std_signal, skewness, kurtosis])
    read_signal_dict = {}
    for input_file in input_data['input_files']:
        logging.info("File name: %s", input_file)

        # Iterate through each read in the file
        with p5.Reader(input_file) as reader:
            for read in reader:
                # Check if the read is in the list of reads to process (if
                # provided)
                # Convert the UUID to a string
                read_id = str(read.read_id)
                if read_id_list and read_id not in read_id_list:
                    logging.info("Skipping read ID: %s", read_id)
                    continue
                logging.info("Processing read ID: %s", read_id)

                # Get the basecall signals
                read_signal = read.signal

                read_signal_dict[read_id] = {}
                # logging.info("Adding read signals")
                read_signal_dict[read_id]['signal'] = read_signal

                # Calculate QC metrics using numpy
                # Mean signal
                mean_signal = np.mean(read_signal)
                read_signal_dict[read_id]['mean'] = mean_signal

                # Median signal
                median_signal = np.median(read_signal)
                read_signal_dict[read_id]['median'] = median_signal

                # Standard deviation
                std_signal = np.std(read_signal)
                read_signal_dict[read_id]['std'] = std_signal

                # Pearson skewness coefficient
                moment3 = np.mean((read_signal - mean_signal) ** 3)
                moment2 = np.mean((read_signal - mean_signal) ** 2)
                skewness = moment3 / (moment2 ** 1.5)
                read_signal_dict[read_id]['skewness'] = skewness

                # Kurtosis
                moment4 = np.mean((read_signal - mean_signal) ** 4)
                kurtosis = moment4 / (moment2 ** 2)
                read_signal_dict[read_id]['kurtosis'] = kurtosis

                
    return read_signal_dict

        
    # Write QC details to the file
    

    

            
    # try {
    #     // Open the file
    #     H5::H5File f5 = H5::H5File(input_file, H5F_ACC_RDONLY);

    #     // Check if it is a multi-read FAST5 file
    #     std::string signal_group_str;
    #     std::string read_name;
    #     if (f5.nameExists("/Raw")) {
    #         std::cout << "Single read FAST5" << std::endl;

    #         // Append the basecall signals to the output structure
    #         signal_group_str = "/Raw/Reads";
    #         read_name = getFileReadName(f5);
    #         std::cout << read_name << std::endl;
    #         Base_Signals basecall_obj = getReadBaseSignalData(f5, read_name, true);
    #         output_data.addReadBaseSignals(basecall_obj);
    #     } else {
    #         std::cout << "Multi-read FAST5" << std::endl;

    #         // Loop through each read
    #         std::cout << "Reading all reads" << std::endl;
    #         H5::Group root_group = f5.openGroup("/");
    #         size_t num_objs = root_group.getNumObjs();
    #         for (size_t i=0; i < num_objs; i++) {
    #             read_name = root_group.getObjnameByIdx(i);

    #             // Check if the read is in the list of reads to process (if provided)
    #             if (read_id_list.size() > 0) {
    #                 // First remove the prefix
    #                 std::string read_id = read_name.substr(5);
    #                 if (std::find(read_id_list.begin(), read_id_list.end(), read_id) == read_id_list.end()) {
    #                     //std::cout << "Skipping read ID: " << read_id << std::endl;
    #                     continue;
    #                 } else {
    #                     std::cout << "Processing read ID: " << read_id << std::endl;
    #                 }
    #             }
    #             // std::cout << "Read: " << read_name << std::endl;

    #             // Get the basecall signals
    #             // std::cout << "Getting basecall signals" << std::endl;
    #             Base_Signals basecall_obj = getReadBaseSignalData(f5, read_name, false);

    #             //std::cout << "Adding basecall signals" << std::endl;
    #             output_data.addReadBaseSignals(basecall_obj);
    #         }
    #     }

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

    