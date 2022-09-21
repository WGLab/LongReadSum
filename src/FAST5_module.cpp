/*
F5_module.cpp:
Class for calling FAST5 statistics modules.

*/
#include <array>
#include <algorithm>    // std::sort, copy
#include <numeric>      // std::accumulate
#include <iostream>
#include <vector>  // std::begin, std::end
#include <string>
#include <sstream>
#include <fstream>  // std::ofstream

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include "FAST5_module.h"
#include "ComFunction.h"
#include "H5Cpp.h"

using namespace H5;


std::vector<std::string> splitString(const std::string& str)
{
    std::vector<std::string> tokens;

    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, '\n')) {
        tokens.push_back(token);
    }

    return tokens;
}


// Read the FASTQ dataset from the file
std::vector<std::string> getFastq(H5::H5File f5, std::string basecall_group)
{
    // Attempt to access the FASTQ dataset
    H5::DataSet dataset;
    dataset = f5.openDataSet("/Analyses/"+basecall_group+"/BaseCalled_template/Fastq");
    H5::DataType  datatype   = dataset.getDataType();

    // Read output
    std::string data;
    dataset.read(data, datatype);
    std::vector<std::string> tokens = splitString(data);

    return tokens;
}


// Add base and read QC to the output data structure adn the output details file
static int writeBaseQCDetails(const char *input_file, Output_FAST5 &output_data, FILE *read_details_fp)
{
    int exit_code = 0;
    const char * read_name;
    char baseq;
    double gc_content_pct;
    std::string basecall_group = "Basecall_1D_000";

    // Access output QC and base quality QC structures
    Basic_Seq_Statistics &long_read_info = output_data.long_read_info;
    Basic_Seq_Quality_Statistics &seq_quality_info = output_data.seq_quality_info;

    // Run QC on the HDF5 file
    H5::Exception::dontPrint();  // Disable error printing
    try {
        // Open the file
        H5::H5File f5 = H5::H5File(input_file, H5F_ACC_RDONLY);

        // Get the FASTQ dataset
        std::vector<std::string> fq = getFastq(f5, basecall_group);

        // Access the read name
        std::string header_str = fq[0];
        std::istringstream iss_header( header_str );
        std::string read_name_str;
        std::getline( iss_header, read_name_str, ' ' );
        read_name = read_name_str.c_str();

        // Access the sequence data
        std::string sequence_data_str = fq[1];

        // Update the total number of bases
        int base_count = sequence_data_str.length();
        long_read_info.total_num_bases += base_count;

        // Store the read length
        long_read_info.read_lengths.push_back(base_count);

        // Access base quality data
        char value;
        std::vector<int> base_quality_values;
        std::string base_quality_str = fq[3];
        std::istringstream iss( base_quality_str );
        while (iss >> value) {
            int base_quality_value = value - '!';  // '!' symbol represent 0-quality score
            base_quality_values.push_back(base_quality_value);
        }

        // Update the base quality and GC content information
        int gc_count = 0;
        double read_mean_base_qual = 0;
        char current_base;
        for (int i = 0; i < base_count; i++)
        {
            current_base = sequence_data_str[i];
            if (current_base == 'A' || current_base == 'a')
            {
                long_read_info.total_a_cnt += 1;
            }
            else if (current_base == 'G' || current_base == 'g')
            {
                long_read_info.total_g_cnt += 1;
                gc_count += 1;
            }
            else if (current_base == 'C' || current_base == 'c')
            {
                long_read_info.total_c_cnt += 1;
                gc_count += 1;
            }
            else if (current_base == 'T' || current_base == 't' || current_base == 'U' || current_base == 'u')
            {
                long_read_info.total_tu_cnt += 1;
            }
            baseq = base_quality_values[i];  // Get the base quality
            seq_quality_info.base_quality_distribution[baseq] += 1;
            read_mean_base_qual += baseq;
        }

        // Calculate percent guanine & cytosine
        gc_content_pct = 100.0 *( (double)gc_count / (double)base_count );

        // Look into this section
        long_read_info.read_gc_content_count[(int)(gc_content_pct + 0.5)] += 1;
        read_mean_base_qual /= (double) base_count;
        seq_quality_info.read_average_base_quality_distribution[(uint)(read_mean_base_qual + 0.5)] += 1;
        fprintf(read_details_fp, "%s\t%d\t%.2f\t%.2f\n", read_name, base_count, gc_content_pct, read_mean_base_qual);

        // Update the total number of reads (1 per FAST5 file)
        output_data.long_read_info.total_num_reads += 1;
    } catch (FileIException &error) {
        // If the dataset is missing, continue and ignore this file
        if (error.getFuncName() == "H5File::openDataSet") {
            std::cout << "No FASTQ sequence dataset found for file: " << input_file << std::endl;
        } else {
            std::cerr << "Error accessing the FAST5 file: " << input_file << std::endl;
            error.printErrorStack();
            exit_code = 2;
        }
    };

    return exit_code;
}


// Add read signal QC to the output data structure and the output details file
static int writeSignalQCDetails(const char *input_file, Output_FAST5 &output_data)
{
    int exit_code = 0;
    std::string basecall_group = "Basecall_1D_000";

//    // Open the CSV files
//    std::ofstream raw_csv;
//    raw_csv.open(signal_raw_csv);
//    std::ofstream qc_csv;
//    qc_csv.open(signal_qc_csv);

    // Run QC on the HDF5 file
    //H5::Exception::dontPrint();  // Disable error printing
    try {
        // Open the file
        H5::H5File f5 = H5::H5File(input_file, H5F_ACC_RDONLY);

        // Get the read name
        std::string signal_group_str = "/Raw/Reads";
        H5::Group signal_group = f5.openGroup(signal_group_str);
        std::string read_name;
        read_name = signal_group.getObjnameByIdx(0);

        // Get the sequence string
        std::vector<std::string> fq = getFastq(f5, basecall_group);
        std::string sequence_data_str = fq[1];

        // Format the read signal dataset name
        std::ostringstream ss;
        ss << signal_group_str << "/" << read_name << "/Signal";
        std::string signal_dataset_str = ss.str();

        // Get the signal dataset
        H5::DataSet signal_ds = f5.openDataSet(signal_dataset_str);
        H5::DataType mdatatype= signal_ds.getDataType();
        H5::DataSpace dataspace  = signal_ds.getSpace();
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, NULL); // rank = 1
        int data_count = dims[0];

        // Store the signals in an array
        int f5signals [data_count];
        signal_ds.read(f5signals, mdatatype);

        // Get the block stride (window length) attribute
        H5::Group group_obj = f5.openGroup("/Analyses/Basecall_1D_000/Summary/basecall_1d_template");
        H5::Attribute block_stride_obj = group_obj.openAttribute("block_stride");
        H5::DataType bs_datatype= block_stride_obj.getDataType();
        int bs_buffer [1];
        block_stride_obj.read(bs_datatype, bs_buffer);
        int block_stride_value = bs_buffer[0];

        // Get the sequence length attribute
        H5::Attribute seq_length_obj = group_obj.openAttribute("sequence_length");
        H5::DataType sl_datatype= seq_length_obj.getDataType();
        int sl_buffer [1];
        seq_length_obj.read(sl_datatype, sl_buffer);
        //int sequence_length = sl_buffer[0];

        // Get the raw signal basecall start position
        H5::Group segm_group_obj = f5.openGroup("/Analyses/Segmentation_000/Summary/segmentation");
        H5::Attribute start_attr_obj = segm_group_obj.openAttribute("first_sample_template");
        H5::DataType st_datatype= start_attr_obj.getDataType();
        int st_buffer [1];
        start_attr_obj.read(st_datatype, st_buffer);
        int start_index = st_buffer[0];

        // Get the boolean array of base calls (move)
        H5::DataSet move_dataset_obj;
        move_dataset_obj = f5.openDataSet("/Analyses/" + basecall_group + "/BaseCalled_template/Move");
        H5::DataType  move_datatype = move_dataset_obj.getDataType();
        H5::DataSpace move_dataspace  = move_dataset_obj.getSpace();
        hsize_t move_dims[2];
        move_dataspace.getSimpleExtentDims(move_dims, NULL); // rank = 1
        int move_data_count = move_dims[0];

        // Read the boolean array
        uint8_t move_bool [move_data_count];
        move_dataset_obj.read(move_bool, move_datatype, move_dataspace);

        // Segment the raw signal by the window length
        std::vector<std::vector<int>> basecall_signals;
        int basecall_index = 0;
        int base_start_index = start_index;
        int base_count = 0;
        for (int i = 0; i < move_data_count; i++)
        {
            int move_value(move_bool[i]);
            if (move_value == 1)
            {
                base_count++;

                // Grab the signal
                std::vector<int> called_base_signal(block_stride_value);
                int end_index = base_start_index + block_stride_value;
                called_base_signal.assign(f5signals + base_start_index, f5signals + end_index);

                // Store in the 2D array
                basecall_signals.push_back(called_base_signal);
            }

            // Update indices
            base_start_index += block_stride_value;
            basecall_index ++;
        }

        // Append the basecall signals to the output structure
        Base_Signals basecall_obj(read_name, sequence_data_str, basecall_signals);
        output_data.addReadBaseSignals(basecall_obj);

    // catch failure caused by the H5File operations
    }
    catch (FileIException &error) {
        error.printErrorStack();
        exit_code = 2;
    }

    // catch failure caused by the DataSet operations
    catch (DataSetIException &error) {
        error.printErrorStack();
        exit_code = 2;
    }

    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException &error) {
        error.printErrorStack();
        exit_code = 2;
    }

    // catch failure caused by the Attribute operations
    catch (AttributeIException &error) {
        error.printErrorStack();
        exit_code = 2;
    }

    // Other
    catch (std::exception& e) {
        std::cerr << "Exception caught : " << e.what() << std::endl;
    }

//    // Close the CSV files
//    raw_csv.close();
//    qc_csv.close();

    return exit_code;
}


// Generate QC data for either base calls or signal data
int generateQCForFAST5(Input_Para &_input_data, Output_FAST5 &output_data)
{
    int exit_code = 0;
    const char *input_file = NULL;
    std::string read_details_file, read_summary_file;
    FILE *read_details_fp, *read_summary_fp;

    // Determine which statistics (base or signal) to generate
    bool signal_mode = (bool)_input_data.other_flags;
    std::cout << "FAST5 mode: " << ( signal_mode == true? "Signal" : "Base" ) << " QC" << std::endl;

    if (signal_mode == true) {
//        // Generate the signal data QC output
//        std::string signal_raw_csv(_input_data.output_folder + "/FAST5_signal_raw.csv");
//        std::string signal_qc_csv(_input_data.output_folder + "/FAST5_signal_QC.csv");

        // Loop through each input file and get the QC data across files
        size_t file_count = _input_data.num_input_files;
        for (size_t i = 0; i < file_count; i++)
        {
            input_file = _input_data.input_files[i].c_str();
            std::cout << "File name: " << input_file << std::endl;

            // Write QC details to the file
            exit_code = writeSignalQCDetails(input_file, output_data);
        }
    } else {
        // Generate the usual read and base QC output
        read_details_file = _input_data.output_folder + "/FAST5_details.txt";
        read_summary_file = _input_data.output_folder + "/FAST5_summary.txt";

        output_data.long_read_info.total_num_reads = ZeroDefault; // total number of long reads
        output_data.long_read_info.total_num_bases = ZeroDefault; // total number of bases

        output_data.long_read_info.longest_read_length = ZeroDefault; // the length of longest reads
        output_data.long_read_info.n50_read_length = MoneDefault;     // N50
        output_data.long_read_info.n95_read_length = MoneDefault;     // N95
        output_data.long_read_info.n05_read_length = MoneDefault;     // N05;
        output_data.long_read_info.mean_read_length = MoneDefault;    // mean of read length

        output_data.long_read_info.NXX_read_length.clear();
        output_data.long_read_info.median_read_length = MoneDefault; // median of read length

        output_data.long_read_info.total_a_cnt = ZeroDefault;  // A content
        output_data.long_read_info.total_c_cnt = ZeroDefault;  // C content
        output_data.long_read_info.total_g_cnt = ZeroDefault;  // G content
        output_data.long_read_info.total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
        output_data.long_read_info.total_n_cnt = ZeroDefault;  // N content
        output_data.long_read_info.gc_cnt = ZeroDefault;       // GC ratio

        //int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

        output_data.long_read_info.read_gc_content_count.clear();
        output_data.long_read_info.read_length_count.clear();
        output_data.seq_quality_info.base_quality_distribution.clear();
        output_data.seq_quality_info.read_average_base_quality_distribution.clear();

        output_data.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
        // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

        output_data.long_read_info.read_gc_content_count.resize(101, 0);
        // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

        output_data.long_read_info.NXX_read_length.resize(101, 0);
        // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;

        output_data.seq_quality_info.base_quality_distribution.resize(256, 0);
        // base_quality_distribution[x] means number of bases that quality = x.

        output_data.seq_quality_info.read_average_base_quality_distribution.resize(256, 0);
        // base_quality_distribution[x] means number of reads that average base quality = x.

        // Set up the output summary text file
        read_details_fp = fopen(read_details_file.c_str(), "w");
        if (NULL == read_details_fp)
        {
            std::cerr << "Failed to set up output file: " << read_details_file << std::endl;
            exit_code = 3;
        } else {
            fprintf(read_details_fp, "#read_name\tlength\tGC\taverage_baseq_quality\n");

            // Loop through each input file and get the QC data across files
            size_t file_count = _input_data.num_input_files;

            // Write QC details to the file
            for (size_t i = 0; i < file_count; i++)
            {
                input_file = _input_data.input_files[i].c_str();
                exit_code = writeBaseQCDetails(input_file, output_data, read_details_fp);
            }
            fclose(read_details_fp);

            // Check if the GC content was calculated successfully
            if (exit_code == 0) {

                // Add the G + C bases
                double g_c = output_data.long_read_info.total_g_cnt + output_data.long_read_info.total_c_cnt;

                // Add all bases
                double a_tu_g_c = g_c + output_data.long_read_info.total_a_cnt + output_data.long_read_info.total_tu_cnt;

                // Check that our total base counts match what was stored (That our code works)
                int total_num_bases = output_data.long_read_info.total_num_bases;
                if (a_tu_g_c != (double)total_num_bases)
                {
                    std::cerr << "Total number of bases is not consistent." << std::endl;
                    exit_code = 4;
                } else {
                    // Calculate GC-content
                    output_data.long_read_info.gc_cnt = g_c / a_tu_g_c;

                    // Sort the read lengths in descending order
                    std::vector<int> read_lengths = output_data.long_read_info.read_lengths;
                    std::sort(read_lengths.begin(), read_lengths.end(), std::greater<int>());

                    // Get the max read length
                    int max_read_length = *std::max_element(read_lengths.begin(), read_lengths.end());
                    output_data.long_read_info.longest_read_length = max_read_length;

                    // Get the median read length
                    int median_read_length = read_lengths[read_lengths.size() / 2];
                    output_data.long_read_info.median_read_length = median_read_length;

                    // Get the mean read length
                    float mean_read_length = (double)total_num_bases / (double)read_lengths.size();
                    output_data.long_read_info.mean_read_length = mean_read_length;

                    // Calculate N50 and other N-scores
                    for (int percent_value = 1; percent_value <= 100; percent_value++)
                    {
                        // Get the base percentage threshold for this N-score
                        double base_threshold = (double)total_num_bases * (percent_value / 100.0);

                        // Calculate the NXX score
                        double current_base_count = 0;
                        int current_read_index = -1;
                        while (current_base_count < base_threshold) {
                            current_read_index ++;
                            current_base_count += read_lengths[current_read_index];
                        }
                        int nxx_read_length = read_lengths[current_read_index];
                        output_data.long_read_info.NXX_read_length[percent_value] = nxx_read_length;
                    }
                    // Set common score variables
                    output_data.long_read_info.n50_read_length = output_data.long_read_info.NXX_read_length[50];
                    output_data.long_read_info.n95_read_length = output_data.long_read_info.NXX_read_length[95];
                    output_data.long_read_info.n05_read_length = output_data.long_read_info.NXX_read_length[5];

                    // Create the summary file
                    std::cout << "Writing summary file: " << read_summary_file.c_str() << std::endl;
                    read_summary_fp = fopen(read_summary_file.c_str(), "w");
                    fprintf(read_summary_fp, "total number of reads\t%ld\n", output_data.long_read_info.total_num_reads);
                    fprintf(read_summary_fp, "total number of bases\t%ld\n", output_data.long_read_info.total_num_bases);
                    fprintf(read_summary_fp, "longest read length\t%lu\n", output_data.long_read_info.longest_read_length);
                    fprintf(read_summary_fp, "N50 read length\t%ld\n", output_data.long_read_info.n50_read_length);
                    fprintf(read_summary_fp, "mean read length\t%.2f\n", output_data.long_read_info.mean_read_length);
                    fprintf(read_summary_fp, "median read length\t%ld\n", output_data.long_read_info.median_read_length);
                    fprintf(read_summary_fp, "GC%%\t%.2f\n", output_data.long_read_info.gc_cnt * 100);
                    fprintf(read_summary_fp, "\n\n");
                    for (int percent = 5; percent < 100; percent += 5)
                    {
                        fprintf(read_summary_fp, "N%02d read length\t%.ld\n", percent, output_data.long_read_info.NXX_read_length[percent]);
                    }

                    fprintf(read_summary_fp, "\n\n");

                    fprintf(read_summary_fp, "GC content\tnumber of reads\n");
                    for (int gc_ratio = 0; gc_ratio < 100; gc_ratio++)
                    {
                        fprintf(read_summary_fp, "GC=%d%%\t%ld\n", gc_ratio, output_data.long_read_info.read_gc_content_count[gc_ratio]);
                    }


                    fprintf(read_summary_fp, "\n\n");
                    fprintf(read_summary_fp, "base quality\tnumber of bases\n");
                    for (int baseq = 0; baseq <= 60; baseq++)
                    {
                        fprintf(read_summary_fp, "%d\t%ld\n", baseq, output_data.seq_quality_info.base_quality_distribution[baseq]);
                    }

                    fprintf(read_summary_fp, "\n\n");
                    fprintf(read_summary_fp, "read average base quality\tnumber of reads\n");
                    for (int baseq = 0; baseq <= 60; baseq++)
                    {
                        fprintf(read_summary_fp, "%d\t%ld\n", baseq, output_data.seq_quality_info.read_average_base_quality_distribution[baseq]);
                    }
                    fclose(read_summary_fp);
                }
            }
        }
    }

    return exit_code;
}
