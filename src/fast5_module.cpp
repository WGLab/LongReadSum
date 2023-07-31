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

#include "fast5_module.h"
#include "ComFunction.h"
#include "H5Cpp.h"

using namespace H5;

// Split a string into tokens
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
    std::vector<std::string> sequence;  // Output header and sequence
    H5::Exception::dontPrint();  // Disable error printing
    try {
        // Attempt to access the FASTQ dataset
        //std::cout << "Opening FASTQ dataset..." << std::endl;
        std::string fastq_path = basecall_group + "/BaseCalled_template/Fastq";
        //std::cout << fastq_path << std::endl;
        H5::DataSet dataset;
        dataset = f5.openDataSet(fastq_path);

        // Read output
        H5::DataType  datatype   = dataset.getDataType();
        std::string data;
        //std::cout << "Reading FASTQ data..." << std::endl;
        dataset.read(data, datatype);
        //std::cout << "FASTQ data read." << std::endl;
        sequence = splitString(data);
        //std::cout << "FASTQ data split." << std::endl;
    } catch (FileIException &error) {
    }

    return sequence;
}

// Read alignment data from the file to get start and end map positions
std::vector<int> getMapPositions(H5::H5File f5, std::string read_name)
{
    std::vector<int> map_positions;  // Start and end map positions
    H5::Exception::dontPrint();  // Disable error printing
    try {
        // Get the alignment start position
        H5::Group align_group_obj;
        if (read_name == "") {
            // Single-read FAST5 file
            align_group_obj = f5.openGroup("/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment");
        } else {
            // Multi-read FAST5 file
            std::string align_group_path = "/" + read_name + "/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment";
            align_group_obj = f5.openGroup(align_group_path);
        }
        H5::Attribute start_map_obj = align_group_obj.openAttribute("mapped_start");
        H5::DataType start_map_datatype= start_map_obj.getDataType();
        int start_map_buffer [1];
        start_map_obj.read(start_map_datatype, start_map_buffer);
        int start_map = start_map_buffer[0];
        map_positions.push_back(start_map);
        std::cout << "Start map = " << start_map << std::endl;

        // Get the alignment end position
        H5::Attribute end_map_obj = align_group_obj.openAttribute("mapped_end");
        H5::DataType end_map_datatype= end_map_obj.getDataType();
        int end_map_buffer [1];
        end_map_obj.read(end_map_datatype, end_map_buffer);
        int end_map = end_map_buffer[0];
        map_positions.push_back(end_map);
        std::cout << "End map = " << end_map << std::endl;

    } catch (FileIException &error) {
    } catch (AttributeIException &error) {
    }

    return map_positions;
}

// Read alignment data from the file to get the mapped chromosome
std::string getChromosome(H5::H5File f5, std::string read_name)
{
    std::string mapped_chrom;  // Start and end map positions
    std::string error_msg("No mapping found.");  // Message when no alignment information is found.
    H5::Exception::dontPrint();  // Disable HDF5 error printing
    try {
        // Get the alignment start position
        std::cout << "Opening mapping group..." << std::endl;

        // Get the alignment start position
        H5::Group align_group_obj;
        if (read_name == "") {
            // Single-read FAST5 file
            align_group_obj = f5.openGroup("/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment");
        } else {
            // Multi-read FAST5 file
            std::string align_group_path = "/" + read_name + "/Analyses/RawGenomeCorrected_000/BaseCalled_template/Alignment";
            align_group_obj = f5.openGroup(align_group_path);
        }
        H5::Attribute chrom_obj = align_group_obj.openAttribute("mapped_chrom");
        H5::DataType chrom_datatype= chrom_obj.getDataType();
        std::string chrom_buffer;
        chrom_obj.read(chrom_datatype, chrom_buffer);
        std::cout << "Chromosome = " << chrom_buffer << std::endl;
        mapped_chrom = chrom_buffer;
    } catch (FileIException &error) {
        std::cout << error_msg << std::endl;
    } catch (AttributeIException &error) {
        std::cout << error_msg << std::endl;
    }

    return mapped_chrom;
}

// Get the read name for a single-read FAST5 file
std::string getFileReadName(H5::H5File f5) {
    std::string read_name;
    H5::Group signal_group = f5.openGroup("/Raw/Reads");
    read_name = signal_group.getObjnameByIdx(0);

    return read_name;
}

// Format read base signal data from the specified signal group
Base_Signals getReadBaseSignalData(H5::H5File f5, std::string read_name, bool single_read)
{
    std::string sequence_data_str("");
    std::vector<std::vector<int>> basecall_signals;  // Holds signals across bases
    std::string mapped_chrom;
    std::vector<int> map_positions;
    std::string read_info = read_name;

    // Access signal data
    std::string signal_group_str;
    std::string signal_dataset_str;
    if (single_read) {
        signal_group_str = "/Raw/Reads";
        signal_dataset_str = signal_group_str + "/" + read_name + "/Signal";
    } else {
        signal_group_str =  "/" + read_name + "/Raw";
        signal_dataset_str = signal_group_str + "/Signal";
    }

    std::cout << "Opening signal group " << signal_group_str << std::endl;
    H5::Group signal_group = f5.openGroup(signal_group_str);

    // Get the signal dataset
    //std::cout << "Opening dataset " << signal_dataset_str << std::endl;
    H5::DataSet signal_ds = f5.openDataSet(signal_dataset_str);

    //std::cout << "Getting data type..." << std::endl;
    H5::DataType mdatatype= signal_ds.getDataType();
    H5::DataSpace dataspace = signal_ds.getSpace();
    hsize_t dims[2];

    //std::cout << "Getting data dimensions..." << std::endl;
    dataspace.getSimpleExtentDims(dims, NULL); // rank = 1
    int data_count = (int) dims[0];
    std::cout << "Data count = " << data_count << std::endl;

    // Read the signals if available
    if (data_count == 0) {
        std::cout << "No signal data found." << std::endl;
        Base_Signals basecall_obj(read_info, sequence_data_str, basecall_signals);
        return basecall_obj;
    }

    // Store the signals in an array
    std::cout << "Reading data..." << std::endl;
    int16_t f5signals_16 [data_count];
    signal_ds.read(f5signals_16, mdatatype);

    // Cast signals to int
    //std::cout << "Creating int vector..." << std::endl;
    std::vector<int> f5signals (data_count);
    //std::cout << "Casting data..." << std::endl;
    for (int i = 0; i < data_count; i++) { f5signals[i] = (int) f5signals_16[i]; };
    //std::cout << "Data cast." << std::endl;

    // Get the sequence string if available
    std::cout << "Getting FASTQ..." << std::endl;
    std::string basecall_group;
    if (single_read) {
        basecall_group = "/Analyses/Basecall_1D_000";
    } else {
        basecall_group = "/" + read_name + "/Analyses/Basecall_1D_000";
    }
    std::vector<std::string> fq = getFastq(f5, basecall_group);
    if (fq.empty()){
        // Return the raw signal only
        std::cout << "No FASTQ found." << std::endl;
        basecall_signals.push_back(f5signals);

    } else {
        // Access signal basecalling information
        sequence_data_str = fq[1];

        // Get the block stride (window length) attribute
        std::cout << "Opening basecall group..." << std::endl;
        H5::Group group_obj;
        H5::Attribute block_stride_obj;
        if (single_read) {
            group_obj = f5.openGroup("/Analyses/Basecall_1D_000/Summary/basecall_1d_template");
            block_stride_obj = group_obj.openAttribute("block_stride");
        } else {
            group_obj = f5.openGroup("/" + read_name + "/Analyses/Basecall_1D_000/Summary/basecall_1d_template");
            block_stride_obj = group_obj.openAttribute("block_stride");
        }

        H5::DataType bs_datatype= block_stride_obj.getDataType();
        int bs_buffer [1];
        block_stride_obj.read(bs_datatype, bs_buffer);
        int block_stride_value = bs_buffer[0];

        // Get the sequence length attribute
        H5::Attribute seq_length_obj = group_obj.openAttribute("sequence_length");
        H5::DataType sl_datatype= seq_length_obj.getDataType();
        int sl_buffer [1];
        seq_length_obj.read(sl_datatype, sl_buffer);

        // Get the raw signal basecall start position
        //std::cout << "Opening segmentation group..." << std::endl;
        H5::Group segm_group_obj;
        if (single_read) {
            segm_group_obj = f5.openGroup("/Analyses/Segmentation_000/Summary/segmentation");
        } else {
            segm_group_obj = f5.openGroup("/" + read_name + "/Analyses/Segmentation_000/Summary/segmentation");
        }

        //std::cout << "Opening start attribute..." << std::endl;
        H5::Attribute start_attr_obj = segm_group_obj.openAttribute("first_sample_template");
        //std::cout << "Getting data type..." << std::endl;
        H5::DataType st_datatype= start_attr_obj.getDataType();
        int st_buffer [1];
        //std::cout << "Reading data..." << std::endl;
        start_attr_obj.read(st_datatype, st_buffer);
        //std::cout << "Data read." << std::endl;
        int start_index = st_buffer[0];

        // Get the boolean array of base calls (move)
        //std::cout << "Opening move dataset..." << std::endl;
        H5::DataSet move_dataset_obj;
        if (single_read) {
            std::string move_group = "/Analyses/Basecall_1D_000/BaseCalled_template/Move";
            std::cout << "Opening move dataset: " << move_group << std::endl;
            move_dataset_obj = f5.openDataSet(move_group);
        } else {
            std::string move_group = "/" + read_name + "/Analyses/Basecall_1D_000/BaseCalled_template/Move";
            std::cout << "Opening move dataset: " << move_group << std::endl;
            move_dataset_obj = f5.openDataSet(move_group);
        }

        //std::cout << "Getting data type..." << std::endl;
        H5::DataType  move_datatype = move_dataset_obj.getDataType();
        H5::DataSpace move_dataspace  = move_dataset_obj.getSpace();
        hsize_t move_dims[2];
        move_dataspace.getSimpleExtentDims(move_dims, NULL); // rank = 1
        int move_data_count = move_dims[0];

        // Get alignment information
        if (single_read) {
            map_positions = getMapPositions(f5, read_name);
            mapped_chrom = getChromosome(f5, read_name);
        } else {
            map_positions = getMapPositions(f5, read_name);
            mapped_chrom = getChromosome(f5, read_name);
        }

        // Read the boolean array
        uint8_t move_bool [move_data_count];
        move_dataset_obj.read(move_bool, move_datatype, move_dataspace);

        // Segment the raw signal by the window length
        std::vector<int> called_base_signal;  // Holds the current base's signal
        int basecall_index = 0;
        int base_start_index = start_index;
        int base_count = 0;
        for (int i = 0; i < move_data_count; i++)
        {
            // Grab the signal
            int end_index = base_start_index + block_stride_value;
//            std::vector<int> block_signal(block_stride_value);

            // Get the signal from base_start_index to end_index
//            block_signal = std::vector<int>(f5signals.begin() + base_start_index, f5signals.begin() + end_index);

            // Append the signal to the current base signal vector
//            called_base_signal.insert( called_base_signal.end(), block_signal.begin(), block_signal.end() );
            called_base_signal.insert( called_base_signal.end(), f5signals.begin() + base_start_index, f5signals.begin() + end_index );

            // Check whether a basecall occurred
            int move_value(move_bool[i]);
            if (move_value == 1)
            {
                // Update the base count
                base_count++;

                // Store the base's signal in the 2D array
                basecall_signals.push_back(called_base_signal);

                // Reset the current base signal vector
                called_base_signal.clear();
            }

            // Update indices
            base_start_index += block_stride_value;
            basecall_index ++;
        }
    }

    // Set up the read information header for the output report
    if (!mapped_chrom.empty()) {
        read_info += " " + mapped_chrom;
        if (!map_positions.empty()) {
            std::string map_start = std::to_string(map_positions[0]);
            read_info += " [" + map_start;
            std::string map_end = std::to_string(map_positions[1]);
            read_info += ", " + map_end + "]";
        }
    }

    // Append the basecall signals to the output structure
    Base_Signals basecall_obj(read_info, sequence_data_str, basecall_signals);

    return basecall_obj;
}

// Add base and read QC to the output data structure adn the output details file
static int writeBaseQCDetails(const char *input_file, Output_FAST5 &output_data, FILE *read_details_fp)
{
    int exit_code = 0;
    const char * read_name;
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

        // Check if it is a multi-read FAST5 file
        std::string signal_group_str;
        std::string read_name;
        std::string basecall_group = "/Analyses/Basecall_1D_000";
        bool single_read = false;
        if (f5.nameExists("/Raw")) {
            single_read = true;
            std::cout << "Single read FAST5" << std::endl;

            // Get the basecall group
            basecall_group = "/Analyses/Basecall_1D_000";

            // Get the sequence data
            std::vector<std::string> fq = getFastq(f5, basecall_group);

            // Append the basecalls to the output structure
            output_data.addReadFastq(fq, read_details_fp);

        } else {
            std::cout << "Multi-read FAST5" << std::endl;

            // Loop through each read
            std::cout << "Reading all reads" << std::endl;
            H5::Group root_group = f5.openGroup("/");
            size_t num_objs = root_group.getNumObjs();
            for (size_t i=0; i < num_objs; i++) {
                read_name = root_group.getObjnameByIdx(i);

                // First remove the prefix
                std::string read_id = read_name.substr(5);
                std::cout << "Processing read ID: " << read_id << std::endl;
                //std::cout << "Read: " << read_name << std::endl;

                // Get the basecall group
                basecall_group = "/" + read_name + "/Analyses/Basecall_1D_000";

                // Get the sequence data
                std::vector<std::string> fq = getFastq(f5, basecall_group);

                // Append the basecalls to the output structure
                output_data.addReadFastq(fq, read_details_fp);
            }
        }

    } catch (FileIException &error) {
        // If the dataset is missing, continue and ignore this file
        if (error.getFuncName() == "H5File::openDataSet") {
            std::cout << "No FASTQ sequence dataset found for file: " << input_file << std::endl;
        } else {
            std::cerr << "Error accessing the FAST5 file: " << input_file << std::endl;
            error.printErrorStack();
        }
        exit_code = 2;
    };

    return exit_code;
}


// Add read signal QC to the output data structure and the output details file
static int writeSignalQCDetails(const char *input_file, Output_FAST5 &output_data, std::vector<std::string> read_id_list)
{
    int exit_code = 0;

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

        // Check if it is a multi-read FAST5 file
        std::string signal_group_str;
        std::string read_name;
        if (f5.nameExists("/Raw")) {
            std::cout << "Single read FAST5" << std::endl;

            // Append the basecall signals to the output structure
            signal_group_str = "/Raw/Reads";
            read_name = getFileReadName(f5);
            std::cout << read_name << std::endl;
            Base_Signals basecall_obj = getReadBaseSignalData(f5, read_name, true);
            output_data.addReadBaseSignals(basecall_obj);
        } else {
            std::cout << "Multi-read FAST5" << std::endl;

            // Loop through each read
            std::cout << "Reading all reads" << std::endl;
            H5::Group root_group = f5.openGroup("/");
            size_t num_objs = root_group.getNumObjs();
            for (size_t i=0; i < num_objs; i++) {
                read_name = root_group.getObjnameByIdx(i);

                // Check if the read is in the list of reads to process (if provided)
                if (read_id_list.size() > 0) {
                    // First remove the prefix
                    std::string read_id = read_name.substr(5);
                    if (std::find(read_id_list.begin(), read_id_list.end(), read_id) == read_id_list.end()) {
                        //std::cout << "Skipping read ID: " << read_id << std::endl;
                        continue;
                    } else {
                        std::cout << "Processing read ID: " << read_id << std::endl;
                    }
                }

                std::cout << "Read: " << read_name << std::endl;

                // Get the basecall signals
                std::cout << "Getting basecall signals" << std::endl;
                Base_Signals basecall_obj = getReadBaseSignalData(f5, read_name, false);

                std::cout << "Adding basecall signals" << std::endl;
                output_data.addReadBaseSignals(basecall_obj);
            }
        }

    // catch failure caused by the H5File operations
    }
    catch (FileIException &error) {
        std::cerr << "FileIException" << std::endl;
        error.printErrorStack();
        exit_code = 2;
    }

    // catch failure caused by the DataSet operations
    catch (DataSetIException &error) {
        std::cerr << "DataSetIException" << std::endl;
        error.printErrorStack();
        exit_code = 2;
    }

    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException &error) {
        std::cerr << "DataSpaceIException" << std::endl;
        error.printErrorStack();
        exit_code = 2;
    }

    // catch failure caused by the Attribute operations
    catch (AttributeIException &error) {
        std::cerr << "AttributeIException" << std::endl;
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

    // Get the list of read IDs to process (if specified)
    std::vector<std::string> read_id_list;
    if (_input_data.read_ids.empty() == false) {
        // Parse the comma-separated list of read IDs
        std::stringstream ss(_input_data.read_ids);
        std::string read_id;
        while (std::getline(ss, read_id, ',')) {
            read_id_list.push_back(read_id);
        }
    }

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
            exit_code = writeSignalQCDetails(input_file, output_data, read_id_list);
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
        //output_data.seq_quality_info.base_quality_distribution.clear();
        output_data.seq_quality_info.read_average_base_quality_distribution.clear();

        output_data.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
        // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

        output_data.long_read_info.read_gc_content_count.resize(101, 0);
        // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

        output_data.long_read_info.NXX_read_length.resize(101, 0);
        // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;

        //output_data.seq_quality_info.base_quality_distribution.resize(256, 0);
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

                // Calculate read length statistics if base counts are not zero
                uint64_t total_num_bases = output_data.long_read_info.total_num_bases;
                if (total_num_bases == 0) {
                    std::cerr << "No bases found in input files." << std::endl;
                    exit_code = 3;
                } else {
                    // Calculate GC-content
                    output_data.long_read_info.gc_cnt = g_c / a_tu_g_c;

                    // Sort the read lengths in descending order
                    std::vector<int> read_lengths = output_data.long_read_info.read_lengths;
                    std::sort(read_lengths.begin(), read_lengths.end(), std::greater<int64_t>());

                    // Get the max read length
                    int64_t max_read_length = read_lengths.at(0);
                    output_data.long_read_info.longest_read_length = max_read_length;

                    // Get the median read length
                    int64_t median_read_length = read_lengths[read_lengths.size() / 2];
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
                            current_base_count += read_lengths.at(current_read_index);
                        }
                        int nxx_read_length = read_lengths.at(current_read_index);
                        output_data.long_read_info.NXX_read_length[percent_value] = nxx_read_length;
                    }

                    // Set common score variables
                    output_data.long_read_info.n50_read_length = output_data.long_read_info.NXX_read_length[50];
                    output_data.long_read_info.n95_read_length = output_data.long_read_info.NXX_read_length[95];
                    output_data.long_read_info.n05_read_length = output_data.long_read_info.NXX_read_length[5];

                    // Create the summary file
                    std::cout << "Writing summary file: " << read_summary_file.c_str() << std::endl;
                    read_summary_fp = fopen(read_summary_file.c_str(), "w");
                    fprintf(read_summary_fp, "total number of reads\t%d\n", output_data.long_read_info.total_num_reads);
                    fprintf(read_summary_fp, "total number of bases\t%ld\n", output_data.long_read_info.total_num_bases);
                    fprintf(read_summary_fp, "longest read length\t%d\n", output_data.long_read_info.longest_read_length);
                    fprintf(read_summary_fp, "N50 read length\t%d\n", output_data.long_read_info.n50_read_length);
                    fprintf(read_summary_fp, "mean read length\t%.2f\n", output_data.long_read_info.mean_read_length);
                    fprintf(read_summary_fp, "median read length\t%d\n", output_data.long_read_info.median_read_length);
                    fprintf(read_summary_fp, "GC%%\t%.2f\n", output_data.long_read_info.gc_cnt * 100);
                    fprintf(read_summary_fp, "\n\n");
                    for (int percent = 5; percent < 100; percent += 5)
                    {
                        fprintf(read_summary_fp, "N%02d read length\t%.d\n", percent, output_data.long_read_info.NXX_read_length[percent]);
                    }

                    fprintf(read_summary_fp, "\n\n");

                    fprintf(read_summary_fp, "GC content\tnumber of reads\n");
                    for (int gc_ratio = 0; gc_ratio < 100; gc_ratio++)
                    {
                        fprintf(read_summary_fp, "GC=%d%%\t%d\n", gc_ratio, output_data.long_read_info.read_gc_content_count[gc_ratio]);
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
                        fprintf(read_summary_fp, "%d\t%d\n", baseq, output_data.seq_quality_info.read_average_base_quality_distribution[baseq]);
                    }
                    fclose(read_summary_fp);
                }
            }
        }
    }

    return exit_code;
}
