/*
BAM_module.cpp:
Class for generating BAM file statistics. Records are accessed using multi-threading.
*/

#include "bam_module.h"

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <iostream>
#include <cmath>
#include <unordered_set>
#include <algorithm>

#include "utils.h"
#include "ref_query.h"  // For reference genome analysis
#include "tin.h"  // For TIN analysis

// Run the BAM module
int BAM_Module::run(Input_Para &input_params, Output_BAM &final_output)
{
    // Determine if RRMS read ID filtering is required
    if (input_params.rrms_csv != ""){
        std::cout << "RRMS read ID filtering enabled" << std::endl;
        std::cout << "RRMS CSV file: " << input_params.rrms_csv << std::endl;

        // Determine if RRMS stats should be generated for accepted or rejected
        // reads
        if (input_params.rrms_filter){
            std::cout << "RRMS stats will be generated for accepted reads" << std::endl;
        } else {
            std::cout << "RRMS stats will be generated for rejected reads" << std::endl;
        }

        // Read the RRMS CSV file and store the read IDs
        std::cout << "Reading RRMS CSV file..." << std::endl;
        std::unordered_set<std::string> rrms_read_ids = readRRMSFile(input_params.rrms_csv, input_params.rrms_filter);
        std::cout << "Number of read IDs = " << rrms_read_ids.size() << std::endl;

        // Store the read IDs in the input parameters
        input_params.rrms_read_ids = rrms_read_ids;
    }

    // Calculate statistics
    int exit_code = calculateStatistics(input_params, final_output);

    return exit_code;
}

int BAM_Module::calculateStatistics(Input_Para &input_params, Output_BAM &final_output)
{
    int exit_code = 0;
    auto relapse_start_time = std::chrono::high_resolution_clock::now();

    // Create mutexes
    std::mutex bam_mutex;
    std::mutex output_mutex;
    std::mutex cout_mutex;

    // Get the modified base flag and threshold
    bool mod_analysis = input_params.mod_analysis;
    double base_mod_threshold = input_params.base_mod_threshold;

    // Get the gene BED file if available
    std::string gene_bed = input_params.gene_bed;
    if (gene_bed != ""){
        std::cout << "Gene BED file: " << gene_bed << std::endl;

        // Calculate TIN scores if the gene BED file is available
        std::cout << "Calculating TIN scores..." << std::endl;
        int sample_size = input_params.tin_sample_size;
        int min_cov = input_params.tin_min_coverage;
        calculateTIN(gene_bed, input_params.input_files[0], min_cov, sample_size, input_params.output_folder);
        std::cout << "TIN scores calculated." << std::endl;
    }

    // Loop through the input files
    int file_count = (int) input_params.num_input_files;
    for (int i = 0; i < file_count; i++) {
        this->file_index = i;

        // Create a BAM reader
        std::string filepath(input_params.input_files[this->file_index]);
        HTSReader reader(filepath);
        std::cout << "Processing file: "<< filepath << std::endl;

        // Get the number of threads (Set to 1 if the number of threads is not specified/invalid)
        int thread_count = input_params.threads;
        if (thread_count < 1) {
            thread_count = 1;
        }

        // Get the number of records in the file using the BAM index, and
        // process base modifications and TINs if available.
        // Note: This section utilizes one thread.
        std::cout << "Getting number of records..." << std::endl;
        int num_records = reader.getNumRecords(filepath, final_output, mod_analysis, base_mod_threshold);
        std::cout << "Number of records = " << num_records << std::endl;

        // Exit if there are no records
        if (num_records == 0){
            std::cerr << "No records found in file: " << filepath << std::endl;
            exit_code = 1;
            return exit_code;
        }

        // Determine the batch sizes if the user-specified thread count is greater than 1
        int batch_size = 0;
        if (thread_count > 1) {
            // If the number of records is less than the number of threads, set the number of threads to the number of records
            if (num_records < thread_count){
                thread_count = num_records;
                batch_size = 1;
                std::cout << "Number of threads set to " << thread_count << std::endl;
                std::cout << "Batch size (records per thread) = " << batch_size << std::endl;

            // Otherwise, the batch size is the number of records divided by the number of threads
            } else {
                // Determine the number of records per thread
                batch_size = (int) ceil((double)num_records / (double)thread_count);
                std::cout << "Batch size (records per thread) = " << batch_size << std::endl;
            }

        } else {
            // Set the batch size to the number of records
            batch_size = num_records;
        }

         // Calculate statistics in batches
         while (reader.hasNextRecord()){
            std::cout << "Generating " << thread_count << " thread(s)..." << std::endl;
            std::vector<std::thread> thread_vector;
            for (int thread_index=0; thread_index<thread_count; thread_index++){

                // Copy the input read IDs to a new vector
                std::unordered_set<std::string> rrms_read_ids_copy = input_params.rrms_read_ids;

                // Create a thread
                std::thread t((BAM_Module::batchStatistics), std::ref(reader), batch_size, rrms_read_ids_copy,std::ref(final_output), std::ref(bam_mutex), std::ref(output_mutex), std::ref(cout_mutex), base_mod_threshold);

                // Add the thread to the vector
                thread_vector.push_back(std::move(t));
            }

            // Join the threads in thread_vector
            std::cout<<"Joining threads..."<<std::endl;
            int thread_index = 0;
            for (auto& t : thread_vector){
                // Join the thread if it is joinable
                if (t.joinable()){
                    t.join();
                }
                thread_index++;
            }
            std::cout << "All threads joined." << std::endl;
        }
    }

    // Summarize the base modifications
    if (final_output.modified_prediction_count > 0) {

        // Determine CpG modification rate
        // First, read the reference genome
        if (input_params.ref_genome != ""){
            std::cout << "Reading reference genome for CpG modification rate calculation..." << std::endl;
            RefQuery ref_query;
            ref_query.setFilepath(input_params.ref_genome);
            std::cout << "Reference genome read." << std::endl;

            // Preprint revision: Loop through the base modifications and find
            // the CpG modifications
            for (auto const &it : final_output.sample_c_modified_positions)
            {
                std::string chr = it.first;
                std::vector<std::pair<int32_t, int>> cpg_mods = it.second;

                // Determine if it resides in a CpG site
                for (auto const &it2 : cpg_mods)
                {
                    int64_t ref_pos = (int64_t) it2.first;
                    int strand = it2.second;

                    // Update the strand-specific CpG modified base count
                    if (strand == 0 && ref_query.getBase(chr, ref_pos) == "C" && ref_query.getBase(chr, ref_pos + 1) == "G")
                    {
                        final_output.sample_cpg_forward_count++;
                    }
                    else if (strand == 1 && ref_query.getBase(chr, ref_pos) == "G" && ref_query.getBase(chr, ref_pos - 1) == "C")
                    {
                        final_output.sample_cpg_reverse_count++;
                    }
                }
            }
        }
    }

    // Calculate the global sums across all records
    std::cout << "Calculating summary QC..." << std::endl;
    final_output.global_sum();
    std::cout << "QC complete" << std::endl;

    // Save the summary statistics to a file
    std::cout << "Saving summary statistics to file..." << std::endl;

    // If in RRMS mode, append RRMS accepted/rejected to the output prefix
    std::string output_prefix = "bam";
    if (input_params.rrms_csv != ""){
        output_prefix += input_params.rrms_filter ? "_rrms_accepted" : "_rrms_rejected";
    } 
    std::string summary_filepath = input_params.output_folder + "/" + output_prefix + "_summary.txt";
    final_output.save_summary(summary_filepath, input_params, final_output);
    std::cout << "Saved file: " << summary_filepath << std::endl;

    auto relapse_end_time = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time (seconds) = " << std::chrono::duration_cast<std::chrono::seconds>(relapse_end_time - relapse_start_time).count() << std::endl;

    return exit_code;
}

void BAM_Module::batchStatistics(HTSReader& reader, int batch_size, std::unordered_set<std::string> read_ids, Output_BAM& final_output, std::mutex& bam_mutex, std::mutex& output_mutex, std::mutex& cout_mutex, double base_mod_threshold)
{
    // Read the next N records
    Output_BAM record_output;
    reader.readNextRecords(batch_size, record_output, bam_mutex, read_ids, base_mod_threshold);

    // Update the final output
    output_mutex.lock();
    final_output.add(record_output);
    output_mutex.unlock();
}

std::unordered_set<std::string> BAM_Module::readRRMSFile(std::string rrms_csv_file, bool accepted_reads)
{
    // Create an unordered set to store the read IDs for fast lookup
    std::unordered_set<std::string> rrms_read_ids;

    // Open the file
    std::ifstream rrms_file(rrms_csv_file);

    // Read the header and find the 'read_id' and 'decision' columns
    std::string header;
    std::vector<std::string> header_fields;
    std::getline(rrms_file, header);
    std::stringstream ss(header);
    std::string field;
    // std::cout << "RRMS CSV header:" << std::endl;
    while (std::getline(ss, field, ',')){
        header_fields.push_back(field);
        // std::cout << field << std::endl;
    }
    
    // Find the 'read_id' and 'decision' columns
    int read_id_index = -1;
    int decision_index = -1;
    for (size_t i=0; i<header_fields.size(); i++){
        if (header_fields[i] == "read_id"){
            read_id_index = i;
        } else if (header_fields[i] == "decision"){
            decision_index = i;
        }
    }

    // Exit if the read_id or decision columns are not found
    if (read_id_index == -1){
        std::cerr << "Error: 'read_id' column not found in RRMS CSV file" << std::endl;
        exit(1);
    }

    if (decision_index == -1){
        std::cerr << "Error: 'decision' column not found in RRMS CSV file" << std::endl;
        exit(1);
    }

    // Read all rows in the file and store the read IDs if the decision is
    // 'stop_receiving' for accepted, or 'unblock' for rejected reads.
    std::string pattern = accepted_reads ? "stop_receiving" : "unblock";
    std::string line;
    while (std::getline(rrms_file, line)){
        std::vector<std::string> fields;
        std::string field;
        std::stringstream ss(line);
        while (std::getline(ss, field, ',')){
            fields.push_back(field);
        }

        // Get the read ID and decision
        std::string read_id = fields[read_id_index];
        std::string decision = fields[decision_index];

        // Store the read ID if the decision matches the pattern
        if (decision == pattern){
            rrms_read_ids.insert(read_id);

            std::string test_id = "65d8befa-eec0-4496-bf2b-aa1a84e6dc5e";
            if (read_id == test_id){
                std::cout << "[TEST 1] Found test ID: " << test_id << std::endl;
            }

            // std::cout << read_id << std::endl;
        }
    }
    
    // Close the file
    rrms_file.close();

    return rrms_read_ids;
}