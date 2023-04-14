/*
BAM_module.cpp:
Class for generating BAM file statistics. Records are accessed using multi-threading.
*/

#include <iostream>
#include <thread>
#include <iostream>
#include <cmath>

#include "BAM_module.h"
#include "ComFunction.h"

BAM_Thread_data::BAM_Thread_data(Input_Para& input_params, int p_thread_id, int p_batch_size=1000){
//   this->record_list.reserve(p_batch_size+1);
   _thread_id = p_thread_id;
   input_params = input_params;
}

BAM_Thread_data::~BAM_Thread_data(){
}


int BAM_Module::calculateStatistics(Input_Para& input_params, Output_BAM& final_output){
    int exit_code = 0;
    auto relapse_start_time = std::chrono::high_resolution_clock::now();

    // Create mutexes
    std::mutex bam_mutex;
    std::mutex output_mutex;
    std::mutex cout_mutex;

    // Loop through the input files
    int file_count = (int) input_params.num_input_files;
    std::cout << "INFO: Number of input files = " << file_count << std::endl;
    for (int i=0; i < file_count; i++){
        this->file_index = i;

        // Create a BAM reader
        std::string filepath(input_params.input_files[this->file_index]);
        HTSReader reader(filepath);
        std::cout<<"INFO: Processing file "<< filepath << std::endl;

        // Determine the number of records per thread
        int thread_count = (double)input_params.threads;
        int batch_size = 10;
        std::cout << "INFO: Number of records per thread = " << batch_size << std::endl;

         // Calculate statistics in batches
         while (reader.hasNextRecord()){
            std::vector<std::thread> thread_vector;
            for (int thread_index=0; thread_index<thread_count; thread_index++){
                std::cout<<"INFO: Generated thread "<< thread_index+1 <<std::endl;

                std::thread t((BAM_Module::batchStatistics), std::ref(reader), batch_size, std::ref(input_params),std::ref(final_output), std::ref(bam_mutex), std::ref(output_mutex), std::ref(cout_mutex));
                thread_vector.push_back(std::move(t));
            }
            std::cout << "file_index = " << i << std::endl;

            // Join the threads in thread_vector
            int thread_index = 0;
            for (auto& t : thread_vector){
                t.join();
                cout_mutex.lock();
                std::cout<<"INFO: Joined thread "<< thread_index+1 << std::endl;
                cout_mutex.unlock();
                thread_index++;
            }
        }
    }
    // Print the output data results
    std::cout << "INFO: Number of primary alignments = " << final_output.num_primary_alignment << std::endl;
    std::cout << "INFO: Number of mismatched bases = " << final_output.num_mismatched_bases << std::endl;
    return exit_code;
}

void BAM_Module::batchStatistics(HTSReader& reader, int batch_size, Input_Para& input_params, Output_BAM& final_output, std::mutex& bam_mutex, std::mutex& output_mutex, std::mutex& cout_mutex)
{
    // Lock the cout mutex
    cout_mutex.lock();
    std::cout << "INFO: Starting batch statistics" << std::endl;
    cout_mutex.unlock();

    // Read the next N records
    Output_BAM record_output;
    int exit_code = reader.readNextRecords(batch_size, record_output, bam_mutex);
    int record_count = record_output.num_primary_alignment;

    if (record_count > 0) {
        // Update the final output
        output_mutex.lock();
        final_output.num_primary_alignment += record_output.num_primary_alignment;
        final_output.num_mismatched_bases  += record_output.num_mismatched_bases;
        output_mutex.unlock();
        std::cout << "Processed " << record_count << " records" << std::endl;
    }
}
