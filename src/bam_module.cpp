/*
BAM_module.cpp:
Class for generating BAM file statistics. Records are accessed using multi-threading.
*/

#include <iostream>
#include <thread>
#include <iostream>
#include <cmath>

#include "bam_module.h"
#include "ComFunction.h"


int BAM_Module::calculateStatistics(Input_Para& input_params, Output_BAM& final_output){
    int exit_code = 0;
    auto relapse_start_time = std::chrono::high_resolution_clock::now();

    // Create mutexes
    std::mutex bam_mutex;
    std::mutex output_mutex;
    std::mutex cout_mutex;

    // Loop through the input files
    int file_count = (int) input_params.num_input_files;
    //std::cout << "Number of input files = " << file_count << std::endl;
    for (int i=0; i < file_count; i++){
        this->file_index = i;

        // Create a BAM reader
        std::string filepath(input_params.input_files[this->file_index]);
        HTSReader reader(filepath);
        std::cout<<"Processing file: "<< filepath << std::endl;

        // Get the number of records in the file using the BAM index
        int num_records = reader.getNumRecords(filepath);
        std::cout << "Number of records = " << num_records << std::endl;

        // Determine the number of records per thread
        int thread_count = input_params.threads;
        int batch_size = (int) ceil((double)num_records / (double)thread_count);
        std::cout << "Batch size (records per thread) = " << batch_size << std::endl;

         // Calculate statistics in batches
         while (reader.hasNextRecord()){
            std::cout << "Generating " << thread_count << " thread(s)..." << std::endl;
            std::vector<std::thread> thread_vector;
            for (int thread_index=0; thread_index<thread_count; thread_index++){
                //cout_mutex.lock();
                //::cout<<"Generated thread "<< thread_index+1 <<std::endl;
                //cout_mutex.unlock();

                // Create a thread
                std::thread t((BAM_Module::batchStatistics), std::ref(reader), batch_size, std::ref(input_params),std::ref(final_output), std::ref(bam_mutex), std::ref(output_mutex), std::ref(cout_mutex));
                thread_vector.push_back(std::move(t));
            }

            // Join the threads in thread_vector
            std::cout<<"Joining threads..."<<std::endl;
            int thread_index = 0;
            for (auto& t : thread_vector){
                t.join();
                cout_mutex.lock();
                //std::cout<<"Joined thread "<< thread_index+1 << std::endl;
                cout_mutex.unlock();
                thread_index++;
            }
            std::cout << "All threads joined." << std::endl;
        }
    }

    // Calculate the global sums across all records
    std::cout << "Calculating summary QC..." << std::endl;
    final_output.global_sum();
    std::cout << "QC complete" << std::endl;

    // Save the summary statistics to a file
    std::cout << "Saving summary statistics to file..." << std::endl;
    std::string summary_filepath = input_params.output_folder + "/bam_summary.txt";
    final_output.save_summary(summary_filepath, input_params, final_output);
    std::cout << "Saved file: " << summary_filepath << std::endl;

    auto relapse_end_time = std::chrono::high_resolution_clock::now();
    std::cout<<"Elapsed time (seconds) = "<< std::chrono::duration_cast<std::chrono::seconds>(relapse_end_time - relapse_start_time).count() << std::endl;

    return exit_code;
}

void BAM_Module::batchStatistics(HTSReader& reader, int batch_size, Input_Para& input_params, Output_BAM& final_output, std::mutex& bam_mutex, std::mutex& output_mutex, std::mutex& cout_mutex)
{
    // Read the next N records
    Output_BAM record_output;  // Output for the current batch
    reader.readNextRecords(batch_size, record_output, bam_mutex);

    // Update the final output
    output_mutex.lock();
    final_output.add(record_output);
    output_mutex.unlock();
//    // Print if records were processed
//    if (record_output.num_primary_alignment > 0){
//        cout_mutex.lock();
//        std::cout << "Thread " << std::this_thread::get_id() << " processed " << record_output.num_primary_alignment << " records" << std::endl;
//        cout_mutex.unlock();
//    }
}
