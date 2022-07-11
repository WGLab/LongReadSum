/*
F5_module.cpp:
Class for calling FAST5 statistics modules.

*/

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "FAST5_module.h"
#include "ComFunction.h"
#include "H5Cpp.h"
using namespace H5;

const H5std_string FILE_NAME("/home/perdomoj/github/LongReadSum/output/SimpleCompound.h5");
const H5std_string DATASET_NAME("dset");
const int          NX   = 4; // dataset dimensions
const int          NY   = 6;
const int          RANK = 2;


size_t FAST5_Module::batch_size_of_record=3000;
std::mutex FAST5_Module::myMutex_readFAST5;
std::mutex FAST5_Module::myMutex_output;
size_t FAST5_Module::file_index = 0;

FAST5_Thread_data::FAST5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   _batch_size = p_batch_size;
   _thread_id = p_thread_id;
   _input_parameters = ref_input_op;
   stored_records.reserve(p_batch_size+1);
   for (int _i_=0; _i_<p_batch_size+1; _i_++){
      stored_records.push_back( FAST5_Record() );
   }
}

FAST5_Thread_data::~FAST5_Thread_data(){
}

size_t FAST5_Thread_data::readRecord(std::ifstream* file_stream){
    thread_index = 0;  // Index where this thread's data will be stored
//    while( std::getline( *file_stream, current_line )) {
//        std::istringstream column_stream( current_line );
//
//        // Read each column value from the record line
//        std::string column_value;
//        std::vector<std::string> column_values;
//    }

    return thread_index;
}

// FAST5_module Constructor
FAST5_Module::FAST5_Module(Input_Para& input_parameters){
    _input_parameters = input_parameters;
    has_error = 0;
    file_index = 0;

    if (file_index >= _input_parameters.num_input_files){
        std::cerr << "Input file list error." << std::endl;
        has_error |= 1;
        return;
    }

    // Check if the first file exists
    const char * first_filepath = _input_parameters.input_files[file_index].c_str();
    std::ifstream *input_file_stream = new std::ifstream(first_filepath);
    if (!(input_file_stream->is_open())){
        std::cerr<< "Cannot open FAST5 file: "<< first_filepath <<std::endl;
        has_error |= 2;
    }else{
        // Close the file
        input_file_stream->close();
        input_file_stream->clear();

        // Set the current active filepath
        input_filepath = first_filepath;
        file_index ++;
        std::cout<< "Opened FAST5 file: "<< first_filepath <<" " << file_index<<"/"<<_input_parameters.num_input_files <<std::endl;
    }
}

FAST5_Module::~FAST5_Module(){
//   if (input_file_stream!=NULL){
//       delete input_file_stream;
//   }
//   input_file_stream = NULL;
}

int FAST5_Module::generateStatistics( Output_FAST5& output_data){
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   output_data.all_long_read_info.long_read_info.resize();
   output_data.passed_long_read_info.long_read_info.resize();
   output_data.failed_long_read_info.long_read_info.resize();

   if (has_error==0){
       all_threads.reserve(_input_parameters.threads+3);

      int _i_t=0;
      FAST5_Thread_data** thread_data_vector = new FAST5_Thread_data*[_input_parameters.threads];
      try{
         for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
             std::cout<<"INFO: generate threads "<<_i_t<<std::endl<<std::flush;
             thread_data_vector[_i_t] = new FAST5_Thread_data(_input_parameters, _i_t, FAST5_Module::batch_size_of_record);
             std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
             all_threads.push_back(std::thread((FAST5_Module::createThread), input_filepath, std::ref(_input_parameters), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(output_data) ));
         }

         std::cout<<"INFO: join threads"<<std::endl<<std::flush;
         for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
             std::cout<<"INFO: join threads "<<_i_t<<std::endl<<std::flush;
             all_threads[_i_t].join();
         }
      }catch(const std::runtime_error& re){
         std::cerr << "Runtime error: " << re.what() << std::endl;
      }catch(const std::exception& ex){
         std::cerr << "Error occurred: " << ex.what() << std::endl;
      }catch(...){
         std::cerr << "Unknown failure occurred. Possible memory corruption" << std::endl;
      }
     
      for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
         delete thread_data_vector[_i_t];
      }
      delete [] thread_data_vector;
   }

   output_data.global_sum();
 
   auto relapse_end_time = std::chrono::high_resolution_clock::now();
   std::cout<<"Total time(Elapsed): "<<round3((relapse_end_time-relapse_start_time).count()/1000000000.0)<<std::endl;

   std::cout<<"FAST5 Module "<< (has_error==0?"completed successfully":"failed") << std::endl;
 
   return has_error;
}

void FAST5_Module::createThread(std::string input_filepath, Input_Para& ref_input_op, int thread_id, FAST5_Thread_data& ref_thread_data, Output_FAST5& ref_output ){
    size_t read_ss_size, read_ss_i;
//    while (true){
//        myMutex_readFAST5.lock();
//        read_ss_size = ref_thread_data.readRecord(file_stream);
//
//        if (read_ss_size == 0 && !(file_index < ref_input_op.num_input_files) ){
//            myMutex_readFAST5.unlock();
//            break;
//        }
//        if ( read_ss_size < batch_size_of_record ){
//            if ( file_index < ref_input_op.num_input_files ){
//               std::cout<< "INFO: Open FAST5 file="<< ref_input_op.input_files[file_index] <<std::endl;
//               file_stream->close();
//               file_stream->clear();
//
//               file_stream->open( ref_input_op.input_files[file_index].c_str() );
//               std::string firstline;
//               std::getline( *file_stream, firstline );
//               file_index++;
//            }
//        }
//        myMutex_readFAST5.unlock();
//        if (read_ss_size == 0 ) { continue; }
//
//        // Columns used for statistics: passes_filtering, sequence_length_template, mean_qscore_template
//        ref_thread_data.output_data.reset();
//        ref_thread_data.output_data.all_long_read_info.long_read_info.resize();
//        ref_thread_data.output_data.passed_long_read_info.long_read_info.resize();
//        ref_thread_data.output_data.failed_long_read_info.long_read_info.resize();
//        for(read_ss_i=0; read_ss_i<read_ss_size; read_ss_i++){
//           FAST5_Statistics* output_statistics = NULL;
//           bool passes_filtering_value = ref_thread_data.stored_records[read_ss_i].passes_filtering;
//           if ( passes_filtering_value == true) {
//               output_statistics = &(ref_thread_data.output_data.passed_long_read_info);
//           } else {
//                output_statistics = &(ref_thread_data.output_data.failed_long_read_info);
//           }
//           output_statistics->long_read_info.total_num_reads++;
//           size_t sequence_base_count = ref_thread_data.stored_records[read_ss_i].sequence_length_template;
//           output_statistics->long_read_info.total_num_bases += sequence_base_count;
//
//           if ( output_statistics->long_read_info.longest_read_length < ref_thread_data.stored_records[read_ss_i].sequence_length_template){
//               output_statistics->long_read_info.longest_read_length = ref_thread_data.stored_records[read_ss_i].sequence_length_template;
//           }
//           output_statistics->long_read_info.read_length_count[ ref_thread_data.stored_records[read_ss_i].sequence_length_template<MAX_READ_LENGTH?ref_thread_data.stored_records[read_ss_i].sequence_length_template:(MAX_READ_LENGTH-1) ] += 1;
//
//           output_statistics->seq_quality_info.read_quality_distribution[ int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template ) ] += 1;
//           if ( output_statistics->seq_quality_info.min_read_quality == MoneDefault ||
//               output_statistics->seq_quality_info.min_read_quality>int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template ) ){
//              output_statistics->seq_quality_info.min_read_quality = int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template );
//           }
//           if ( output_statistics->seq_quality_info.max_read_quality < int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template) ){
//              output_statistics->seq_quality_info.max_read_quality = int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template);
//           }
//        }
//
//        myMutex_output.lock();
//        ref_output.add( ref_thread_data.output_data );
//        myMutex_output.unlock();
//    }
}
