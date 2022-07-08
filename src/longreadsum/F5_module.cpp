/*
F5_module.cpp:
Class for calling FAST5 statistics modules.

*/

#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

#include "F5_module.h"
#include "ComFunction.h"


size_t F5_Module::batch_size_of_record=3000;
std::mutex F5_Module::myMutex_readF5;
std::mutex F5_Module::myMutex_output;
size_t F5_Module::file_index = 0;

F5_Thread_data::F5_Thread_data(Input_Para& ref_input_op, std::map<std::string, int> header_columns, int p_thread_id, int p_batch_size=1000){
   _batch_size = p_batch_size;
   _thread_id = p_thread_id;
   _input_parameters = ref_input_op;
   stored_records.reserve(p_batch_size+1);
   for (int _i_=0; _i_<p_batch_size+1; _i_++){
      stored_records.push_back( F5_SS_Record() );
   }
   _header_columns = header_columns;
}

F5_Thread_data::~F5_Thread_data(){ 
}

std::map<std::string, int> F5_Thread_data::getHeaderColumns()
{
    return _header_columns;
}

size_t F5_Thread_data::read_ss_record(std::ifstream* file_stream, std::map<std::string, int> header_columns){
    std::cout << "Type 1." << std::endl;
    thread_index = 0;  // Index where this thread's data will be stored
    while( std::getline( *file_stream, current_line )) {
        std::istringstream column_stream( current_line );

        // Read each column value from the record line
        std::string column_value;
        std::vector<std::string> column_values;
        while (std::getline( column_stream, column_value, '\t' )) {
            column_values.push_back(column_value);
        }

        // Store the required values from the record
        F5_SS_Record record_data;

        // passes_filtering
        std::string passes_filtering_str = column_values[header_columns.at("passes_filtering")];

        // Make lowercase
        std::transform(passes_filtering_str.begin(), passes_filtering_str.end(), passes_filtering_str.begin(), [](unsigned char c){ return std::tolower(c); });

        // Convert to boolean
        bool passes_filtering;
        std::istringstream(passes_filtering_str) >>  std::boolalpha >> passes_filtering;
        record_data.passes_filtering = passes_filtering;

        // sequence_length_template
        std::stringstream sstream(column_values[header_columns.at("sequence_length_template")]);
        size_t sequence_length_template;
        sstream >> sequence_length_template;
        record_data.sequence_length_template = sequence_length_template;

        // mean_qscore_template
        float mean_qscore_template;
        mean_qscore_template = std::stof(column_values[header_columns.at("mean_qscore_template")]);
        record_data.mean_qscore_template = mean_qscore_template;

        // Store the records
        stored_records[thread_index].passes_filtering = passes_filtering;
        stored_records[thread_index].sequence_length_template = sequence_length_template;
        stored_records[thread_index].mean_qscore_template = mean_qscore_template;

        // Update the thread index
        thread_index++;
        if ( thread_index >= _batch_size){ break; }
    }

    return thread_index;
}

// F5_Module Constructor
F5_Module::F5_Module(Input_Para& input_parameters){
    _input_parameters = input_parameters;
    has_error = 0;
    file_index = 0;

    input_file_stream = NULL;
    if (file_index >= _input_parameters.num_input_files){
        std::cerr << "Input file list error." << std::endl;
        has_error |= 1;
        return;
    }

    // Open the first file in the list
    const char * first_filepath = _input_parameters.input_files[file_index].c_str();
    input_file_stream = new std::ifstream(first_filepath);
    if (!(input_file_stream->is_open())){
        std::cout<< "Error!!! Cannot open FAST5 file="<< first_filepath <<std::endl;
        has_error |= 2;
    }else{
        file_index ++;
        std::cout<< "INFO: Open FAST5 file = "<< first_filepath <<" " << file_index<<"/"<<_input_parameters.num_input_files <<std::endl;

        // Ensure that we have the columns we need for statistics
        std::string column_line;
        std::getline( *input_file_stream, column_line );
        if (requiredHeadersFound(column_line))
        {
//            // Print the column names
//            std::cout << "\n=====\nHeader names: " << std::endl;
//            for (std::string i: header_names)
//            {
//                std::cout << i << std::endl;
//            }
//            std::cout << "=====\n" << std::endl;
        } else {
            has_error = 4;
        }
    }
}

bool F5_Module::requiredHeadersFound(std::string header_string) {
    // Ensure that we have the headers we need for statistics: passes_filtering, sequence_length_template, mean_qscore_template
    std::vector<std::string> required_headers = { "passes_filtering", "sequence_length_template", "mean_qscore_template" };
    std::istringstream header_stream(header_string);
    std::string header_name;
    int current_column = 0;
    while (std::getline( header_stream, header_name, '\t' )) {
        header_names.push_back(header_name);

        // Check if this column contains a required header
        std::vector<std::string>::iterator header_it;
        for (header_it = required_headers.begin(); header_it != required_headers.end();)
        {
            if (*header_it == header_name)
            {
                // Store the header column for later use in parsing records
                _header_columns[header_name] = current_column;

                // Remove the header from the search vector
                required_headers.erase(header_it);

            } else {
                ++header_it;
            }
        }
        current_column++;
    }

    // Store the total column count
    column_count = current_column;

    // Return true if all required headers were found
    bool headers_found = false;
    if (required_headers.empty()) {
        headers_found = true;
    } else {
        // Print the missing headers
        for (std::string missing_header: required_headers)
            std::cerr << "Required header '" << missing_header << "' not found." << std::endl;
    }
    return headers_found;
}

F5_Module::~F5_Module(){
   if (input_file_stream!=NULL){
       delete input_file_stream;
   }
   input_file_stream = NULL;
}

int F5_Module::F5_st( Output_F5& t_output_F5_info){  
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   t_output_F5_info.f5_long_read_info.long_read_info.resize();
   t_output_F5_info.f5_passed_long_read_info.long_read_info.resize();
   t_output_F5_info.f5_failed_long_read_info.long_read_info.resize();

   if (has_error==0){
       m_threads.reserve(_input_parameters.threads+3);

      int _i_t=0;
      F5_Thread_data** thread_data_vector = new F5_Thread_data*[_input_parameters.threads];
      try{
         for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
             std::cout<<"INFO: generate threads "<<_i_t<<std::endl<<std::flush;
             thread_data_vector[_i_t] = new F5_Thread_data(_input_parameters, _header_columns, _i_t, F5_Module::batch_size_of_record);
             std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
             m_threads.push_back(std::thread((F5_Module::F5_do_thread), input_file_stream, std::ref(_input_parameters), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(t_output_F5_info) ));
         }

         std::cout<<"INFO: join threads"<<std::endl<<std::flush;
         for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
             std::cout<<"INFO: join threads "<<_i_t<<std::endl<<std::flush;
             m_threads[_i_t].join();
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

   t_output_F5_info.global_sum();
 
   auto relapse_end_time = std::chrono::high_resolution_clock::now();
   std::cout<<"Total time(Elapsed): "<<round3((relapse_end_time-relapse_start_time).count()/1000000000.0)<<std::endl;

   std::cout<<"FAST5 Module "<< (has_error==0?"completed successfully":"failed") << std::endl;
 
   return has_error;
}

void F5_Module::F5_do_thread(std::ifstream* file_stream, Input_Para& ref_input_op, int thread_id, F5_Thread_data& ref_thread_data, Output_F5& ref_output ){
    size_t read_ss_size, read_ss_i;
    while (true){
        myMutex_readF5.lock();
        std::map<std::string, int> header_column_data = ref_thread_data.getHeaderColumns();
        read_ss_size = ref_thread_data.read_ss_record(file_stream, header_column_data);

        if (read_ss_size == 0 && !(file_index < ref_input_op.num_input_files) ){
            myMutex_readF5.unlock();
            break;
        }
        if ( read_ss_size < batch_size_of_record ){
            if ( file_index < ref_input_op.num_input_files ){ 
               std::cout<< "INFO: Open F5 file="<< ref_input_op.input_files[file_index] <<std::endl;
               file_stream->close();
               file_stream->clear();

               file_stream->open( ref_input_op.input_files[file_index].c_str() );
               std::string firstline;
               std::getline( *file_stream, firstline );
               file_index++;
            }
        }
        myMutex_readF5.unlock();
        if (read_ss_size == 0 ) { continue; }

        // Columns used for statistics: passes_filtering, sequence_length_template, mean_qscore_template
        ref_thread_data.t_output_F5_.reset();
        ref_thread_data.t_output_F5_.f5_long_read_info.long_read_info.resize();
        ref_thread_data.t_output_F5_.f5_passed_long_read_info.long_read_info.resize();
        ref_thread_data.t_output_F5_.f5_failed_long_read_info.long_read_info.resize();
        for(read_ss_i=0; read_ss_i<read_ss_size; read_ss_i++){
           Basic_F5_Statistics* _f5_st = NULL;
           bool passes_filtering_value = ref_thread_data.stored_records[read_ss_i].passes_filtering;
           if ( passes_filtering_value == true) {
               _f5_st = &(ref_thread_data.t_output_F5_.f5_passed_long_read_info);
           } else {
                _f5_st = &(ref_thread_data.t_output_F5_.f5_failed_long_read_info);
           }
           _f5_st->long_read_info.total_num_reads++;
           size_t sequence_base_count = ref_thread_data.stored_records[read_ss_i].sequence_length_template;
           _f5_st->long_read_info.total_num_bases += sequence_base_count;

           if ( _f5_st->long_read_info.longest_read_length < ref_thread_data.stored_records[read_ss_i].sequence_length_template){
               _f5_st->long_read_info.longest_read_length = ref_thread_data.stored_records[read_ss_i].sequence_length_template;
           }
           _f5_st->long_read_info.read_length_count[ ref_thread_data.stored_records[read_ss_i].sequence_length_template<MAX_READ_LENGTH?ref_thread_data.stored_records[read_ss_i].sequence_length_template:(MAX_READ_LENGTH-1) ] += 1;

           _f5_st->seq_quality_info.read_quality_distribution[ int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template ) ] += 1;
           if ( _f5_st->seq_quality_info.min_read_quality == MoneDefault ||
               _f5_st->seq_quality_info.min_read_quality>int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template ) ){
              _f5_st->seq_quality_info.min_read_quality = int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template );
           }
           if ( _f5_st->seq_quality_info.max_read_quality < int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template) ){
              _f5_st->seq_quality_info.max_read_quality = int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template);
           }
        }

        myMutex_output.lock();
        ref_output.add( ref_thread_data.t_output_F5_ );
        myMutex_output.unlock();
    }
}
