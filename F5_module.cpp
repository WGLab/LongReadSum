/*
F5_module.cpp:
Class for calling FAST5 statistics modules.

*/

#include <iostream>
#include <string>
#include <algorithm>

#include "F5_module.h"
#include "ComFunction.h"


size_t F5_Module::batch_size_of_record=3000;
std::mutex F5_Module::myMutex_readF5;
std::mutex F5_Module::myMutex_output;
size_t F5_Module::file_index = 0;

F5_Thread_data::F5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   _batch_size = p_batch_size;
   _thread_id = p_thread_id;
   _input_parameters = ref_input_op;
   stored_records.reserve(p_batch_size+1);
   for (int _i_=0; _i_<p_batch_size+1; _i_++){
      stored_records.push_back( F5_SS_Record() );
   }
}

F5_Thread_data::~F5_Thread_data(){ 
}

size_t F5_Thread_data::read_ss_record(std::ifstream* file_stream){
    std::cout << "Type 1." << std::endl;
    thread_index = 0;  // Index where this thread's data will be stored
    while( std::getline( *file_stream, current_line )) {
        std::istringstream current_line_str( current_line );

        // Read variables from the record line
        F5_SS_Record record_data;
        current_line_str >> record_data.filename;
        std::cout << "Stored: " << record_data.filename << std::endl;
//        try {
//            current_line_str >> record_data.filename;
//            current_line_str >> record_data.read_id;
//            current_line_str >> record_data.run_id;
//            current_line_str >> record_data.channel;
//            current_line_str >> record_data.start_time;
//            current_line_str >> record_data.duration;
//            current_line_str >> record_data.num_events;
//            current_line_str >> record_data.passes_filtering;
//            current_line_str >> record_data.template_start;
//            current_line_str >> record_data.num_events_template;
//            current_line_str >> record_data.template_duration;
//            current_line_str >> record_data.num_called_template;
//            current_line_str >> record_data.sequence_length_template;
//            current_line_str >> record_data.mean_qscore_template;
//            current_line_str >> record_data.strand_score_template;
//            current_line_str >> record_data.calibration_strand_genome_template;
//            current_line_str >> record_data.calibration_strand_identity_template;
//            current_line_str >> record_data.calibration_strand_accuracy_template;
//            current_line_str >> record_data.calibration_strand_speed_bps_template;
//
//            // Store the record
//            stored_records[thread_index] = record_data;
//        } catch (const std::exception& e) {
//            std::cerr << e.what();
//        }



//       if (!(current_line_str >> record_data[thread_index].filename >> record_data[thread_index].read_id
//                 >> record_data[thread_index].run_id >> record_data[thread_index].channel
//                 >> record_data[thread_index].start_time >> record_data[thread_index].duration
//                 >> record_data[thread_index].num_events >> record_data[thread_index].passes_filtering
//                 >> record_data[thread_index].template_start >> record_data[thread_index].num_events_template
//                 >> record_data[thread_index].template_duration >> record_data[thread_index].num_called_template
//                 >> record_data[thread_index].sequence_length_template >> record_data[thread_index].mean_qscore_template >> record_data[thread_index].strand_score_template
//                 >> record_data[thread_index].calibration_strand_genome_template >> record_data[thread_index].calibration_strand_identity_template
//                 >> record_data[thread_index].calibration_strand_accuracy_template >> record_data[thread_index].calibration_strand_speed_bps_template)) {
//             std::cerr<<"Error! for <"<<current_line<<">"<<std::endl;
//             break;
//       }

       thread_index++;
       if ( thread_index >= _batch_size){ break; }
    }

    return thread_index;
}
size_t F5_Thread_data::read_ss_record_2(std::ifstream* file_stream){
    std::cout << "Type 2." << std::endl;
    thread_index = 0;
    std::string fq_file;
    int muxv;
    float median_template, mad_template;
    while( std::getline( *file_stream, current_line )) {
       std::istringstream current_line_str( current_line );
       F5_SS_Record _t_f5_ss_record;
       if (!(current_line_str >> fq_file >> stored_records[thread_index].filename >> stored_records[thread_index].read_id
                 >> stored_records[thread_index].run_id >> stored_records[thread_index].channel >> muxv
                 >> stored_records[thread_index].start_time >> stored_records[thread_index].duration
                 >> stored_records[thread_index].num_events >> stored_records[thread_index].passes_filtering
                 >> stored_records[thread_index].template_start >> stored_records[thread_index].num_events_template
                 >> stored_records[thread_index].template_duration
                 >> stored_records[thread_index].sequence_length_template
                 >> stored_records[thread_index].mean_qscore_template >> stored_records[thread_index].strand_score_template
                 >> median_template >> mad_template)) {
             // TODO: Make the warning below more informative, here it prints out the whole read record  w/o info about the line parsing error.
             std::cerr<<"Error! for <"<<current_line<<">"<<std::endl;
             break;
       }

       thread_index++;
       if ( thread_index >= _batch_size){ break; }
    }

    return thread_index;
}
size_t F5_Thread_data::read_ss_record_3(std::ifstream* file_stream){
    std::cout << "Type 3." << std::endl;
    thread_index = 0;
    std::string fq_file;
    int muxv;
    int batch_id;
    float median_template, mad_template;
    while( std::getline( *file_stream, current_line )) {
       std::istringstream current_line_str( current_line );
       F5_SS_Record _t_f5_ss_record;
       // TODO: How many of these stored variables are used for generating statistics?
       if (!(current_line_str >> fq_file >> stored_records[thread_index].filename >> stored_records[thread_index].read_id
                 >> stored_records[thread_index].run_id >> batch_id >> stored_records[thread_index].channel >> muxv
                 >> stored_records[thread_index].start_time >> stored_records[thread_index].duration
                 >> stored_records[thread_index].num_events >> stored_records[thread_index].passes_filtering
                 >> stored_records[thread_index].template_start >> stored_records[thread_index].num_events_template
                 >> stored_records[thread_index].template_duration
                 >> stored_records[thread_index].sequence_length_template
                 >> stored_records[thread_index].mean_qscore_template >> stored_records[thread_index].strand_score_template
                 >> median_template >> mad_template)) {
             std::cout<<"Error! for <"<<current_line<<">"<<std::endl;
             break;
       }

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

        // Determine the columns in this file (tab-delimited)
        std::string column_line;
        std::getline( *input_file_stream, column_line );
        std::istringstream column_stream(column_line);
        std::string column_name;
        std::vector<std::string> parsed_column_names;
        while (std::getline( column_stream, column_name, '\t' )) {
            parsed_column_names.push_back(column_name);
        }

        // Ensure that we found the columns we need for statistics: passes_filtering, sequence_length_template, mean_qscore_template
        if (requiredHeadersFound(parsed_column_names))
        {
            // Update the column names
            column_names = parsed_column_names;

            // Print the column names
            std::cout << "\n=====\nColumn names: " << std::endl;
            for (std::string i: column_names)
            {
                std::cout << i << std::endl;
            }
            std::cout << "=====\n" << std::endl;
        } else {
            has_error = 4;
        }
    }
}

bool F5_Module::requiredHeadersFound(std::vector<std::string> header_vector) {
    // Ensure that we have the headers we need for statistics: passes_filtering, sequence_length_template, mean_qscore_template
    bool headers_found = true;
    std::vector<std::string> required_headers = { "passes_filtering", "sequence_length_template", "mean_qscore_template" };
    for (std::string required_header: required_headers)
    {
        if (std::find(header_vector.begin(), header_vector.end(), required_header) == header_vector.end())
        {
          headers_found = false;
          std::cerr << "Required header '" << required_header << "' not found." << std::endl;
        }
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
             thread_data_vector[_i_t] = new F5_Thread_data(_input_parameters, _i_t, F5_Module::batch_size_of_record);
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

    // Determine the format of the FAST5 file
    int64_t format_type = ((ref_input_op.other_flags) & 15);
    while (true){
        myMutex_readF5.lock();

        std::cout << "Format type: " << std::to_string(format_type) << std::endl;
        if ( format_type ==1 ){
           read_ss_size = ref_thread_data.read_ss_record(file_stream);
        }else if ( format_type ==2 ){
           read_ss_size = ref_thread_data.read_ss_record_2(file_stream);
        }else if ( format_type ==3 ){
           read_ss_size = ref_thread_data.read_ss_record_3(file_stream);
        }else{
            std::cerr << "No records read." << std::endl;
            read_ss_size = 0;
        }

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
           std::transform( ref_thread_data.stored_records[read_ss_i].passes_filtering.begin(),
                           ref_thread_data.stored_records[read_ss_i].passes_filtering.end(),
                           ref_thread_data.stored_records[read_ss_i].passes_filtering.begin(), [](unsigned char c) -> unsigned char { return std::tolower(c); } );
           if ( ref_thread_data.stored_records[read_ss_i].passes_filtering.compare("true")==0){
               _f5_st = &(ref_thread_data.t_output_F5_.f5_passed_long_read_info);
           }else{
               _f5_st = &(ref_thread_data.t_output_F5_.f5_failed_long_read_info);
           }
           _f5_st->long_read_info.total_num_reads++;
           _f5_st->long_read_info.total_num_bases += ref_thread_data.stored_records[read_ss_i].sequence_length_template;
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
