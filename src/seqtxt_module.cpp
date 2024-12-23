/*
SeqTxt_module.cpp:
Class for calling FAST5 statistics modules.

*/

#include "seqtxt_module.h"

/// @cond
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
/// @endcond

#include "utils.h"


size_t SeqTxt_Module::batch_size_of_record=3000;
std::mutex SeqTxt_Module::myMutex_readSeqTxt;
std::mutex SeqTxt_Module::myMutex_output;
size_t SeqTxt_Module::file_index = 0;

SeqTxt_Thread_data::SeqTxt_Thread_data(Input_Para& ref_input_op, std::map<std::string, int> header_columns, int p_thread_id, int p_batch_size=1000){
   _batch_size = p_batch_size;
   _thread_id = p_thread_id;
   _input_parameters = ref_input_op;
   stored_records.reserve(p_batch_size+1);
   for (int _i_=0; _i_<p_batch_size+1; _i_++){
      stored_records.push_back( SeqTxt_SS_Record() );
   }
   _header_columns = header_columns;
}

SeqTxt_Thread_data::~SeqTxt_Thread_data(){ 
}

std::map<std::string, int> SeqTxt_Thread_data::getHeaderColumns()
{
    return _header_columns;
}

size_t SeqTxt_Thread_data::read_ss_record(std::ifstream& file_stream, std::map<std::string, int> header_columns){
    //std::cout << "Type 1." << std::endl;
    thread_index = 0;  // Index where this thread's data will be stored
    while( std::getline( file_stream, current_line ) )
    {
        std::istringstream column_stream( current_line );

        // Read each column value from the record line
        std::string column_value;
        std::vector<std::string> column_values;
        while (std::getline( column_stream, column_value, '\t' )) {
            column_values.push_back(column_value);
        }

        // Store the required values from the record
        SeqTxt_SS_Record record_data;

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

// SeqTxt_Module Constructor
SeqTxt_Module::SeqTxt_Module(Input_Para& input_parameters){
    _input_parameters = input_parameters;
    has_error = 0;
    file_index = 0;

    // input_file_stream = NULL;
    if (file_index >= _input_parameters.num_input_files){
        std::cerr << "Input file list error." << std::endl;
        has_error |= 1;
        return;
    }

    // Open the first file in the list
    const char * first_filepath = _input_parameters.input_files[file_index].c_str();
    // input_file_stream = new std::ifstream(first_filepath);
    input_file_stream.open(first_filepath);
    if (!(input_file_stream.is_open())){
        std::cerr << "Cannot open sequencing_summary.txt file="<< first_filepath <<std::endl;
        has_error |= 2;
    }else{
        file_index ++;
        std::cout<< "INFO: Open sequencing_summary.txt file = "<< first_filepath <<" " << file_index<<"/"<<_input_parameters.num_input_files <<std::endl;

        // Ensure that we have the columns we need for statistics
        std::string column_line;
        std::getline( input_file_stream, column_line );
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

bool SeqTxt_Module::requiredHeadersFound(std::string header_string) {
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

SeqTxt_Module::~SeqTxt_Module(){
    if (input_file_stream.is_open()){
        input_file_stream.close();
    }
}

int SeqTxt_Module::generateStatistics( Output_SeqTxt& t_output_SeqTxt_info){  
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

    t_output_SeqTxt_info.all_long_read_info.long_read_info.resize();
    t_output_SeqTxt_info.passed_long_read_info.long_read_info.resize();
    t_output_SeqTxt_info.failed_long_read_info.long_read_info.resize();
    printMemoryUsage("Before generating statistics");

    if (has_error==0) {
        m_threads.reserve(_input_parameters.threads+3);

        int _i_t=0;
        printMessage("Generating statistics...");
        try {
            for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
                m_threads.push_back(std::thread((SeqTxt_Module::SeqTxt_do_thread), std::ref(input_file_stream), std::ref(_input_parameters), _i_t, std::ref(t_output_SeqTxt_info), _header_columns, SeqTxt_Module::batch_size_of_record ));
            }
            printMessage("Joining " + std::to_string(_input_parameters.threads) + " threads...");
            for (_i_t=0; _i_t<_input_parameters.threads; _i_t++){
                m_threads[_i_t].join();
            }

        }catch(const std::runtime_error& re){
            printError("Runtime error: " + std::string(re.what()));
        }catch(const std::exception& ex){
            std::cerr << "Error occurred: " << ex.what() << std::endl;
            printError("Error: " + std::string(ex.what()));
        }catch(...){
            printError("Unknown error occurred in thread " + std::to_string(_i_t));
        }
    }
    t_output_SeqTxt_info.global_sum();
    printMemoryUsage("After generating statistics");
 
    auto relapse_end_time = std::chrono::high_resolution_clock::now();
    std::cout<<"Elapsed time (seconds): "<< std::chrono::duration_cast<std::chrono::seconds>(relapse_end_time - relapse_start_time).count() << std::endl;

    std::cout<<"sequencing_summary.txt QC "<< (has_error==0?"generated":"failed") << std::endl;
 
    return has_error;
}

void SeqTxt_Module::SeqTxt_do_thread(std::ifstream& file_stream, Input_Para& ref_input_op, int thread_id, Output_SeqTxt& ref_output, std::map<std::string, int> header_columns, size_t batch_size_of_record){
    size_t read_ss_size, read_ss_i;
    int total_read_count = 0;
    while (true) {
        SeqTxt_Thread_data ref_thread_data(ref_input_op, header_columns, thread_id, batch_size_of_record);
        {
            std::lock_guard<std::mutex> lock(myMutex_readSeqTxt);
            std::map<std::string, int> header_column_data = ref_thread_data.getHeaderColumns();
            read_ss_size = ref_thread_data.read_ss_record(file_stream, header_column_data);

            if (read_ss_size == 0 && !(file_index < ref_input_op.num_input_files) ){
                break;
            }

            if ( read_ss_size < batch_size_of_record ){
                if ( file_index < ref_input_op.num_input_files ){
                    file_stream.close();
                    file_stream.clear();

                    file_stream.open( ref_input_op.input_files[file_index].c_str() );
                    std::string firstline;
                    std::getline( file_stream, firstline );
                    file_index++;
                }
            }
        }
        if (read_ss_size == 0 ) {
            continue;
        } else {
            total_read_count += read_ss_size;
            printMessage("Thread " + std::to_string(thread_id+1) + " read " + std::to_string(read_ss_size) + " records (total " + std::to_string(total_read_count) + ")");
        }

        // Columns used for statistics: passes_filtering, sequence_length_template, mean_qscore_template
        ref_thread_data.t_output_SeqTxt_.all_long_read_info.long_read_info.resize();
        ref_thread_data.t_output_SeqTxt_.passed_long_read_info.long_read_info.resize();
        ref_thread_data.t_output_SeqTxt_.failed_long_read_info.long_read_info.resize();
        for(read_ss_i=0; read_ss_i<read_ss_size; read_ss_i++){
            bool passes_filtering_value = ref_thread_data.stored_records[read_ss_i].passes_filtering;
            Basic_SeqTxt_Statistics& seqtxt_statistics = (passes_filtering_value == true) ? ref_thread_data.t_output_SeqTxt_.passed_long_read_info : ref_thread_data.t_output_SeqTxt_.failed_long_read_info;
            
            seqtxt_statistics.long_read_info.total_num_reads++;
            size_t sequence_base_count = ref_thread_data.stored_records[read_ss_i].sequence_length_template;
            seqtxt_statistics.long_read_info.total_num_bases += sequence_base_count;

            // Store the read length
            seqtxt_statistics.long_read_info.read_lengths.push_back(sequence_base_count);

            // Update the longest read length
            int64_t current_read_length = (int64_t) ref_thread_data.stored_records[read_ss_i].sequence_length_template;
            if ( seqtxt_statistics.long_read_info.longest_read_length < current_read_length){
                seqtxt_statistics.long_read_info.longest_read_length = current_read_length;
            }
            seqtxt_statistics.long_read_info.read_length_count[ ref_thread_data.stored_records[read_ss_i].sequence_length_template<MAX_READ_LENGTH?ref_thread_data.stored_records[read_ss_i].sequence_length_template:(MAX_READ_LENGTH-1) ] += 1;

            seqtxt_statistics.seq_quality_info.read_quality_distribution[ int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template ) ] += 1;
            if ( seqtxt_statistics.seq_quality_info.min_read_quality == MoneDefault ||
                seqtxt_statistics.seq_quality_info.min_read_quality>int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template ) ){
                seqtxt_statistics.seq_quality_info.min_read_quality = int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template );
            }
            if ( seqtxt_statistics.seq_quality_info.max_read_quality < int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template) ){
                seqtxt_statistics.seq_quality_info.max_read_quality = int( ref_thread_data.stored_records[read_ss_i].mean_qscore_template);
            }
        }
        std::lock_guard<std::mutex> lock(myMutex_output);
        ref_output.add( ref_thread_data.t_output_SeqTxt_ );
    }
}
