#ifndef HDF5_MODULE_H_
#define HDF5_MODULE_H_

#include "CXX_to_python_class.h"
#include "Python_to_CXX_class.h"

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <thread>
#include <mutex>
#include <map>

int runExample(void);

class FAST5_SS_Record{
    public:
        std::string filename;
        std::string read_id;
        std::string run_id;
        bool passes_filtering;
        size_t sequence_length_template;
        float mean_qscore_template;
};

class FAST5_Thread_data {
private:
    int thread_index;  // Index where this thread's data will be stored

    public:
        int _thread_id;
        int _batch_size;
        Input_Para _input_parameters;
        std::vector<FAST5_SS_Record> stored_records;
        std::map<std::string, int> _header_columns;

        Output_FAST5 output_data;
        std::string current_line;  // Current line being read from the file

        size_t read_ss_record(std::ifstream* file_stream);


        FAST5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size);
        ~FAST5_Thread_data();
};

// Main class for FAST5 file statistics generation
class FAST5_Module{
    private:
        static size_t file_index;  // Tracks the current file index being analyzed
    public:
        static std::mutex myMutex_readFAST5;
        static std::mutex myMutex_output;
        static size_t batch_size_of_record;

        Input_Para _input_parameters;

        std::ifstream *input_file_stream;  // Stream for the input text file
        std::vector<std::thread> m_threads;


        int has_error;

        // Methods
        // Assign threads
        static void createThread(std::ifstream* file_stream, Input_Para& ref_input_op, int thread_id, FAST5_Thread_data& ref_thread_data, Output_FAST5& ref_output);

        // Generate statistics
        int generateStatistics( Output_FAST5& output_datainfo);

        FAST5_Module(Input_Para& _m_input);
        ~FAST5_Module();
};

#endif
