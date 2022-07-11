#ifndef HDF5_MODULE_H_
#define HDF5_MODULE_H_

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <thread>
#include <mutex>

#include "CXX_to_python_class.h"
#include "Python_to_CXX_class.h"
#include "H5Cpp.h"
using namespace H5;


class FAST5_Record{
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
        std::vector<FAST5_Record> stored_records;
        Output_FAST5 output_data;
        size_t readRecord(std::ifstream* file_stream);

        FAST5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size);
        ~FAST5_Thread_data();
};

// Main class for FAST5 file statistics generation
class FAST5_Module{
    private:
        static size_t file_index;  // Tracks the current file index being analyzed
    public:
        H5File input_fast5;  // HDF5 file object
        int has_error;  // Error code
        static std::mutex myMutex_readFAST5;
        static std::mutex myMutex_output;
        static size_t batch_size_of_record;
        Input_Para _input_parameters;
        std::string input_filepath;  // The current input file
        std::vector<std::thread> all_threads;

        // Methods
        // Assign threads
        static void createThread(std::string input_filepath, Input_Para& ref_input_op, int thread_id, FAST5_Thread_data& ref_thread_data, Output_FAST5& ref_output);

        // Generate statistics
        int generateStatistics( Output_FAST5& output_datainfo);

        FAST5_Module(Input_Para& _m_input);
        ~FAST5_Module();
};

#endif
