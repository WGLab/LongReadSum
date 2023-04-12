#ifndef BAM_MODULE_H_
#define BAM_MODULE_H_

#include <thread>
#include <mutex>
#include <map>

#include "InputStructure.h"
#include "OutputStructures.h"
#include "BamReader.h"


class BAM_Thread_data {
   public:
      int _thread_id;
//      std::vector<Bam1Record> record_list;  // Vector of individual BAM file records
      Input_Para m_input_op;
      std::map<std::string, bool> t_secondary_alignment;
      std::map<std::string, bool> t_supplementary_alignment;

      Output_BAM t_output_bam_;

      BAM_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size);
      ~BAM_Thread_data();
};

class BAM_Module{
public:
    std::mutex bam_mutex;
    std::mutex output_mutex;

    static const size_t records_per_thread = 3000;
    size_t file_index = 0;

    std::map<std::string, bool> secondary_alignment;
    std::map<std::string, bool> supplementary_alignment;

    Input_Para m_input_op;

    BamReader* _bam_reader_ptr;
    std::vector<std::thread> m_threads;


    int exit_code;  // Exit code for error handling

    static void BAM_do_thread(BamReader* ref_bam_reader_ptr, Input_Para& ref_input_op, int thread_id, BAM_Thread_data& ref_thread_data, Output_BAM& ref_output, std::map<std::string, bool>& ref_secondary_alignment, std::map<std::string, bool>& ref_supplementary_alignment, std::mutex& bam_mutex, std::mutex& output_mutex, size_t &file_index);

    int calculateStatistics( Output_BAM& t_output_bam_info);
    static void batchStatistics(BamReader* bam_reader, Input_Para& input_params, Output_BAM& ref_output, std::mutex& bam_mutex, std::mutex& output_mutex);

    BAM_Module(Input_Para& _m_input);
    ~BAM_Module();
};


#endif

