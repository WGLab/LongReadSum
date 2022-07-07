#ifndef BAM_MODULE_H_
#define BAM_MODULE_H_

#include "CXX_to_python_class.h"
#include "Python_to_CXX_class.h"

#include <thread>
#include <mutex>
#include "BamReader.h"
#include <map>

class BAM_Thread_data {
   public:
      int _thread_id;
      std::vector<Bam1Record> br_list;
      Input_Para m_input_op;
      std::map<std::string, bool> t_secondary_alignment;
      std::map<std::string, bool> t_supplementary_alignment;

      Output_BAM t_output_bam_;

      BAM_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size);
      ~BAM_Thread_data();
};

class BAM_Module{
private:
      static size_t read_i_bam;

public:
      static std::mutex myMutex_readBam;
      static std::mutex myMutex_output;
      static size_t batch_size_of_record;

      std::map<std::string, bool> secondary_alignment;
      std::map<std::string, bool> supplementary_alignment;

      Input_Para m_input_op;

      BamReader* _bam_reader_ptr;
      std::vector<std::thread> m_threads;


   int exit_code;  // Exit code for error handling

   static void BAM_do_thread(BamReader* ref_bam_reader_ptr, Input_Para& ref_input_op, int thread_id, BAM_Thread_data& ref_thread_data, Output_BAM& ref_output, std::map<std::string, bool>& ref_secondary_alignment, std::map<std::string, bool>& ref_supplementary_alignment);   

   int calculateStatistics( Output_BAM& t_output_bam_info);

   BAM_Module(Input_Para& _m_input);
   ~BAM_Module();
};


#endif

