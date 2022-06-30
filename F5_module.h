#ifndef F5_MODULE_H_
#define F5_MODULE_H_

#include "CXX_to_python_class.h"
#include "Python_to_CXX_class.h"

#include <fstream>
#include <sstream>
#include <vector>

#include <thread>
#include <mutex>
#include <map>

class F5_SS_Record{
public:
      std::string filename;
      std::string read_id;
      std::string run_id;
      int channel;
      float start_time;
      float duration;
      int num_events;
      std::string passes_filtering;
      float template_start;
      size_t num_events_template;
      float template_duration;
      size_t num_called_template;
      size_t sequence_length_template;
      float mean_qscore_template;
      float strand_score_template;
      std::string calibration_strand_genome_template;
      float calibration_strand_identity_template;
      float calibration_strand_accuracy_template;
      float calibration_strand_speed_bps_template;
};

class F5_Thread_data {
private:
   int num_ss_read_record;

   public:
      int _thread_id;
      int _batch_size;
      Input_Para m_input_op;
      std::vector<F5_SS_Record> _F5_ss_records;

      Output_F5 t_output_F5_;
      //char* _line_seq_sum[4096];
      std::string _line_seq_sum;

      size_t read_ss_record(std::ifstream* ref_F5_reader_ss);
      size_t read_ss_record_2(std::ifstream* ref_F5_reader_ss);
      size_t read_ss_record_3(std::ifstream* ref_F5_reader_ss);

      F5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size);
      ~F5_Thread_data();
};

// Main class for FAST5 file statistics generation
class F5_Module{
private:
      static size_t read_i_F5;

public:
      static std::mutex myMutex_readF5;
      static std::mutex myMutex_output;
      static size_t batch_size_of_record;

      Input_Para m_input_op;

      /// need correct
      std::ifstream *_F5_reader_ss;
      std::vector<std::thread> m_threads;


   int has_error;


   static void F5_do_thread(std::ifstream* ref_F5_reader_ss, Input_Para& ref_input_op, int thread_id, F5_Thread_data& ref_thread_data, Output_F5& ref_output);   

   // Function for generating statistics
   int F5_st( Output_F5& t_output_F5_info);

   F5_Module(Input_Para& _m_input);
   ~F5_Module();
};


#endif

