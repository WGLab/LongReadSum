#ifndef F5_MODULE_H_
#define F5_MODULE_H_

#include "CXX_to_python_class.h"
#include "Python_to_CXX_class.h"

#include <thread>
#include <mutex>
#include <map>

class F5_Thread_data {
   public:
      int _thread_id;
      Input_Para m_input_op;

      Output_F5 t_output_F5_;

      F5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size);
      ~F5_Thread_data();
};

class F5_Module{
private:
      static size_t read_i_F5;

public:
      static std::mutex myMutex_readF5;
      static std::mutex myMutex_output;
      static size_t batch_size_of_record;

      Input_Para m_input_op;

      /// need correct
      F5Reader* _F5_reader_ptr;
      std::vector<std::thread> m_threads;


   int has_error;

   static void F5_do_thread(F5Reader* ref_F5_reader_ptr, Input_Para& ref_input_op, int thread_id, F5_Thread_data& ref_thread_data, Output_F5& ref_output);   

   int F5_st( Output_F5& t_output_F5_info);

   F5_Module(Input_Para& _m_input);
   ~F5_Module();
};


#endif

