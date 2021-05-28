#include "F5_module.h"
#include "ComFunction.h"

#include <iostream>

size_t F5_Module::batch_size_of_record=3000;
std::mutex F5_Module::myMutex_readF5;
std::mutex F5_Module::myMutex_output;
size_t F5_Module::read_i_F5 = 0;

F5_Thread_data::F5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   _batch_size = p_batch_size;
   _thread_id = p_thread_id;
   m_input_op = ref_input_op;
   _F5_ss_records.reserve(p_batch_size+1);
   for (int _i_=0; _i_<p_batch_size+1; _i_++){
      _F5_ss_records.push_back( F5_SS_Record() );
   }
}

F5_Thread_data::~F5_Thread_data(){ 
   ; 
   //std::cout<<" ~F5_Thread_data " << _thread_id <<std::endl<<std::flush;
}

size_t F5_Thread_data::read_ss_record(std::ifstream* ref_F5_reader_ss){
   num_ss_read_record = 0;
   //_F5_ss_records.clear();

   while( std::getline( *ref_F5_reader_ss, _line_seq_sum )) {
       std::istringstream iss( _line_seq_sum );
       F5_SS_Record _t_f5_ss_record;
       if (!(iss >> _F5_ss_records[num_ss_read_record].filename >> _F5_ss_records[num_ss_read_record].read_id 
                 >> _F5_ss_records[num_ss_read_record].run_id >> _F5_ss_records[num_ss_read_record].channel 
                 >> _F5_ss_records[num_ss_read_record].start_time >> _F5_ss_records[num_ss_read_record].duration 
                 >> _F5_ss_records[num_ss_read_record].num_events >> _F5_ss_records[num_ss_read_record].passes_filtering 
                 >> _F5_ss_records[num_ss_read_record].template_start >> _F5_ss_records[num_ss_read_record].num_events_template 
                 >> _F5_ss_records[num_ss_read_record].template_duration >> _F5_ss_records[num_ss_read_record].num_called_template
                 >> _F5_ss_records[num_ss_read_record].sequence_length_template >> _F5_ss_records[num_ss_read_record].mean_qscore_template >> _F5_ss_records[num_ss_read_record].strand_score_template
                 >> _F5_ss_records[num_ss_read_record].calibration_strand_genome_template >> _F5_ss_records[num_ss_read_record].calibration_strand_identity_template
                 >> _F5_ss_records[num_ss_read_record].calibration_strand_accuracy_template >> _F5_ss_records[num_ss_read_record].calibration_strand_speed_bps_template)) {
             std::cout<<"Error!!! for <"<<_line_seq_sum<<">"<<std::endl;
             break;
       }

       num_ss_read_record++;
       if ( num_ss_read_record >= _batch_size){ break; }
   }

   return num_ss_read_record;
}

F5_Module::F5_Module(Input_Para& _m_input){
   m_input_op = _m_input;

   has_error = 0;
   read_i_F5 = 0;

   _F5_reader_ss = NULL;   
   if ( read_i_F5 >= _m_input.num_input_files ){
      std::cout<<"Error!!! No input F5 are find. #input_files="<< _m_input.num_input_files <<std::endl;
      has_error |= 1;
      return;
   }
   
    _F5_reader_ss = new std::ifstream ( _m_input.input_files[read_i_F5].c_str() );
   if (!(_F5_reader_ss->is_open())){
       std::cout<< "Error!!! Cannot open F5 file="<< _m_input.input_files[read_i_F5] <<std::endl;
       has_error |= 2;
   }else{
       std::cout<< "INFO: Open F5 file="<< _m_input.input_files[read_i_F5] <<" " << read_i_F5<<"/"<<_m_input.num_input_files <<std::endl;
       std::string firstline;
       std::getline( *_F5_reader_ss, firstline );
       read_i_F5 = 1;
   }
}

F5_Module::~F5_Module(){
   //std::cout<<" ~F5_Module begin" <<std::endl<<std::flush;
   if (_F5_reader_ss!=NULL){ 
       delete _F5_reader_ss; 
   }
   _F5_reader_ss = NULL;
   //std::cout<<" ~F5_Module end" <<std::endl<<std::flush;
}

/*
Output_F5 F5_Module::F5_st(){
   Output_F5 t_output_F5_info;*/
int F5_Module::F5_st( Output_F5& t_output_F5_info){  
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   if (has_error==0){
       m_threads.reserve(m_input_op.threads+3);

      int _i_t=0;
      F5_Thread_data** thread_data_vector = new F5_Thread_data*[m_input_op.threads];
      try{
         for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
             std::cout<<"INFO: generate threads "<<_i_t<<std::endl<<std::flush;
             thread_data_vector[_i_t] = new F5_Thread_data(m_input_op, _i_t, F5_Module::batch_size_of_record);
             std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
             m_threads.push_back(std::thread((F5_Module::F5_do_thread), _F5_reader_ss, std::ref(m_input_op), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(t_output_F5_info) ));
         }

         std::cout<<"INFO: join threads"<<std::endl<<std::flush;
         for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
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
     
      // std::cout<<" del "<<std::endl<<std::flush; 
      for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
         // std::cout<<" del >"<< _i_t <<std::endl<<std::flush;
         delete thread_data_vector[_i_t];
         // std::cout<<" del <"<< _i_t <<std::endl<<std::flush;
      }
      // std::cout<<" del2 "<<std::endl<<std::flush;
      delete [] thread_data_vector;
   }

   t_output_F5_info.global_sum();
 
   auto relapse_end_time = std::chrono::high_resolution_clock::now();
   std::cout<<"Total time(replase): "<<round3((relapse_end_time-relapse_start_time).count()/1000000000.0)<<std::endl;

   std::cout<<"<Stastitics on F5>: "<< (has_error==0?"successfully":"Failed") <<"."<<std::endl;
 
   //return t_output_F5_info; 
   return has_error;
}

void F5_Module::F5_do_thread(std::ifstream* ref_F5_reader_ss, Input_Para& ref_input_op, int thread_id, F5_Thread_data& ref_thread_data, Output_F5& ref_output ){
    size_t read_ss_size, read_ss_i;
    while (true){
        //ref_thread_data._F5_ss_records.clear();
        myMutex_readF5.lock();
        //
        read_ss_size = ref_thread_data.read_ss_record(ref_F5_reader_ss);

        if (read_ss_size == 0 && !(read_i_F5 < ref_input_op.num_input_files) ){
            myMutex_readF5.unlock();
            break;
        }
        if ( read_ss_size < batch_size_of_record ){
            if ( read_i_F5 < ref_input_op.num_input_files ){ 
               std::cout<< "INFO: Open F5 file="<< ref_input_op.input_files[read_i_F5] <<std::endl;
               ref_F5_reader_ss->close();
               ref_F5_reader_ss->clear();

               ref_F5_reader_ss->open( ref_input_op.input_files[read_i_F5].c_str() );
               std::string firstline;
               std::getline( *ref_F5_reader_ss, firstline );
               read_i_F5++;
            }
        }
        myMutex_readF5.unlock();

        if (read_ss_size == 0 ) { continue; }       
 
        ref_thread_data.t_output_F5_.reset();
        for(read_ss_i=0; read_ss_i<read_ss_size; read_ss_i++){
           Basic_F5_Statistics* _f5_st = NULL;
           if ( ref_thread_data._F5_ss_records[read_ss_i].passes_filtering.compare("True")!=0){
               _f5_st = &(ref_thread_data.t_output_F5_.f5_passed_long_read_info);
           }else{
               _f5_st = &(ref_thread_data.t_output_F5_.f5_failed_long_read_info);
           }
           _f5_st->long_read_info.total_num_reads++;
           _f5_st->long_read_info.total_num_bases += ref_thread_data._F5_ss_records[read_ss_i].sequence_length_template;
           if ( _f5_st->long_read_info.longest_read_length < ref_thread_data._F5_ss_records[read_ss_i].sequence_length_template){
               _f5_st->long_read_info.longest_read_length = ref_thread_data._F5_ss_records[read_ss_i].sequence_length_template;
           }
           _f5_st->long_read_info.read_length_count[ ref_thread_data._F5_ss_records[read_ss_i].sequence_length_template<MAX_READ_LENGTH?ref_thread_data._F5_ss_records[read_ss_i].sequence_length_template:(MAX_READ_LENGTH-1) ] += 1;
           
           _f5_st->long_read_info.read_quality_distribution[ int( ref_thread_data._F5_ss_records[read_ss_i].mean_qscore_template ) ] += 1;
           if ( _f5_st->long_read_info.min_read_quality == MoneDefault || 
               _f5_st->long_read_info.min_read_quality>int( ref_thread_data._F5_ss_records[read_ss_i].mean_qscore_template ) ){
              _f5_st->long_read_info.min_read_quality = int( ref_thread_data._F5_ss_records[read_ss_i].mean_qscore_template );
           }
           if ( _f5_st->long_read_info.max_read_quality < int( ref_thread_data._F5_ss_records[read_ss_i].mean_qscore_template) ){
              _f5_st->long_read_info.max_read_quality = int( ref_thread_data._F5_ss_records[read_ss_i].mean_qscore_template); 
           }
        }

        myMutex_output.lock();
       
        ref_output.add( ref_thread_data.t_output_F5_ );

        myMutex_output.unlock();
    }
}


