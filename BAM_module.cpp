#include "BAM_module.h"

BAM_Thread_data::BAM_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   br_list.reserve(p_batch_size+1);
   _thread_id = p_thread_id;
   m_input_op = ref_input_op;
}

BAM_Thread_data::~BAM_Thread_data(){ ; }

int BAM_Module::batch_size_of_record=3000;
std::mutex BAM_Module::myMutex_readBam;
std::mutex BAM_Module::myMutex_output;

BAM_Module::BAM_Module(Input_Para _m_input){
   m_input_op = _m_input;

   has_error = 0;
   read_i_bam = 0;

   _bam_reader_ptr = NULL;   
   if ( read_i_bam >= _m_input.num_input_files ){
      std::cout<<"Error!!! No input BAM are find. #input_files="<<num_input_files<<std::endl;
      has_error |= 1;
      return;
   }
   
   _bam_reader_ptr = new BamReader( input_files[read_i_bam] );
   if (!_bam_reader_ptr->check_bam_status()){
       std::cout<< "Error!!! Cannot open bam file="<< input_files[read_i_bam] <<std::endl;
       has_error |= 2;
   }else{
       std::cout<< "INFO: Open bam file="<< input_files[read_i_bam] <<std::endl;
       _bam_reader_ptr->bamReadOp.set_w_supplementary(true);
       _bam_reader_ptr->bamReadOp.set_w_secondary(true);
       _bam_reader_ptr->bamReadOp.set_min_read_len(1);
       _bam_reader_ptr->bamReadOp.set_w_qry_seq(true);
       _bam_reader_ptr->bamReadOp.set_w_unmap(true);
       _bam_reader_ptr->bamReadOp.set_w_qry_qual(true);
       read_i_bam = 1;
   }
}

BAM_Module::~BAM_Module(){
   if (_bam_reader_ptr!=NULL){ delete _bam_reader_ptr; }
   _bam_reader_ptr = NULL;
}

Output_BAM BAM_Module::bam_st(){
   Output_BAM t_output_bam_info;
  
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   if (has_error==0){
       m_threads.reserve(m_input_op.threads+3);

      int _i_t=0;
      BAM_Thread_data** thread_data_vector = new BAM_Thread_data*[deepmod_op.threads];
      try{
         for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
             std::cout<<"INFO: generate threads "<<_i_t<<std::endl<<std::flush;
             thread_data_vector[_i_t] = new BAM_Thread_data(m_input_op, _i_t);
             std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
             m_threads.push_back(std::thread((DeepMod_Run::BAM_do_thread), _bam_reader_ptr, std::ref(m_input_op), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(t_output_bam_info) ));
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
      
      for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
         delete thread_data_vector[_i_t];
      }
      delete [] thread_data_vector;
   }
 
   t_output_bam_info.global_sum();
 
   auto relapse_end_time = std::chrono::high_resolution_clock::now();
   std::cout<<"Total time(replase): "<<round3((relapse_end_time-relapse_start_time).count()/1000000000.0)<<std::endl;

   std::cout<<"<Stastitics on BAM>: "<< (has_error==0?"successfully":"Failed") <<"."<<std::endl;
 
   return t_output_bam_info; 
}

void BAM_Module::BAM_do_thread(BamReader& ref_bam_reader, Input_Para& ref_input_op, int thread_id, BAM_Thread_data& ref_bam_data, Output_BAM& ref_output){
    std::vector<Bam1Record>::iterator br_it;
    while (true){
        //auto start_time = std::chrono::high_resolution_clock::now();

        ref_bam_data.br_lists.clear();
        myMutex_readBam.lock();
        _bam_reader_ptr->readBam( ref_bam_data.br_lists );
        if (ref_bam_data.br_lists.size() == 0 && !(read_i_bam < ref_input_op.num_input_files) ){
            break;
        }
        if ( ref_bam_data.br_lists.size() < batch_size_of_record ){
            if ( read_i_bam < ref_input_op.num_input_files ){ 
               std::cout<< "INFO: Open bam file="<< input_files[read_i_bam] <<std::endl;
               _bam_reader_ptr->resetBam( input_files[read_i_bam] );
               read_i_bam++;
            }
        }
        myMutex_readBam.unlock();

      Output_BAM _output_bam_;

        ref_thread_data.t_secondary_alignment.clear();
        ref_thread_data.t_supplementary_alignment.clear();
        for(br_it=ref_thread_data.br_list.begin(); br_it!=ref_thread_data.br_list.end(); br_it++){
           if ( (br_it->map_flag)& BAM_FUNMAP ) { ; }
           else if ( (br_it->map_flag)& BAM_FSUPPLEMENTARY ) { ; }
           else if ( (br_it->map_flag)& BAM_FSECONDARY ) { ; }

           if (br_it->map_strand==0){ ; }
           else { ; }

           

        }

    }
}


