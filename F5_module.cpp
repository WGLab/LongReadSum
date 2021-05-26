#include "F5_module.h"
#include "ComFunction.h"

F5_Thread_data::F5_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   br_list.reserve(p_batch_size+1);
   _thread_id = p_thread_id;
   m_input_op = ref_input_op;
}

F5_Thread_data::~F5_Thread_data(){ 
   ; 
   //std::cout<<" ~F5_Thread_data " << _thread_id <<std::endl<<std::flush;
}

size_t F5_Module::batch_size_of_record=3000;
std::mutex F5_Module::myMutex_readF5;
std::mutex F5_Module::myMutex_output;
size_t F5_Module::read_i_F5 = 0;

F5_Module::F5_Module(Input_Para& _m_input){
   m_input_op = _m_input;

   has_error = 0;
   read_i_F5 = 0;

   _F5_reader_ptr = NULL;   
   if ( read_i_F5 >= _m_input.num_input_files ){
      std::cout<<"Error!!! No input F5 are find. #input_files="<< _m_input.num_input_files <<std::endl;
      has_error |= 1;
      return;
   }
   
    _F5_reader_ptr = new F5Reader( _m_input.input_files[read_i_F5].c_str() );
   if (!_F5_reader_ptr->check_F5_status()){
       std::cout<< "Error!!! Cannot open F5 file="<< _m_input.input_files[read_i_F5] <<std::endl;
       has_error |= 2;
   }else{
       std::cout<< "INFO: Open F5 file="<< _m_input.input_files[read_i_F5] <<" " << read_i_F5<<"/"<<_m_input.num_input_files <<std::endl;
       read_i_F5 = 1;
   }
}

F5_Module::~F5_Module(){
   //std::cout<<" ~F5_Module begin" <<std::endl<<std::flush;
   if (_F5_reader_ptr!=NULL){ 
       delete _F5_reader_ptr; 
   }
   _F5_reader_ptr = NULL;
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
             thread_data_vector[_i_t] = new F5_Thread_data(m_input_op, _i_t);
             std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
             m_threads.push_back(std::thread((F5_Module::F5_do_thread), _F5_reader_ptr, std::ref(m_input_op), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(t_output_F5_info) ));
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

void F5_Module::F5_do_thread(F5Reader* ref_F5_reader_ptr, Input_Para& ref_input_op, int thread_id, F5_Thread_data& ref_thread_data, Output_F5& ref_output ){
    std::vector<F51Record>::iterator br_it;
    uint64_t match_this_read;
    double accuracy_perc_this_read;
    while (true){
        //auto start_time = std::chrono::high_resolution_clock::now();

        ref_thread_data.br_list.clear();
        myMutex_readF5.lock();
        ref_F5_reader_ptr->readF5( ref_thread_data.br_list, F5_Module::batch_size_of_record );
        if (ref_thread_data.br_list.size() == 0 && !(read_i_F5 < ref_input_op.num_input_files) ){
            // std::cout<<" " << ref_thread_data.br_list.size()<<" Reads." << read_i_F5 <<"/"<< ref_input_op.num_input_files <<std::endl;
            myMutex_readF5.unlock();
            break;
        }
        if ( ref_thread_data.br_list.size() < batch_size_of_record ){
            if ( read_i_F5 < ref_input_op.num_input_files ){ 
               std::cout<< "INFO: Open F5 file="<< ref_input_op.input_files[read_i_F5] <<std::endl;
               ref_F5_reader_ptr->resetF5( ref_input_op.input_files[read_i_F5].c_str() );
               //ref_F5_reader_ptr->resetF5( ref_input_op.input_files[read_i_F5] );
               if ( ref_thread_data.br_list.size() == 0 ){
                  ref_F5_reader_ptr->readF5( ref_thread_data.br_list, F5_Module::batch_size_of_record );
               }
               read_i_F5++;
            }
        }
        myMutex_readF5.unlock();

        if (ref_thread_data.br_list.size() == 0 ) { continue; }       
 
        ref_thread_data.t_output_F5_.reset();
        //std::cout<<" " << ref_thread_data.br_list.size()<<" Reads." <<std::endl; 
        for(br_it=ref_thread_data.br_list.begin(); br_it!=ref_thread_data.br_list.end(); br_it++){
           // std::cout<< br_it->qry_seq_len <<" " << br_it->qry_seq.length() << " " << br_it->qry_name << " " << int(br_it->map_flag) << std::endl <<std::flush;
           Basic_Seq_Statistics* _seq_st = NULL;
           Basic_Seq_Quality_Statistics* _seq_qual_st = NULL;
           if ( (br_it->map_flag)& F5_FUNMAP ) { 
               // ref_thread_data.t_output_F5_.unmapped_long_read_info.total_num_reads += 1; 
               // ref_thread_data.t_output_F5_.unmapped_long_read_info.total_num_bases += br_it->qry_seq_len;
               _seq_st = &( ref_thread_data.t_output_F5_.unmapped_long_read_info);
               _seq_qual_st = &(ref_thread_data.t_output_F5_.unmapped_seq_quality_info);
           }else if ( (br_it->map_flag)& F5_FSUPPLEMENTARY ) { 
               ref_thread_data.t_output_F5_.num_supplementary_alignment += 1; 
               ref_thread_data.t_supplementary_alignment[ br_it->qry_name ] = true;
           }else if ( (br_it->map_flag)& F5_FSECONDARY ) { 
               ref_thread_data.t_output_F5_.num_secondary_alignment += 1;
               ref_thread_data.t_secondary_alignment[ br_it->qry_name ] = true; 
           } else {
               ref_thread_data.t_output_F5_.num_primary_alignment += 1;
               // ref_thread_data.t_output_F5_.mapped_long_read_info.total_num_reads += 1;
               // ref_thread_data.t_output_F5_.mapped_long_read_info.total_num_bases += br_it->qry_seq_len;
               _seq_st = &( ref_thread_data.t_output_F5_.mapped_long_read_info);
               _seq_qual_st = &(ref_thread_data.t_output_F5_.mapped_seq_quality_info);
           }
           
           if ( !( (br_it->map_flag)& F5_FUNMAP ) ){
              if (br_it->map_strand==0){ ref_thread_data.t_output_F5_.forward_alignment += 1; }
              else { ref_thread_data.t_output_F5_.reverse_alignment += 1; }
           }

           if ( br_it->map_qual < ref_thread_data.t_output_F5_.min_map_quality || ref_thread_data.t_output_F5_.min_map_quality==MoneDefault){ 
               ref_thread_data.t_output_F5_.min_map_quality = br_it->map_qual; }
           if ( br_it->map_qual > ref_thread_data.t_output_F5_.max_map_quality ){ ref_thread_data.t_output_F5_.max_map_quality = br_it->map_qual; }
           ref_thread_data.t_output_F5_.map_quality_distribution[ int(br_it->map_qual) ] += 1;
           
           if (_seq_st!= NULL && !( (br_it->map_flag)& F5_FSECONDARY ) && !((br_it->map_flag)& F5_FSUPPLEMENTARY ) ){
              _seq_st->total_num_reads += 1;
              _seq_st->total_num_bases += br_it->qry_seq_len;
              if( br_it->qry_seq_len > _seq_st->longest_read_length){ _seq_st->longest_read_length = br_it->qry_seq_len; }
              if (  br_it->qry_seq_len < MAX_READ_LENGTH ){ _seq_st->read_length_count[br_it->qry_seq_len] += 1; }
              
              for (size_t _si=0; _si<br_it->qry_seq_len; _si++){
                 int _base_qual_int = int(br_it->qry_qual[_si]);
                 if ( _base_qual_int > MAX_BASE_QUALITY || _base_qual_int <0){
                    std::cout<<"WARNINING!!! from "<< thread_id << " The base quality is not in [0,"<< MAX_BASE_QUALITY <<"]: " << int(br_it->qry_qual[_si]) << " at " << _si << " for " <<  br_it->qry_name  <<std::endl;
                 }else{
                    _seq_qual_st->base_quality_distribution[_base_qual_int ] += 1;
                    if ( _seq_qual_st->min_base_quality==MoneDefault || _base_qual_int < _seq_qual_st->min_base_quality){
                        _seq_qual_st->min_base_quality = _base_qual_int; }
                    if ( _base_qual_int > _seq_qual_st->max_base_quality) { _seq_qual_st->max_base_quality = _base_qual_int; }
                 }

                 switch( br_it->qry_seq[ _si ] ){
                     case 'A': case 'a':
                        _seq_st->total_a_cnt += 1;                           break;
                     case 'C': case 'c':
                        _seq_st->total_c_cnt += 1;                           break;
                     case 'G': case 'g':
                        _seq_st->total_g_cnt += 1;                           break;
                     case 'T': case 't': case 'U': case 'u':
                        _seq_st->total_tu_cnt += 1;                           break;
                     case 'N': case 'n':
                        _seq_st->total_n_cnt += 1;                           break;
                     default:
                        std::cout<<"Error!!! from "<< thread_id << " Unknown type of nucletides: "<< br_it->qry_seq[ _si ] << " for " <<  br_it->qry_name <<std::endl;
                 }
              }
           }
           if ( !( ( (br_it->map_flag)& F5_FUNMAP ) || ((br_it->map_flag)&F5_FSECONDARY) ) ){
              match_this_read = 0;
              accuracy_perc_this_read = 0;
              for(size_t _ci=0; _ci<br_it->cigar_type.size(); _ci++){
                  switch (br_it->cigar_type[_ci]){
                     case F5_CEQUAL:
                          // ref_thread_data.t_output_F5_.num_matched_bases += br_it->cigar_len[_ci]; 
                     case F5_CMATCH: // M
                          match_this_read += br_it->cigar_len[_ci];
                          ref_thread_data.t_output_F5_.num_matched_bases += br_it->cigar_len[_ci];
                          break;
                     case F5_CINS:  // I 
                          ref_thread_data.t_output_F5_.num_ins_bases += br_it->cigar_len[_ci];
                          break;
                     case F5_CDEL:   // D
                          ref_thread_data.t_output_F5_.num_del_bases += br_it->cigar_len[_ci];
                          break;
                     case F5_CREF_SKIP: // R
                          break;
                     case F5_CSOFT_CLIP:  // S
                          ref_thread_data.t_output_F5_.num_clip_bases += br_it->cigar_len[_ci];
                          break;
                     case F5_CHARD_CLIP: // H
                          ref_thread_data.t_output_F5_.num_clip_bases += br_it->cigar_len[_ci]; 
                          break;
                     case F5_CPAD:  // P
                          break;
                     case F5_CDIFF: // X 
                          ref_thread_data.t_output_F5_.num_mismatched_bases += br_it->cigar_len[_ci];
                          break;
                     default:
                          std::cout<<"ERROR!!! from "<< thread_id << " Unknown cigar "<< br_it->cigar_type[_ci] << br_it->cigar_len[_ci] << " for " <<  br_it->qry_name <<std::endl;
                  }
              }
              accuracy_perc_this_read = int( (double(match_this_read)/br_it->qry_seq_len)*100+0.5);
              if ( accuracy_perc_this_read <0 || accuracy_perc_this_read>=PERCENTAGE_ARRAY_SIZE) {
                  std::cout<<"ERROR!!! #matched("<<match_this_read<<") is smaller than 0 or larger than total("<<br_it->qry_seq_len <<") for "<<  br_it->qry_name <<std::endl;
              }else{
                  ref_thread_data.t_output_F5_.accuracy_per_read[ int(accuracy_perc_this_read) ] += 1;
              }
           }
        }

        myMutex_output.lock();
       
        ref_output.add( ref_thread_data.t_output_F5_ );

        myMutex_output.unlock();
    }
}


