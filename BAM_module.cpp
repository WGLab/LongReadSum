#include "BAM_module.h"
#include "ComFunction.h"

BAM_Thread_data::BAM_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   br_list.reserve(p_batch_size+1);
   _thread_id = p_thread_id;
   m_input_op = ref_input_op;
}

BAM_Thread_data::~BAM_Thread_data(){ 
   ; 
   //std::cout<<" ~BAM_Thread_data " << _thread_id <<std::endl<<std::flush;
}

size_t BAM_Module::batch_size_of_record=3000;
std::mutex BAM_Module::myMutex_readBam;
std::mutex BAM_Module::myMutex_output;
size_t BAM_Module::read_i_bam = 0;

BAM_Module::BAM_Module(Input_Para& _m_input){
   m_input_op = _m_input;

   has_error = 0;
   read_i_bam = 0;

   _bam_reader_ptr = NULL;   
   if ( read_i_bam >= _m_input.num_input_files ){
      std::cout<<"Error!!! No input BAM are find. #input_files="<< _m_input.num_input_files <<std::endl;
      has_error |= 1;
      return;
   }
   
    _bam_reader_ptr = new BamReader( _m_input.input_files[read_i_bam].c_str() );
   //_bam_reader_ptr = new BamReader( _m_input.input_files[read_i_bam] );
   if (!_bam_reader_ptr->check_bam_status()){
       std::cout<< "Error!!! Cannot open bam file="<< _m_input.input_files[read_i_bam] <<std::endl;
       has_error |= 2;
   }else{
       std::cout<< "INFO: Open bam file="<< _m_input.input_files[read_i_bam] <<" " << read_i_bam<<"/"<<_m_input.num_input_files <<std::endl;
       _bam_reader_ptr->bamReadOp.set_w_supplementary(true);
       _bam_reader_ptr->bamReadOp.set_w_secondary(true);
       _bam_reader_ptr->bamReadOp.set_min_read_len(1);
       _bam_reader_ptr->bamReadOp.set_w_qry_seq(true);
       _bam_reader_ptr->bamReadOp.set_w_unmap(true);
       _bam_reader_ptr->bamReadOp.set_w_qry_qual(true);
       _bam_reader_ptr->bamReadOp.set_w_map_detail(false);
       read_i_bam = 1;
   }
}

BAM_Module::~BAM_Module(){
   //std::cout<<" ~BAM_Module begin" <<std::endl<<std::flush;
   if (_bam_reader_ptr!=NULL){ 
       delete _bam_reader_ptr; 
   }
   _bam_reader_ptr = NULL;
   //std::cout<<" ~BAM_Module end" <<std::endl<<std::flush;
}

/*
Output_BAM BAM_Module::bam_st(){
   Output_BAM t_output_bam_info;*/
int BAM_Module::bam_st( Output_BAM& t_output_bam_info){  
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   if (has_error==0){
       m_threads.reserve(m_input_op.threads+3);

      int _i_t=0;
      BAM_Thread_data** thread_data_vector = new BAM_Thread_data*[m_input_op.threads];
      try{
         for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
             std::cout<<"INFO: generate threads "<<_i_t<<std::endl<<std::flush;
             thread_data_vector[_i_t] = new BAM_Thread_data(m_input_op, _i_t);
             std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
             m_threads.push_back(std::thread((BAM_Module::BAM_do_thread), _bam_reader_ptr, std::ref(m_input_op), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(t_output_bam_info), std::ref(secondary_alignment), std::ref(supplementary_alignment) ));
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

   // std::cout<<" summ. "<<std::endl<<std::flush;
   t_output_bam_info.num_reads_with_secondary_alignment = secondary_alignment.size();
   t_output_bam_info.num_reads_with_supplementary_alignment = supplementary_alignment.size();
   // std::cout<<" summ. 2 "<<std::endl<<std::flush;
   for ( std::map<std::string, bool>::iterator sec_it=secondary_alignment.begin(),sup_it=supplementary_alignment.begin(); (sec_it!=secondary_alignment.end()) && ( sup_it!=supplementary_alignment.end()); ){
      if ( *sec_it < *sup_it ){ ++sec_it; }
      else if ( *sec_it > *sup_it ){ ++sup_it; }
      else{ t_output_bam_info.num_reads_with_both_secondary_supplementary_alignment += 1; ++sec_it; ++sup_it; }
   }
   // std::cout<<" summ. G"<<std::endl<<std::flush;
   t_output_bam_info.global_sum();
 
   auto relapse_end_time = std::chrono::high_resolution_clock::now();
   std::cout<<"Total time(replase): "<<round3((relapse_end_time-relapse_start_time).count()/1000000000.0)<<std::endl;

   std::cout<<"<Stastitics on BAM>: "<< (has_error==0?"successfully":"Failed") <<"."<<std::endl;
 
   //return t_output_bam_info; 
   return has_error;
}

void BAM_Module::BAM_do_thread(BamReader* ref_bam_reader_ptr, Input_Para& ref_input_op, int thread_id, BAM_Thread_data& ref_thread_data, Output_BAM& ref_output, std::map<std::string, bool>& ref_secondary_alignment, std::map<std::string, bool>& ref_supplementary_alignment){
    std::vector<Bam1Record>::iterator br_it;
    while (true){
        //auto start_time = std::chrono::high_resolution_clock::now();

        ref_thread_data.br_list.clear();
        myMutex_readBam.lock();
        ref_bam_reader_ptr->readBam( ref_thread_data.br_list, BAM_Module::batch_size_of_record );
        if (ref_thread_data.br_list.size() == 0 && !(read_i_bam < ref_input_op.num_input_files) ){
            std::cout<<" " << ref_thread_data.br_list.size()<<" Reads." << read_i_bam <<"/"<< ref_input_op.num_input_files <<std::endl;
            break;
        }
        if ( ref_thread_data.br_list.size() < batch_size_of_record ){
            if ( read_i_bam < ref_input_op.num_input_files ){ 
               std::cout<< "INFO: Open bam file="<< ref_input_op.input_files[read_i_bam] <<std::endl;
               ref_bam_reader_ptr->resetBam( ref_input_op.input_files[read_i_bam].c_str() );
               //ref_bam_reader_ptr->resetBam( ref_input_op.input_files[read_i_bam] );
               if ( ref_thread_data.br_list.size() == 0 ){
                  ref_bam_reader_ptr->readBam( ref_thread_data.br_list, BAM_Module::batch_size_of_record );
               }
               read_i_bam++;
            }
        }
        myMutex_readBam.unlock();
        
        ref_thread_data.t_secondary_alignment.clear();
        ref_thread_data.t_supplementary_alignment.clear();
        //std::cout<<" " << ref_thread_data.br_list.size()<<" Reads." <<std::endl; 
        for(br_it=ref_thread_data.br_list.begin(); br_it!=ref_thread_data.br_list.end(); br_it++){
           // std::cout<< br_it->qry_seq_len <<" " << br_it->qry_seq.length() << " " << br_it->qry_name << " " << int(br_it->map_flag) << std::endl <<std::flush;
           Basic_Seq_Statistics* _seq_st = NULL;
           Basic_Seq_Quality_Statistics* _seq_qual_st = NULL;
           if ( (br_it->map_flag)& BAM_FUNMAP ) { 
               // ref_thread_data.t_output_bam_.unmapped_long_read_info.total_num_reads += 1; 
               // ref_thread_data.t_output_bam_.unmapped_long_read_info.total_num_bases += br_it->qry_seq_len;
               _seq_st = &( ref_thread_data.t_output_bam_.unmapped_long_read_info);
               _seq_qual_st = &(ref_thread_data.t_output_bam_.unmapped_seq_quality_info);
           }else if ( (br_it->map_flag)& BAM_FSUPPLEMENTARY ) { 
               ref_thread_data.t_output_bam_.num_supplementary_alignment += 1; 
               ref_thread_data.t_supplementary_alignment[ br_it->qry_name ] = true;
           }else if ( (br_it->map_flag)& BAM_FSECONDARY ) { 
               ref_thread_data.t_output_bam_.num_secondary_alignment += 1;
               ref_thread_data.t_secondary_alignment[ br_it->qry_name ] = true; 
           } else {
               ref_thread_data.t_output_bam_.num_primary_alignment += 1;
               // ref_thread_data.t_output_bam_.mapped_long_read_info.total_num_reads += 1;
               // ref_thread_data.t_output_bam_.mapped_long_read_info.total_num_bases += br_it->qry_seq_len;
               _seq_st = &( ref_thread_data.t_output_bam_.mapped_long_read_info);
               _seq_qual_st = &(ref_thread_data.t_output_bam_.mapped_seq_quality_info);
           }
           
           if ( !( (br_it->map_flag)& BAM_FUNMAP ) ){
              if (br_it->map_strand==0){ ref_thread_data.t_output_bam_.forward_alignment += 1; }
              else { ref_thread_data.t_output_bam_.reverse_alignment += 1; }
           }

           if ( br_it->map_qual < ref_thread_data.t_output_bam_.min_map_quality || ref_thread_data.t_output_bam_.min_map_quality==MoneDefault){ 
               ref_thread_data.t_output_bam_.min_map_quality = br_it->map_qual; }
           if ( br_it->map_qual > ref_thread_data.t_output_bam_.max_map_quality ){ ref_thread_data.t_output_bam_.max_map_quality = br_it->map_qual; }
           ref_thread_data.t_output_bam_.map_quality_distribution[ int(br_it->map_qual) ] += 1;
           
           if (_seq_st!= NULL && !( (br_it->map_flag)& BAM_FSECONDARY ) && !((br_it->map_flag)& BAM_FSUPPLEMENTARY ) ){
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
           if ( !( ( (br_it->map_flag)& BAM_FUNMAP ) || ((br_it->map_flag)&BAM_FSECONDARY) ) ){
              for(size_t _ci=0; _ci<br_it->cigar_type.size(); _ci++){
                  switch (br_it->cigar_type[_ci]){
                     case BAM_CEQUAL:
                          // ref_thread_data.t_output_bam_.num_matched_bases += br_it->cigar_len[_ci]; 
                     case BAM_CMATCH: // M
                          ref_thread_data.t_output_bam_.num_matched_bases += br_it->cigar_len[_ci];
                          break;
                     case BAM_CINS:  // I 
                          ref_thread_data.t_output_bam_.num_ins_bases += br_it->cigar_len[_ci];
                          break;
                     case BAM_CDEL:   // D
                          ref_thread_data.t_output_bam_.num_del_bases += br_it->cigar_len[_ci];
                          break;
                     case BAM_CREF_SKIP: // R
                          break;
                     case BAM_CSOFT_CLIP:  // S
                          ref_thread_data.t_output_bam_.num_clip_bases += br_it->cigar_len[_ci];
                          break;
                     case BAM_CHARD_CLIP: // H
                          ref_thread_data.t_output_bam_.num_clip_bases += br_it->cigar_len[_ci]; 
                          break;
                     case BAM_CPAD:  // P
                          break;
                     case BAM_CDIFF: // X 
                          ref_thread_data.t_output_bam_.num_mismatched_bases += br_it->cigar_len[_ci];
                          break;
                     default:
                          std::cout<<"ERROR!!! from "<< thread_id << " Unknown cigar "<< br_it->cigar_type[_ci] << br_it->cigar_len[_ci] << " for " <<  br_it->qry_name <<std::endl;
                  }
              }
           }
        }

        myMutex_output.lock();
       
        ref_secondary_alignment.insert(ref_thread_data.t_secondary_alignment.begin(), ref_thread_data.t_secondary_alignment.end());
        ref_supplementary_alignment.insert(ref_thread_data.t_supplementary_alignment.begin(), ref_thread_data.t_supplementary_alignment.end());
        // std::cout<< ref_secondary_alignment.size()<<std::endl;
        // std::cout<< ref_supplementary_alignment.size()<<std::endl;
        ref_output.add( ref_thread_data.t_output_bam_ );

        myMutex_output.unlock();
    }
}


