/*
BAM_module.cpp:
Class for generating BAM file statistics. Records are accessed using multi-threading.
*/

#include <iostream>

#include "BAM_module.h"
#include "ComFunction.h"

BAM_Thread_data::BAM_Thread_data(Input_Para& ref_input_op, int p_thread_id, int p_batch_size=1000){
   br_list.reserve(p_batch_size+1);
   _thread_id = p_thread_id;
   m_input_op = ref_input_op;
}

BAM_Thread_data::~BAM_Thread_data(){
}

size_t BAM_Module::batch_size_of_record=3000;  // Number of records for each thread
std::mutex BAM_Module::myMutex_readBam;  // Mutex allows you to compile data across threads
std::mutex BAM_Module::myMutex_output;
size_t BAM_Module::read_i_bam = 0;

BAM_Module::BAM_Module(Input_Para& _m_input){
   m_input_op = _m_input;

   exit_code = 0;
   read_i_bam = 0;

   _bam_reader_ptr = NULL;   
   if ( read_i_bam >= _m_input.num_input_files ){
      std::cerr << "No input BAM are find. #input_files="<< _m_input.num_input_files <<std::endl;
      exit_code = 2;
   } else {
        _bam_reader_ptr = new BamReader( _m_input.input_files[read_i_bam].c_str() );
       if (!_bam_reader_ptr->check_bam_status()){
           std::cerr << "Cannot open bam file="<< _m_input.input_files[read_i_bam] <<std::endl;
           exit_code = 3;
       } else {
           std::cout << "INFO: Open bam file="<< _m_input.input_files[read_i_bam] <<" " << read_i_bam<<"/"<<_m_input.num_input_files <<std::endl;
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
}

BAM_Module::~BAM_Module(){
   if (_bam_reader_ptr!=NULL){
       delete _bam_reader_ptr; 
   }
   _bam_reader_ptr = NULL;
}

int BAM_Module::calculateStatistics( Output_BAM& t_output_bam_info){
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   if (exit_code==0){
        m_threads.reserve(m_input_op.threads+3);

        t_output_bam_info.unmapped_long_read_info.resize();
        t_output_bam_info.mapped_long_read_info.resize();
        t_output_bam_info.long_read_info.resize();

        int _i_t=0;
        BAM_Thread_data** thread_data_vector = new BAM_Thread_data*[m_input_op.threads];
        try{

            for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
                std::cout<<"INFO: generate threads "<<_i_t<<std::endl<<std::flush;
                thread_data_vector[_i_t] = new BAM_Thread_data(m_input_op, _i_t, BAM_Module::batch_size_of_record);
                std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;

                // Push the next N records onto a thread for statistics computations
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
     
      for (_i_t=0; _i_t<m_input_op.threads; _i_t++){
         delete thread_data_vector[_i_t];
      }
      delete [] thread_data_vector;

      // Compile statistics
       t_output_bam_info.num_reads_with_secondary_alignment = secondary_alignment.size();
       t_output_bam_info.num_reads_with_supplementary_alignment = supplementary_alignment.size();
       for ( std::map<std::string, bool>::iterator sec_it=secondary_alignment.begin(),sup_it=supplementary_alignment.begin(); (sec_it!=secondary_alignment.end()) && ( sup_it!=supplementary_alignment.end()); ){
          if ( *sec_it < *sup_it ){ ++sec_it; }
          else if ( *sec_it > *sup_it ){ ++sup_it; }
          else{ t_output_bam_info.num_reads_with_both_secondary_supplementary_alignment += 1; ++sec_it; ++sup_it; }
       }

       // Calculate NXX scores and GC content
       t_output_bam_info.global_sum();

       auto relapse_end_time = std::chrono::high_resolution_clock::now();
       std::cout<<"Total time(Elapsed): "<<round3((relapse_end_time-relapse_start_time).count()/1000000000.0)<<std::endl;
   }
   std::cout<<"<Statistics on BAM>: "<< (exit_code==0?"Completed successfully":"Failed") <<std::endl;

   return exit_code;
}

// Calculate statistics computations on a thread for the next N set of records from the BAM file
void BAM_Module::BAM_do_thread(BamReader* ref_bam_reader_ptr, Input_Para& ref_input_op, int thread_id, BAM_Thread_data& ref_thread_data, Output_BAM& ref_output, std::map<std::string, bool>& ref_secondary_alignment, std::map<std::string, bool>& ref_supplementary_alignment){
    std::vector<Bam1Record>::iterator record_data;
    uint64_t match_this_read;
    double accuracy_perc_this_read;
    uint64_t _t_a, _t_c, _t_g, _t_tu;
    double _t_gc_per;
    while (true){
        ref_thread_data.br_list.clear();
        myMutex_readBam.lock();

        // Read the record into the thread pointer
        ref_bam_reader_ptr->readBam( ref_thread_data.br_list, BAM_Module::batch_size_of_record );
        if (ref_thread_data.br_list.size() == 0 && !(read_i_bam < ref_input_op.num_input_files) ){
            myMutex_readBam.unlock();
            break;
        }
        if ( ref_thread_data.br_list.size() < batch_size_of_record ){
            if ( read_i_bam < ref_input_op.num_input_files ){ 
               std::cout<< "INFO: Open bam file="<< ref_input_op.input_files[read_i_bam] <<std::endl;
               ref_bam_reader_ptr->resetBam( ref_input_op.input_files[read_i_bam].c_str() );
               if ( ref_thread_data.br_list.size() == 0 ){
                  ref_bam_reader_ptr->readBam( ref_thread_data.br_list, BAM_Module::batch_size_of_record );
               }
               read_i_bam++;
            }
        }
        myMutex_readBam.unlock();

        if (ref_thread_data.br_list.size() == 0 ) { continue; }       

        // Count individual bases (mapped and unmapped) in BAM file records
        ref_thread_data.t_secondary_alignment.clear();
        ref_thread_data.t_supplementary_alignment.clear();
        ref_thread_data.t_output_bam_.reset();
        ref_thread_data.t_output_bam_.unmapped_long_read_info.resize();
        ref_thread_data.t_output_bam_.mapped_long_read_info.resize();
        ref_thread_data.t_output_bam_.long_read_info.resize();
        for(record_data=ref_thread_data.br_list.begin(); record_data!=ref_thread_data.br_list.end(); record_data++){
            Basic_Seq_Statistics* output_qc = NULL;  // Output base QC data for either mapped or unmapped reads
            Basic_Seq_Quality_Statistics* output_quality_qc = NULL;  // Output base quality QC data

            // The map flag indicates whether this record refers to unmapped or mapped data.
            //
            if ( (record_data->map_flag)& BAM_FUNMAP ) {
                // Unmapped data
                output_qc = &( ref_thread_data.t_output_bam_.unmapped_long_read_info);
                output_quality_qc = &(ref_thread_data.t_output_bam_.unmapped_seq_quality_info);

                // Update the read lengths
                int read_length = record_data->qry_seq_len;
                output_qc->read_lengths.push_back(read_length);

            }else if ( (record_data->map_flag)& BAM_FSUPPLEMENTARY ) {
                // Supplementary alignments (Only reads are counted for these)
                ref_thread_data.t_output_bam_.num_supplementary_alignment += 1;
                ref_thread_data.t_supplementary_alignment[ record_data->qry_name ] = true;

            }else if ( (record_data->map_flag)& BAM_FSECONDARY ) {
                // Secondary alignments (Only reads are counted for these)
                ref_thread_data.t_output_bam_.num_secondary_alignment += 1;
                ref_thread_data.t_secondary_alignment[ record_data->qry_name ] = true;

            } else if (! ((record_data->map_flag & BAM_FSECONDARY) || (record_data->map_flag & BAM_FSUPPLEMENTARY)) ) {
                // To get primary alignments, remove supplementary (Bit #12) and non-primary alignments (Bit #9)
                ref_thread_data.t_output_bam_.num_primary_alignment += 1;
                output_qc = &( ref_thread_data.t_output_bam_.mapped_long_read_info);
                output_quality_qc = &(ref_thread_data.t_output_bam_.mapped_seq_quality_info);

                // Update the read lengths
                int read_length = record_data->qry_seq_len;
                output_qc->read_lengths.push_back(read_length);
            }

            if ( !( (record_data->map_flag)& BAM_FUNMAP ) ){
              if (record_data->map_strand==0){ ref_thread_data.t_output_bam_.forward_alignment += 1; }
              else { ref_thread_data.t_output_bam_.reverse_alignment += 1; }
            }

            if ( record_data->map_qual < ref_thread_data.t_output_bam_.min_map_quality || ref_thread_data.t_output_bam_.min_map_quality==MoneDefault){
               ref_thread_data.t_output_bam_.min_map_quality = record_data->map_qual; }
            if ( record_data->map_qual > ref_thread_data.t_output_bam_.max_map_quality ){ ref_thread_data.t_output_bam_.max_map_quality = record_data->map_qual; }
            ref_thread_data.t_output_bam_.map_quality_distribution[ int(record_data->map_qual) ] += 1;
           
            if (output_qc!= NULL && !( (record_data->map_flag)& BAM_FSECONDARY ) && !((record_data->map_flag)& BAM_FSUPPLEMENTARY ) ){
                output_qc->total_num_reads += 1;
                output_qc->total_num_bases += record_data->qry_seq_len;

                if( record_data->qry_seq_len > output_qc->longest_read_length){ output_qc->longest_read_length = record_data->qry_seq_len; }
                if (  record_data->qry_seq_len < MAX_READ_LENGTH ){
                    output_qc->read_length_count[record_data->qry_seq_len] += 1;
                }

                _t_a=0;
                _t_c=0;
                _t_g=0;
                _t_tu=0;
                _t_gc_per = 0;

                for (size_t _si=0; _si<record_data->qry_seq_len; _si++){
                    int _base_qual_int = int(record_data->qry_qual[_si]);
                    if ( _base_qual_int > MAX_BASE_QUALITY || _base_qual_int <0){
                    std::cout<<"WARNINING!!! from "<< thread_id << " The base quality is not in [0,"<< MAX_BASE_QUALITY <<"]: " << int(record_data->qry_qual[_si]) << " at " << _si << " for " <<  record_data->qry_name  <<std::endl;
                    }else{
                        output_quality_qc->base_quality_distribution[_base_qual_int ] += 1;
                        if ( output_quality_qc->min_base_quality==MoneDefault || _base_qual_int < output_quality_qc->min_base_quality){
                        output_quality_qc->min_base_quality = _base_qual_int; }
                        if ( _base_qual_int > output_quality_qc->max_base_quality) { output_quality_qc->max_base_quality = _base_qual_int; }
                    }

                    switch( record_data->qry_seq[ _si ] ){
                        case 'A': case 'a':
                            output_qc->total_a_cnt += 1;    _t_a += 1;                       break;
                        case 'C': case 'c':
                            output_qc->total_c_cnt += 1;    _t_c += 1;                       break;
                        case 'G': case 'g':
                            output_qc->total_g_cnt += 1;    _t_g += 1;                       break;
                        case 'T': case 't': case 'U': case 'u':
                            output_qc->total_tu_cnt += 1;   _t_tu += 1;                        break;
                        case 'N': case 'n':
                            output_qc->total_n_cnt += 1;                           break;
                        default:
                            std::cout<<"Error!!! from "<< thread_id << " Unknown type of nucletides: "<< record_data->qry_seq[ _si ] << " for " <<  record_data->qry_name <<std::endl;
                    }
                }
                _t_gc_per = ( _t_c + _t_g)/( (record_data->qry_seq_len>0?record_data->qry_seq_len:-1)/100.0 );
                if ( record_data->qry_seq_len>0 ){
                    if ( int(_t_gc_per+0.5) < PERCENTAGE_ARRAY_SIZE){
                     output_qc->read_gc_content_count[ int(_t_gc_per+0.5) ] += 1;
                    }else{
                     std::cout<<"ERROR!!! GC (%) content("<<_t_gc_per<<") is larger than total("<<PERCENTAGE_ARRAY_SIZE <<") for "<<  record_data->qry_name <<std::endl;
                    }
                }
            }
            if ( !( ( (record_data->map_flag)& BAM_FUNMAP ) || ((record_data->map_flag)&BAM_FSECONDARY) ) ){
                match_this_read = 0;
                accuracy_perc_this_read = 0;
                for(size_t _ci=0; _ci<record_data->cigar_type.size(); _ci++){
                    switch (record_data->cigar_type[_ci]){
                         case BAM_CEQUAL:
                         case BAM_CMATCH: // M
                              match_this_read += record_data->cigar_len[_ci];
                              ref_thread_data.t_output_bam_.num_matched_bases += record_data->cigar_len[_ci];
                              break;
                         case BAM_CINS:  // I
                              ref_thread_data.t_output_bam_.num_ins_bases += record_data->cigar_len[_ci];
                              break;
                         case BAM_CDEL:   // D
                              ref_thread_data.t_output_bam_.num_del_bases += record_data->cigar_len[_ci];
                              break;
                         case BAM_CREF_SKIP: // R
                              break;
                         case BAM_CSOFT_CLIP:  // S
                              ref_thread_data.t_output_bam_.num_clip_bases += record_data->cigar_len[_ci];
                              break;
                         case BAM_CHARD_CLIP: // H
                              ref_thread_data.t_output_bam_.num_clip_bases += record_data->cigar_len[_ci];
                              break;
                         case BAM_CPAD:  // P
                              break;
                         case BAM_CDIFF: // X
                              ref_thread_data.t_output_bam_.num_mismatched_bases += record_data->cigar_len[_ci];
                              break;
                         default:
                              std::cout<<"ERROR!!! from "<< thread_id << " Unknown cigar "<< record_data->cigar_type[_ci] << record_data->cigar_len[_ci] << " for " <<  record_data->qry_name <<std::endl;
                    }
                }
                accuracy_perc_this_read = int( (double(match_this_read)/record_data->qry_seq_len)*100+0.5);
                if ( accuracy_perc_this_read <0 || accuracy_perc_this_read>=PERCENTAGE_ARRAY_SIZE) {
                    std::cout<<"ERROR!!! #matched("<<match_this_read<<") is smaller than 0 or larger than total("<<record_data->qry_seq_len <<") for "<<  record_data->qry_name <<std::endl;
                } else {
                    ref_thread_data.t_output_bam_.accuracy_per_read[ int(accuracy_perc_this_read) ] += 1;
                }
            }
        }
        myMutex_output.lock();
       
        ref_secondary_alignment.insert(ref_thread_data.t_secondary_alignment.begin(), ref_thread_data.t_secondary_alignment.end());
        ref_supplementary_alignment.insert(ref_thread_data.t_supplementary_alignment.begin(), ref_thread_data.t_supplementary_alignment.end());
        ref_output.add( ref_thread_data.t_output_bam_ );

        myMutex_output.unlock();
    }
}


