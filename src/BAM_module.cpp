/*
BAM_module.cpp:
Class for generating BAM file statistics. Records are accessed using multi-threading.
*/

#include <iostream>
#include <thread>
#include <iostream>
#include <cmath>

#include "BAM_module.h"
#include "ComFunction.h"

BAM_Thread_data::BAM_Thread_data(Input_Para& input_params, int p_thread_id, int p_batch_size=1000){
//   this->record_list.reserve(p_batch_size+1);
   _thread_id = p_thread_id;
   m_input_op = input_params;
}

BAM_Thread_data::~BAM_Thread_data(){
}


// Calculate statistics computations on a thread for the next N set of records from the BAM file
void BAM_Module::BAM_do_thread(BamReader& bam_reader, Input_Para& input_params, int thread_id, BAM_Thread_data& ref_thread_data, Output_BAM& ref_output, std::map<std::string, bool>& ref_secondary_alignment, std::map<std::string, bool>& ref_supplementary_alignment, std::mutex& bam_mutex, std::mutex& output_mutex, size_t &file_index){
    std::vector<Bam1Record>::iterator record_data;
    uint64_t match_this_read;
    double accuracy_perc_this_read;
    uint64_t _t_a, _t_c, _t_g, _t_tu;
    double _t_gc_per;
//    while (true){
//        ref_thread_data.record_list.clear();
//
//        // Lock the mutex before reading the records
//        bam_mutex.lock();
//
//        // Read the record into the thread pointer
//        size_t file_count = input_params.num_input_files;
//        size_t records_per_thread = BAM_Module::records_per_thread;
//
//        // Read a batch of records in the BAM file into the record list
//        refbam_reader->readBam( ref_thread_data.record_list, records_per_thread );
//
//
//
//        // If the record list is empty and there are no more files to read, then break out of the loop
//        if (ref_thread_data.record_list.size() == 0 && !(file_index < file_count) ){
//
//            // Unlock the mutex after reading the records
//            bam_mutex.unlock();
//            break;
//        }
//
//        // Print the record list size and records per thread
//        std::cout << "Thread " << thread_id << " record list size: " << ref_thread_data.record_list.size() << " records per thread: " << records_per_thread << std::endl;
//
//        // Check if the record list is less than the batch size
//        if ( ref_thread_data.record_list.size() < records_per_thread ){
//            std::cout << "File index: " << file_index << " File count: " << file_count << std::endl;
//            std::cout << "Completed processing file: " << input_params.input_files[file_index] << std::endl;
//            file_index += 1;
//
//            // Read a new file if the current file is exhausted of records
//            if ( file_index < file_count ){
//                std::cout<< "Open BAM file: " << input_params.input_files[file_index] <<std::endl;
//                refbam_reader->resetBam( input_params.input_files[file_index].c_str() );
//                if ( ref_thread_data.record_list.size() == 0 ){
//                    refbam_reader->readBam( ref_thread_data.record_list, records_per_thread );
//
//                    // Continue if there are no records for this file
//                    if (ref_thread_data.record_list.size() == 0 )
//                    {
//                        std::cerr << "No records in file: " << input_params.input_files[file_index] << std::endl;
//
//                        // Unlock the mutex after reading the records
//                        bam_mutex.unlock();
//                        continue;
//                    }
//                }
//            }
//        }
//
//        // Unlock the mutex after reading the records
//        bam_mutex.unlock();
//
//        // Count individual bases (mapped and unmapped) in BAM file records
//        ref_thread_data.t_secondary_alignment.clear();
//        ref_thread_data.t_supplementary_alignment.clear();
//        ref_thread_data.t_output_bam_.reset();
//        ref_thread_data.t_output_bam_.unmapped_long_read_info.resize();
//        ref_thread_data.t_output_bam_.mapped_long_read_info.resize();
//        ref_thread_data.t_output_bam_.long_read_info.resize();
//        for(record_data=ref_thread_data.record_list.begin(); record_data!=ref_thread_data.record_list.end(); record_data++){
//            Basic_Seq_Statistics* output_qc = NULL;  // Output base QC data for either mapped or unmapped reads
//            Basic_Seq_Quality_Statistics* output_quality_qc = NULL;  // Output base quality QC data
//
//            // Set the number of mismatches using the NM tag
//            ref_thread_data.t_output_bam_.num_mismatched_bases = record_data->num_mismatch;
//
//            // The map flag indicates whether this record refers to unmapped or mapped data.
//            if ( (record_data->map_flag)& BAM_FUNMAP ) {
//                // Unmapped data
//                output_qc = &( ref_thread_data.t_output_bam_.unmapped_long_read_info);
//                output_quality_qc = &(ref_thread_data.t_output_bam_.unmapped_seq_quality_info);
//
//                // Update the read lengths
//                int read_length = record_data->qry_seq_len;
//                output_qc->read_lengths.push_back(read_length);
//
//            }else if ( (record_data->map_flag)& BAM_FSUPPLEMENTARY ) {
//                // Supplementary alignments (Only reads are counted for these)
//                ref_thread_data.t_output_bam_.num_supplementary_alignment += 1;
//                ref_thread_data.t_supplementary_alignment[ record_data->qry_name ] = true;
//
//            }else if ( (record_data->map_flag)& BAM_FSECONDARY ) {
//                // Secondary alignments (Only reads are counted for these)
//                ref_thread_data.t_output_bam_.num_secondary_alignment += 1;
//                ref_thread_data.t_secondary_alignment[ record_data->qry_name ] = true;
//
//            } else if (! ((record_data->map_flag & BAM_FSECONDARY) || (record_data->map_flag & BAM_FSUPPLEMENTARY)) ) {
//                // To get primary alignments, remove supplementary (Bit #12) and non-primary alignments (Bit #9)
//                ref_thread_data.t_output_bam_.num_primary_alignment += 1;
//                output_qc = &( ref_thread_data.t_output_bam_.mapped_long_read_info);
//                output_quality_qc = &(ref_thread_data.t_output_bam_.mapped_seq_quality_info);
//
//                // Update the read lengths
//                int read_length = record_data->qry_seq_len;
//                output_qc->read_lengths.push_back(read_length);
//            }
//
//            if ( !( (record_data->map_flag)& BAM_FUNMAP ) ){
//              if (record_data->map_strand==0){ ref_thread_data.t_output_bam_.forward_alignment += 1; }
//              else { ref_thread_data.t_output_bam_.reverse_alignment += 1; }
//            }
//
//            if ( record_data->map_qual < ref_thread_data.t_output_bam_.min_map_quality || ref_thread_data.t_output_bam_.min_map_quality==MoneDefault){
//               ref_thread_data.t_output_bam_.min_map_quality = record_data->map_qual; }
//            if ( record_data->map_qual > ref_thread_data.t_output_bam_.max_map_quality ){ ref_thread_data.t_output_bam_.max_map_quality = record_data->map_qual; }
//            ref_thread_data.t_output_bam_.map_quality_distribution[ int(record_data->map_qual) ] += 1;
//
//            if (output_qc!= NULL && !( (record_data->map_flag)& BAM_FSECONDARY ) && !((record_data->map_flag)& BAM_FSUPPLEMENTARY ) ){
//                output_qc->total_num_reads += 1;
//                output_qc->total_num_bases += record_data->qry_seq_len;
//
//                if( record_data->qry_seq_len > output_qc->longest_read_length){ output_qc->longest_read_length = record_data->qry_seq_len; }
//                if (  record_data->qry_seq_len < MAX_READ_LENGTH ){
//                    output_qc->read_length_count[record_data->qry_seq_len] += 1;
//                }
//
//                _t_a=0;
//                _t_c=0;
//                _t_g=0;
//                _t_tu=0;
//                _t_gc_per = 0;
//
//                for (size_t _si=0; _si<record_data->qry_seq_len; _si++){
//                    int _base_qual_int = int(record_data->qry_qual[_si]);
//                    if ( _base_qual_int > MAX_BASE_QUALITY || _base_qual_int <0){
//                    std::cerr << "Thread "<< thread_id << ": The base quality is not in [0,"<< MAX_BASE_QUALITY <<"]: " << int(record_data->qry_qual[_si]) << " at " << _si << " for " <<  record_data->qry_name  <<std::endl;
//                    }else{
//                        output_quality_qc->base_quality_distribution[_base_qual_int ] += 1;
//                        if ( output_quality_qc->min_base_quality==MoneDefault || _base_qual_int < output_quality_qc->min_base_quality){
//                        output_quality_qc->min_base_quality = _base_qual_int; }
//                        if ( _base_qual_int > output_quality_qc->max_base_quality) { output_quality_qc->max_base_quality = _base_qual_int; }
//                    }
//
//                    switch( record_data->qry_seq[ _si ] ){
//                        case 'A': case 'a':
//                            output_qc->total_a_cnt += 1;    _t_a += 1;                       break;
//                        case 'C': case 'c':
//                            output_qc->total_c_cnt += 1;    _t_c += 1;                       break;
//                        case 'G': case 'g':
//                            output_qc->total_g_cnt += 1;    _t_g += 1;                       break;
//                        case 'T': case 't': case 'U': case 'u':
//                            output_qc->total_tu_cnt += 1;   _t_tu += 1;                        break;
//                        case 'N': case 'n':
//                            output_qc->total_n_cnt += 1;                           break;
//                        default:
//                            std::cerr<<"Thread "<< thread_id << ": Unknown type of nucletides: "<< record_data->qry_seq[ _si ] << " for " <<  record_data->qry_name <<std::endl;
//                    }
//                }
//                _t_gc_per = ( _t_c + _t_g)/( (record_data->qry_seq_len>0?record_data->qry_seq_len:-1)/100.0 );
//                if ( record_data->qry_seq_len>0 ){
//                    if ( int(_t_gc_per+0.5) < PERCENTAGE_ARRAY_SIZE){
//                     output_qc->read_gc_content_count[ int(_t_gc_per+0.5) ] += 1;
//                    }else{
//                     std::cerr<<"GC (%) content("<<_t_gc_per<<") is larger than total("<<PERCENTAGE_ARRAY_SIZE <<") for "<<  record_data->qry_name <<std::endl;
//                    }
//                }
//            }
//            if ( !( ( (record_data->map_flag)& BAM_FUNMAP ) || ((record_data->map_flag)&BAM_FSECONDARY) ) ){
//                match_this_read = 0;
//                accuracy_perc_this_read = 0;
//                for(size_t _ci=0; _ci<record_data->cigar_type.size(); _ci++){
//                    switch (record_data->cigar_type[_ci]){
//                         case BAM_CEQUAL:
//                         case BAM_CMATCH: // M
//                              match_this_read += record_data->cigar_len[_ci];
//                              ref_thread_data.t_output_bam_.num_matched_bases += record_data->cigar_len[_ci];
//
//                              // Update the number of columns for the percent identity computation
//                              ref_thread_data.t_output_bam_.num_columns += record_data->cigar_len[_ci];
//                              break;
//                         case BAM_CINS:  // I
//                              ref_thread_data.t_output_bam_.num_ins_bases += record_data->cigar_len[_ci];
//
//                              // Update the number of columns for the percent identity computation
//                              ref_thread_data.t_output_bam_.num_columns += record_data->cigar_len[_ci];
//                              break;
//                         case BAM_CDEL:   // D
//                              ref_thread_data.t_output_bam_.num_del_bases += record_data->cigar_len[_ci];
//
//                              // Update the number of columns for the percent identity computation
//                              ref_thread_data.t_output_bam_.num_columns += record_data->cigar_len[_ci];
//                              break;
//                         case BAM_CREF_SKIP: // R
//                              break;
//                         case BAM_CSOFT_CLIP:  // S
//                              ref_thread_data.t_output_bam_.num_clip_bases += record_data->cigar_len[_ci];
//                              break;
//                         case BAM_CHARD_CLIP: // H
//                              ref_thread_data.t_output_bam_.num_clip_bases += record_data->cigar_len[_ci];
//                              break;
//                         case BAM_CPAD:  // P
//                              break;
//                         case BAM_CDIFF: // X
//                              // Below is replaced with the NM tag in BAMReader.cpp
//                              //ref_thread_data.t_output_bam_.num_mismatched_bases += record_data->cigar_len[_ci];
//                              break;
//                         default:
//                              std::cerr<<"Thread "<< thread_id << ": Unknown cigar "<< record_data->cigar_type[_ci] << record_data->cigar_len[_ci] << " for " <<  record_data->qry_name <<std::endl;
//                    }
//                }
//
//                // Calculate percent identity
//                // = Number of matching bases / number of alignment columns
//                // = (num columns - NM) / num columns
//                // Get the number of columns
//                double num_columns_ = ref_thread_data.t_output_bam_.num_columns;
//                double num_mismatch_ = record_data->num_mismatch;
//                ref_thread_data.t_output_bam_.percent_identity = ((num_columns_ - num_mismatch_) / num_columns_) * 100.0;
//                //std::cout << "Percent identity = " << ref_thread_data.t_output_bam_.percent_identity << std::endl;
//
//                // Calculate accuracy percentage if non-zero matches
//                if (match_this_read > 0) {
//                    accuracy_perc_this_read = int( (double(match_this_read)/record_data->qry_seq_len) * 100+0.5);
//                } else {
//                    accuracy_perc_this_read = 0;
//                }
//                ref_thread_data.t_output_bam_.accuracy_per_read[ int(accuracy_perc_this_read) ] += 1;
//            }
//        }
//        output_mutex.lock();
//
//        ref_secondary_alignment.insert(ref_thread_data.t_secondary_alignment.begin(), ref_thread_data.t_secondary_alignment.end());
//        ref_supplementary_alignment.insert(ref_thread_data.t_supplementary_alignment.begin(), ref_thread_data.t_supplementary_alignment.end());
//        ref_output.add( ref_thread_data.t_output_bam_ );
//
//        output_mutex.unlock();
//    }
}


int BAM_Module::calculateStatistics( Output_BAM& t_output_bam_info){
   auto relapse_start_time = std::chrono::high_resolution_clock::now();

   if (exit_code==0){
        m_threads.reserve(m_input_op.threads+3);

        t_output_bam_info.unmapped_long_read_info.resize();
        t_output_bam_info.mapped_long_read_info.resize();
        t_output_bam_info.long_read_info.resize();

//        BAM_Thread_data** thread_data_vector = new BAM_Thread_data*[m_input_op.threads];


        // Loop through the input files
        int file_count = (int)m_input_op.num_input_files;
        for (int i=0; i < file_count; i++){
            this->file_index = i;

            // Print the number of records
            std::string filepath(m_input_op.input_files[this->file_index]);
            int record_count = bam_reader->getNumberOfRecords(filepath.c_str());
            std::cout << "INFO: Number of records = " << record_count << std::endl;

            std::cout<<"INFO: Processing file "<< filepath << std::endl;

            // Determine the number of records per thread
            int thread_count = (double)m_input_op.threads;
            int batch_size = std::ceil((double)record_count / (double)thread_count);
            std::cout << "INFO: Number of records per thread = " << batch_size << std::endl;

            // Create a vector of threads
            std::vector<std::thread> thread_vector;

            // Update the BAM file
            bam_reader->resetBam(filepath.c_str());

            // Calculate statistics in batches
            for (int thread_index=0; thread_index<thread_count; thread_index++){
                std::cout<<"INFO: generate thread "<<thread_index<<std::endl<<std::flush;
//                thread_data_vector.push_back(new BAM_Thread_data(m_input_op, _i_t, batch_size));
//                std::cout<<"INFO: Thread = "<< _i_t+1  <<std::endl<<std::flush;
                std::thread t((BAM_Module::batchStatistics), std::ref(bam_reader), batch_size, std::ref(m_input_op),std::ref(t_output_bam_info), std::ref(this->bam_mutex), std::ref(this->output_mutex));
                thread_vector.push_back(std::move(t));
                // Push the next N records onto a thread for statistics computations
//                m_threads.push_back(std::thread((BAM_Module::BAM_do_thread), bam_reader, std::ref(m_input_op), _i_t, std::ref(*(thread_data_vector[_i_t])), std::ref(t_output_bam_info), std::ref(secondary_alignment), std::ref(supplementary_alignment)));
            }

            // Join the threads in thread_vector
            for (auto& t : thread_vector){
                t.join();
                std::cout<<"INFO: Joined thread "<<std::endl<<std::flush;
            }
        }
    }
   return exit_code;
}

void BAM_Module::batchStatistics(BamReader& bam_reader, int batch_size, Input_Para& input_params, Output_BAM& ref_output, std::mutex& bam_mutex, std::mutex& output_mutex)
{

    // Create a vector of BAM records
    std::vector<Bam1Record> record_list;

    // Read the next N records in the BAM file
    record_list = bam_reader->readNextNRecords(batch_size);

    // Create the output structure
    Output_BAM thread_data;

//    // Loop through each record
//    for (auto& record : record_list){
//        // Get the record data
//        Bam1RecordData* record_data = record.getRecordData();
//
//        // Get the record alignment
//        Bam1Alignment* record_alignment = record.getAlignment();
//
//        // Get the record sequence
//        Bam1Sequence* record_sequence = record.getSequence();
//
//        // Get the record cigar
//        Bam1Cigar* record_cigar = record.getCigar();
//
//        // Get the record quality
//        Bam1Quality* record_quality = record.getQuality();
//
//        // Get the record tags
//        Bam1Tags* record_tags = record.getTags();
//
//        // Get the re
//    }
}

BAM_Module::BAM_Module(Input_Para& _m_input){
   m_input_op = _m_input;

   exit_code = 0;
   this->file_index = 0;

   bam_reader = NULL;
   if ( this->file_index >= _m_input.num_input_files ){
      std::cerr << "No input BAM are find. #input_files="<< _m_input.num_input_files <<std::endl;
      exit_code = 2;
   } else {
        // Open the BAM file
        std::string filepath(_m_input.input_files[this->file_index]);

        // Create a new reader object
        this->bam_reader = BamReader(filepath.c_str());
//        bam_reader = new BamReader(filepath.c_str() );
        if (!this->bam_reader->check_bam_status()){
            std::cerr << "Cannot open bam file="<< filepath <<std::endl;
            exit_code = 3;
        } else {
            std::cout << "Opened file " << this->file_index+1 << "/" <<_m_input.num_input_files << ": " << _m_input.input_files[this->file_index] << std::endl;
            bam_reader->bamReadOp.set_w_supplementary(true);
            bam_reader->bamReadOp.set_w_secondary(true);
            bam_reader->bamReadOp.set_min_read_len(1);
            bam_reader->bamReadOp.set_w_qry_seq(true);
            bam_reader->bamReadOp.set_w_unmap(true);
            bam_reader->bamReadOp.set_w_qry_qual(true);
            bam_reader->bamReadOp.set_w_map_detail(false);
        }
   }
}

BAM_Module::~BAM_Module(){
   if (bam_reader!=NULL){
       delete bam_reader;
   }
   bam_reader = NULL;
}
