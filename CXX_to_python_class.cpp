#include "CXX_to_python_class.h"

Output_Info::Output_Info(){
   error_flag = 0;
   //error_str = "";
}


//// function for Basic_Seq_Statistics
//
//
Basic_Seq_Statistics::Basic_Seq_Statistics(){
   read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = ZeroDefault;
   }
   nx_read_length[9] = ZeroDefault;
}

Basic_Seq_Statistics::~Basic_Seq_Statistics(){
   delete [] read_length_count;
}

Basic_Seq_Statistics::Basic_Seq_Statistics( const Basic_Seq_Statistics& _bss){
   read_length_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = _bss.read_length_count[ _i_ ];
   }
   for(int _i_=0; _i_< 10 ; _i_++){
      nx_read_length[ _i_ ] = _bss.nx_read_length[ _i_ ];
   }
   total_num_reads = _bss.total_num_reads ;
   total_num_bases = _bss.total_num_bases ;

   longest_read_length = _bss.longest_read_length ;
   n50_read_length = _bss.n50_read_length ;
   n95_read_length = _bss.n95_read_length  ;
   n05_read_length = _bss.n05_read_length ;
   mean_read_length = _bss.mean_read_length ;
   median_read_length = _bss.median_read_length ;

   total_a_cnt = _bss.total_a_cnt ;
   total_c_cnt = _bss.total_c_cnt ;
   total_g_cnt = _bss.total_g_cnt ;
   total_tu_cnt = _bss.total_tu_cnt ;
   total_n_cnt = _bss.total_n_cnt ;
   gc_cnt = _bss.gc_cnt ;
}

void Basic_Seq_Statistics::reset(){
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_< 10 ; _i_++){
      nx_read_length[ _i_ ] = ZeroDefault;
   }
   
   total_num_reads = ZeroDefault; 
   total_num_bases = ZeroDefault; 

   longest_read_length = ZeroDefault;
   n50_read_length = MoneDefault; 
   n95_read_length = MoneDefault; 
   n05_read_length = MoneDefault; 
   mean_read_length = MoneDefault; 
   median_read_length = MoneDefault; 

   total_a_cnt = ZeroDefault; 
   total_c_cnt = ZeroDefault; 
   total_g_cnt = ZeroDefault; 
   total_tu_cnt = ZeroDefault; 
   total_n_cnt = ZeroDefault; 
   gc_cnt = ZeroDefault; 
}

void Basic_Seq_Statistics::add(Basic_Seq_Statistics& t_seq_st){
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] += t_seq_st. read_length_count[ _i_ ];
   }
   /*for(int _i_=0; _i_< 10 ; _i_++){
      nx_read_length[ _i_ ] = t_seq_st. ;
   }
   n50_read_length = t_seq_st. ;
   n95_read_length = t_seq_st. ;
   n05_read_length = t_seq_st. ;
   mean_read_length = t_seq_st. ;
   median_read_length = t_seq_st. ;*/

   total_num_reads += t_seq_st.total_num_reads ;
   total_num_bases += t_seq_st.total_num_bases ;

   if ( longest_read_length < t_seq_st.longest_read_length){
      longest_read_length = t_seq_st.longest_read_length;
   }

   total_a_cnt += t_seq_st.total_a_cnt ;
   total_c_cnt += t_seq_st.total_c_cnt ;
   total_g_cnt += t_seq_st.total_g_cnt ;
   total_tu_cnt += t_seq_st.total_tu_cnt ;
   total_n_cnt += t_seq_st.total_n_cnt ;
   //gc_cnt = t_seq_st. ;
}

void Basic_Seq_Statistics::global_sum(){
   gc_cnt = float( total_c_cnt + total_g_cnt)/(total_num_bases>0?total_num_bases:1);

   mean_read_length = total_num_bases / float(total_num_reads>0?total_num_reads:1);
   uint64_t t_base=0, t_read=0; 
   bool has_cal_median = false;
   double base_perc = 0;
   int start_perct = 0;
   uint64_t t_max_length = (MAX_READ_LENGTH<longest_read_length?MAX_READ_LENGTH:longest_read_length);
   for(size_t _i_=0; _i_<t_max_length; _i_++){
      if ( read_length_count[ _i_ ] == 0) { continue; }

      t_base += _i_ * read_length_count[ _i_ ];
      t_read += read_length_count[ _i_ ];
      if ( t_read > 0){
         if ( t_read/float( total_num_reads )< 0.5 ) { median_read_length = _i_; }
         else if ( t_read/float( total_num_reads ) == 0.5 && !has_cal_median ){ median_read_length = _i_; has_cal_median = true; }
         else if ( t_read/float( total_num_reads ) > 0.5 && !has_cal_median ){
            median_read_length = ( median_read_length + _i_ )/2;
            has_cal_median = true;
         }
         
         base_perc = double( t_base )/total_num_bases;
         start_perct = int(base_perc*10+0.5);
         for (int _t_sp=start_perct; _t_sp<start_perct+1; _t_sp++){
            if ( nx_read_length[ _t_sp ]==ZeroDefault && base_perc>=_t_sp/float(10) ){
               nx_read_length[ _t_sp ] = _i_;
            }
            if ( base_perc<0.05){ n05_read_length = _i_; }
            if ( base_perc<0.5){ n50_read_length = _i_; }
            if ( base_perc<0.95){ n95_read_length = _i_; }
         }
      }
   }
}

//// function for Basic_Seq_Quality_Statistics
//
//
Basic_Seq_Quality_Statistics::Basic_Seq_Quality_Statistics(){
   pos_quality_distribution = new int[MAX_READ_LENGTH];
   pos_quality_distribution_dev = new float[MAX_READ_LENGTH];
   pos_quality_distribution_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = ZeroDefault;
      pos_quality_distribution_dev[ _i_ ] = ZeroDefault;
      pos_quality_distribution_count[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] = ZeroDefault;
   }
}

Basic_Seq_Quality_Statistics::~Basic_Seq_Quality_Statistics(){
   delete [] pos_quality_distribution;
   delete [] pos_quality_distribution_dev;
   delete [] pos_quality_distribution_count;
}

Basic_Seq_Quality_Statistics::Basic_Seq_Quality_Statistics( const Basic_Seq_Quality_Statistics& _bsqs){
   pos_quality_distribution = new int[MAX_READ_LENGTH];
   pos_quality_distribution_dev = new float[MAX_READ_LENGTH];
   pos_quality_distribution_count = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = _bsqs.pos_quality_distribution[ _i_ ];
      pos_quality_distribution_dev[ _i_ ] = _bsqs.pos_quality_distribution_dev[ _i_ ];
      pos_quality_distribution_count[ _i_ ] = _bsqs.pos_quality_distribution_count[ _i_ ];
   }
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] = _bsqs.base_quality_distribution[ _i_ ];
   }

   min_base_quality = _bsqs.min_base_quality;
   max_base_quality = _bsqs.max_base_quality;

   max_length = _bsqs.max_length;
}

void Basic_Seq_Quality_Statistics::reset(){
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = ZeroDefault;
      pos_quality_distribution_dev[ _i_ ] = ZeroDefault;
      pos_quality_distribution_count[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] = ZeroDefault;
   }

   min_base_quality = MoneDefault; 
   max_base_quality = MoneDefault; 

   max_length = ZeroDefault;
}
void Basic_Seq_Quality_Statistics::add(Basic_Seq_Quality_Statistics& t_qual_st){
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = t_qual_st.pos_quality_distribution[ _i_ ];
      pos_quality_distribution_dev[ _i_ ] = t_qual_st.pos_quality_distribution_dev[ _i_ ] ;
      pos_quality_distribution_count[ _i_ ] += t_qual_st.pos_quality_distribution_count[ _i_ ] ;
   }
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] += t_qual_st.base_quality_distribution[ _i_ ];
   }
   
   if ( min_base_quality < 0 || min_base_quality > t_qual_st.min_base_quality){
      min_base_quality = t_qual_st.min_base_quality;
   }
   if ( max_base_quality < t_qual_st.max_base_quality){
      max_base_quality = t_qual_st.max_base_quality;
   }

   if ( max_length < t_qual_st.max_length){
      max_length = t_qual_st.max_length;
   }
}

void Basic_Seq_Quality_Statistics::global_sum(){
   ;
}

//// function for Output_BAM
//
//
Output_BAM::Output_BAM(){
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] = ZeroDefault;
   }
   
}

Output_BAM::~Output_BAM(){
}

void Output_BAM::reset(){
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] = ZeroDefault;
   }

   num_primary_alignment = ZeroDefault; 
   num_secondary_alignment = ZeroDefault; 
   num_reads_with_secondary_alignment = ZeroDefault; 
   num_supplementary_alignment = ZeroDefault; 
   num_reads_with_supplementary_alignment = ZeroDefault; 
   num_reads_with_both_secondary_supplementary_alignment = ZeroDefault; 

   forward_alignment = ZeroDefault;
   reverse_alignment = ZeroDefault;

   min_map_quality = MoneDefault; 
   max_map_quality = MoneDefault; 

   num_matched_bases = ZeroDefault; 
   num_mismatched_bases = ZeroDefault; 
   num_ins_bases = ZeroDefault; 
   num_del_bases = ZeroDefault; 
   num_clip_bases = ZeroDefault; 

   mapped_long_read_info.reset();
   unmapped_long_read_info.reset();
   mapped_seq_quality_info.reset();
   unmapped_seq_quality_info.reset();
}

void Output_BAM::add(Output_BAM& t_output_bam){
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] += t_output_bam.map_quality_distribution[ _i_ ];
   }

   num_primary_alignment += t_output_bam.num_primary_alignment ;
   num_secondary_alignment += t_output_bam.num_secondary_alignment ;
   num_reads_with_secondary_alignment += t_output_bam.num_reads_with_secondary_alignment ;
   num_supplementary_alignment += t_output_bam.num_supplementary_alignment;
   num_reads_with_supplementary_alignment += t_output_bam.num_reads_with_supplementary_alignment;
   num_reads_with_both_secondary_supplementary_alignment += t_output_bam.num_reads_with_both_secondary_supplementary_alignment;

   forward_alignment += t_output_bam.forward_alignment;
   reverse_alignment += t_output_bam.reverse_alignment;

   if ( min_map_quality < 0 || min_map_quality > t_output_bam.min_map_quality){
      min_map_quality = t_output_bam.min_map_quality;
   }
   if ( max_map_quality < t_output_bam.max_map_quality ){
      max_map_quality =  t_output_bam.max_map_quality;
   }

   num_matched_bases += t_output_bam.num_matched_bases;
   num_mismatched_bases += t_output_bam.num_mismatched_bases;
   num_ins_bases += t_output_bam.num_ins_bases;
   num_del_bases += t_output_bam.num_del_bases;
   num_clip_bases += t_output_bam.num_clip_bases;
   
   mapped_long_read_info.add(t_output_bam.mapped_long_read_info);
   unmapped_long_read_info.add(t_output_bam.unmapped_long_read_info);
   mapped_seq_quality_info.add(t_output_bam.mapped_seq_quality_info);
   unmapped_seq_quality_info.add(t_output_bam.unmapped_seq_quality_info);
}

void Output_BAM::global_sum(){
   mapped_long_read_info.global_sum();
   unmapped_long_read_info.global_sum();
   mapped_seq_quality_info.global_sum();
   unmapped_seq_quality_info.global_sum();
}

//// function for Output_F5
//
//
//
Basic_F5_Statistics::Basic_F5_Statistics(){
   /*read_length_list = new int64_t[MAX_READ_LENGTH];
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_list[ _i_ ] = ZeroDefault;
   }*/
   for(int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] = ZeroDefault;
   }
}

Basic_F5_Statistics::~Basic_F5_Statistics(){
   // delete [] read_length_list;
}


