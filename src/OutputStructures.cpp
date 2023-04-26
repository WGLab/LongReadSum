#include <numeric>  // std::accumulate
#include<algorithm>  // std::foreach
#include <math.h>  // sqrt
#include <iostream>

#include "OutputStructures.h"
#include "BasicStatistics.h"

// Base class for storing error output.
Output_Info::Output_Info(){
   error_flag = 0;
   error_str = "";
}


// Base class for storing basic QC data
Basic_Seq_Statistics::Basic_Seq_Statistics(){
   read_length_count.resize(MAX_READ_LENGTH);
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = ZeroDefault;
   }
   nx_read_length.resize(10);
   for(int _i_=0; _i_< 10 ; _i_++){
      nx_read_length[ _i_ ] = ZeroDefault;
   }
}

Basic_Seq_Statistics::~Basic_Seq_Statistics(){
}

void Basic_Seq_Statistics::resize(){
   read_gc_content_count.resize(PERCENTAGE_ARRAY_SIZE);
   for(int _i_=0; _i_<PERCENTAGE_ARRAY_SIZE; _i_++){
      read_gc_content_count[ _i_ ] = ZeroDefault;
   }
}

// Base class for storing base quality data
Basic_Seq_Statistics::Basic_Seq_Statistics( const Basic_Seq_Statistics& _bss){
   read_gc_content_count = _bss.read_gc_content_count;
   read_length_count.resize(MAX_READ_LENGTH);
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = _bss.read_length_count[ _i_ ];
   }
   nx_read_length.resize(10);
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

   for(size_t _i_=0; _i_<read_gc_content_count.size(); _i_++){
      read_gc_content_count[ _i_ ] = ZeroDefault;
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

void Basic_Seq_Statistics::add(Basic_Seq_Statistics& basic_qc){
    // Update the longest read length
   if ( longest_read_length < basic_qc.longest_read_length){
      longest_read_length = basic_qc.longest_read_length;
   }

    // Update base counts
    total_a_cnt += basic_qc.total_a_cnt;
    total_c_cnt += basic_qc.total_c_cnt;
    total_g_cnt += basic_qc.total_g_cnt;
    total_tu_cnt += basic_qc.total_tu_cnt;
    total_n_cnt += basic_qc.total_n_cnt;

    // Update total number of bases
    this->total_num_bases += basic_qc.total_num_bases;

    // Update read counts
    this->total_num_reads += basic_qc.total_num_reads;

    // Add read lengths
    this->read_lengths.insert(this->read_lengths.end(), basic_qc.read_lengths.begin(), basic_qc.read_lengths.end());

    // Add GC content
    this->read_gc_content_count.insert(this->read_gc_content_count.end(), basic_qc.read_gc_content_count.begin(), basic_qc.read_gc_content_count.end());
}

// Calculates NXX scores and GC content for BAM files
void Basic_Seq_Statistics::global_sum(){
    // Print the size of the read length vector
    std::cout << "GLOBALSUM: Read length vector size: " << this->read_lengths.size() << std::endl;

    if (this->read_lengths.empty()) {
        this->gc_cnt = 0;
        this->longest_read_length = 0;
        this->median_read_length = 0;
        this->mean_read_length = 0;
        this->n50_read_length = 0;
        this->n95_read_length = 0;
        this->n05_read_length = 0;

    } else {
        // Add the G + C bases
        double g_c = this->total_g_cnt + this->total_c_cnt;

        // Add all bases
        double a_tu_g_c = g_c + this->total_a_cnt + this->total_tu_cnt;

        // Check that our total base counts match what was stored (That our code works)
        int _total_num_bases = this->total_num_bases;
        if (a_tu_g_c != (double)_total_num_bases)
        {
            std::cerr << "Total number of bases is not consistent." << std::endl;
            std::cout << _total_num_bases << std::endl;
            std::cout << a_tu_g_c << std::endl;
        } else {
            // Calculate GC-content
            this->gc_cnt = g_c / a_tu_g_c;

            // Sort the read lengths in descending order
            std::vector<int64_t> _read_lengths = this->read_lengths;
            std::sort(_read_lengths.begin(), _read_lengths.end(), std::greater<int>());

            // Get the max read length
            int max_read_length = _read_lengths.at(0);
            this->longest_read_length = max_read_length;

            // Get the median read length
            double _median_read_length = _read_lengths[_read_lengths.size() / 2];
            this->median_read_length = _median_read_length;

            // Get the mean read length
            float _mean_read_length = (double)_total_num_bases / (double)_read_lengths.size();
            this->mean_read_length = _mean_read_length;

            // Calculate N50 and other N-scores
            this->NXX_read_length.resize(101, 0);
            for (int percent_value = 1; percent_value <= 100; percent_value++)
            {
                // Get the base percentage threshold for this N-score
                double base_threshold = (double)_total_num_bases * (percent_value / 100.0);

                // Calculate the NXX score
                double current_base_count = 0;
                int current_read_index = -1;
                while (current_base_count < base_threshold) {
                    current_read_index ++;
                    current_base_count += _read_lengths.at(current_read_index);
                }
                int nxx_read_length = _read_lengths.at(current_read_index);
                this->NXX_read_length[percent_value] = nxx_read_length;
            }

            // Set common score variables
            this->n50_read_length = this->NXX_read_length[50];
            this->n95_read_length = this->NXX_read_length[95];
            this->n05_read_length = this->NXX_read_length[5];
        }
    }
}

// Calculates NXX scores for sequencing_summary.txt files
void Basic_Seq_Statistics::global_sum_no_gc(){
    if (this->read_lengths.empty()) {
        this->longest_read_length = 0;
        this->median_read_length = 0;
        this->mean_read_length = 0;
        this->n50_read_length = 0;
        this->n95_read_length = 0;
        this->n05_read_length = 0;

    } else {
        // Check that our total base counts match what was stored (That our code works)
        int _total_num_bases = this->total_num_bases;

        // Sort the read lengths in descending order
        std::vector<int64_t> _read_lengths = this->read_lengths;
        std::sort(_read_lengths.begin(), _read_lengths.end(), std::greater<int64_t>());

        // Get the max read length
        int64_t max_read_length = _read_lengths.at(0);
        this->longest_read_length = max_read_length;

        // Get the median read length
        double _median_read_length = _read_lengths[_read_lengths.size() / 2];
        this->median_read_length = _median_read_length;

        // Get the mean read length
        float _mean_read_length = (double)_total_num_bases / (double)_read_lengths.size();
        this->mean_read_length = _mean_read_length;

        // Calculate N50 and other N-scores
        this->NXX_read_length.resize(101, 0);
        for (int percent_value = 1; percent_value <= 100; percent_value++)
        {
            // Get the base percentage threshold for this N-score
            double base_threshold = (double)_total_num_bases * (percent_value / 100.0);

            // Calculate the NXX score
            double current_base_count = 0;
            int current_read_index = -1;
            while (current_base_count < base_threshold) {
                current_read_index ++;
                current_base_count += _read_lengths.at(current_read_index);
            }
            int nxx_read_length = _read_lengths.at(current_read_index);
            this->NXX_read_length[percent_value] = nxx_read_length;
        }

        // Set common score variables
        this->n50_read_length = this->NXX_read_length[50];
        this->n95_read_length = this->NXX_read_length[95];
        this->n05_read_length = this->NXX_read_length[5];
    }
}

// Constructor
Basic_Seq_Quality_Statistics::Basic_Seq_Quality_Statistics(){
   pos_quality_distribution.resize(MAX_READ_LENGTH);
   pos_quality_distribution_dev.resize(MAX_READ_LENGTH);
   pos_quality_distribution_count.resize(MAX_READ_LENGTH);
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = ZeroDefault;
      pos_quality_distribution_dev[ _i_ ] = ZeroDefault;
      pos_quality_distribution_count[ _i_ ] = ZeroDefault;
   }
   base_quality_distribution.resize(MAX_BASE_QUALITY);
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] = ZeroDefault;
   }

   read_quality_distribution.resize( MAX_READ_QUALITY );
   for(int _i_=0; _i_<MAX_READ_QUALITY; _i_++){
      read_quality_distribution[ _i_ ] = ZeroDefault;
   }
}

Basic_Seq_Quality_Statistics::~Basic_Seq_Quality_Statistics(){
}

Basic_Seq_Quality_Statistics::Basic_Seq_Quality_Statistics( const Basic_Seq_Quality_Statistics& _bsqs){
    pos_quality_distribution.resize(MAX_READ_LENGTH);
   pos_quality_distribution_dev.resize(MAX_READ_LENGTH);
   pos_quality_distribution_count.resize(MAX_READ_LENGTH);
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      pos_quality_distribution[ _i_ ] = _bsqs.pos_quality_distribution[ _i_ ];
      pos_quality_distribution_dev[ _i_ ] = _bsqs.pos_quality_distribution_dev[ _i_ ];
      pos_quality_distribution_count[ _i_ ] = _bsqs.pos_quality_distribution_count[ _i_ ];
   }
    base_quality_distribution.resize(MAX_BASE_QUALITY);
   for(int _i_=0; _i_<MAX_BASE_QUALITY; _i_++){
      base_quality_distribution[ _i_ ] = _bsqs.base_quality_distribution[ _i_ ];
   }

   min_base_quality = _bsqs.min_base_quality;
   max_base_quality = _bsqs.max_base_quality;

   max_length = _bsqs.max_length;

   read_quality_distribution.resize( MAX_READ_QUALITY );
   for(int _i_=0; _i_<MAX_READ_QUALITY; _i_++){
      read_quality_distribution[ _i_ ] += _bsqs.read_quality_distribution[ _i_ ];
   }
   min_read_quality = _bsqs.min_read_quality;
   max_read_quality = _bsqs.max_read_quality;
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

   for(int _i_=0; _i_<MAX_READ_QUALITY; _i_++){
      read_quality_distribution[ _i_ ] = ZeroDefault;
   }
   min_read_quality = MoneDefault;
   max_read_quality = MoneDefault;
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

   for(int _i_=0; _i_<MAX_READ_QUALITY; _i_++){
      read_quality_distribution[ _i_ ] += t_qual_st.read_quality_distribution[ _i_ ];
   }
   if ( min_read_quality==MoneDefault || min_read_quality > t_qual_st.min_read_quality){
      min_read_quality = t_qual_st.min_read_quality;
   }
   if ( max_read_quality < t_qual_st.max_read_quality){
      max_read_quality = t_qual_st.max_read_quality;
   }
}

void Basic_Seq_Quality_Statistics::global_sum(){
   if ( min_base_quality==MoneDefault){ min_base_quality=ZeroDefault; }
   if ( max_base_quality==MoneDefault){ max_base_quality=ZeroDefault; }
   if ( min_read_quality ==MoneDefault){ min_read_quality=ZeroDefault; }
   if ( max_read_quality ==MoneDefault){ max_read_quality=ZeroDefault; }
}

// BAM output constructor
Output_BAM::Output_BAM(){
   map_quality_distribution.resize( MAX_MAP_QUALITY );
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] = ZeroDefault;
   }
   accuracy_per_read.resize( PERCENTAGE_ARRAY_SIZE );
   for(int _i_=0; _i_<PERCENTAGE_ARRAY_SIZE; _i_++){
      accuracy_per_read[ _i_ ] = ZeroDefault;
   }
}

Output_BAM::~Output_BAM(){
}

void Output_BAM::reset(){
   for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
      map_quality_distribution[ _i_ ] = ZeroDefault;
   }
   for(int _i_=0; _i_<PERCENTAGE_ARRAY_SIZE; _i_++){
      accuracy_per_read[ _i_ ] = ZeroDefault;
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

   long_read_info.reset();
   seq_quality_info.reset();
}

// TODO: Complete setting all metrics
void Output_BAM::add(Output_BAM& output_data){
//    for(int _i_=0; _i_<MAX_MAP_QUALITY; _i_++){
//        map_quality_distribution[ _i_ ] += output_data.map_quality_distribution[ _i_ ];
//    }
//    for(int _i_=0; _i_<PERCENTAGE_ARRAY_SIZE; _i_++){
//        accuracy_per_read[ _i_ ] += output_data.accuracy_per_read[ _i_ ];
//    }

    num_primary_alignment += output_data.num_primary_alignment;
    num_secondary_alignment += output_data.num_secondary_alignment;
//    num_reads_with_secondary_alignment += output_data.num_reads_with_secondary_alignment ;
    num_supplementary_alignment += output_data.num_supplementary_alignment;
    
    // Update the supplementary alignment information
    this->reads_with_supplementary.insert( output_data.reads_with_supplementary.begin(), output_data.reads_with_supplementary.end() );
    this->num_reads_with_supplementary_alignment = this->reads_with_supplementary.size();

    // Update the secondary alignment information
    this->reads_with_secondary.insert( output_data.reads_with_secondary.begin(), output_data.reads_with_secondary.end() );
    this->num_reads_with_secondary_alignment = this->reads_with_secondary.size();

//    num_reads_with_supplementary_alignment += output_data.num_reads_with_supplementary_alignment;
//    num_reads_with_both_secondary_supplementary_alignment += output_data.num_reads_with_both_secondary_supplementary_alignment;

    forward_alignment += output_data.forward_alignment;
    reverse_alignment += output_data.reverse_alignment;
//
//    if ( min_map_quality < 0 || min_map_quality > output_data.min_map_quality){
//      min_map_quality = output_data.min_map_quality;
//    }
//    if ( max_map_quality < output_data.max_map_quality ){
//      max_map_quality =  output_data.max_map_quality;
//    }
//
    num_matched_bases += output_data.num_matched_bases;
    num_mismatched_bases += output_data.num_mismatched_bases;
    num_ins_bases += output_data.num_ins_bases;
    num_del_bases += output_data.num_del_bases;
    num_clip_bases += output_data.num_clip_bases;

    mapped_long_read_info.add(output_data.mapped_long_read_info);
    unmapped_long_read_info.add(output_data.unmapped_long_read_info);
//    mapped_seq_quality_info.add(output_data.mapped_seq_quality_info);
//    unmapped_seq_quality_info.add(output_data.unmapped_seq_quality_info);
//
    long_read_info.add(output_data.mapped_long_read_info);
    long_read_info.add(output_data.unmapped_long_read_info);
//    seq_quality_info.add(output_data.mapped_seq_quality_info);
//    seq_quality_info.add(output_data.unmapped_seq_quality_info);
}

void Output_BAM::global_sum(){
    mapped_long_read_info.global_sum();
    unmapped_long_read_info.global_sum();
    mapped_seq_quality_info.global_sum();
    unmapped_seq_quality_info.global_sum();

    long_read_info.global_sum();
    seq_quality_info.global_sum();

    // Loop through each read and check if it is in both the secondary and supplementary sets
    for ( auto it = this->reads_with_secondary.begin(); it != this->reads_with_secondary.end(); ++it ){
        std::string read_id = it->first;
        if ( this->reads_with_supplementary.find( read_id ) != this->reads_with_supplementary.end() ){
            this->num_reads_with_both_secondary_supplementary_alignment++;
        }
    }

    if ( min_map_quality==MoneDefault){ min_map_quality=ZeroDefault; }
    if ( max_map_quality==MoneDefault){ max_map_quality=ZeroDefault; }
}


// sequencing_summary.txt output
Basic_SeqTxt_Statistics::Basic_SeqTxt_Statistics(){
   signal_range.resize( MAX_SIGNAL_VALUE );
   for(int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] = ZeroDefault;
   }
}

Basic_SeqTxt_Statistics::~Basic_SeqTxt_Statistics(){
}

void Basic_SeqTxt_Statistics::reset(){
   long_read_info.reset();
   seq_quality_info.reset();

   for (int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
       signal_range[ _i_ ] = ZeroDefault;
   }
   min_signal = MoneDefault;
   max_signal = MoneDefault;
}

void Basic_SeqTxt_Statistics::add(Basic_SeqTxt_Statistics& t_output_bSeqTxt){
   long_read_info.add( t_output_bSeqTxt.long_read_info );
   seq_quality_info.add( t_output_bSeqTxt.seq_quality_info );
   for (int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] += t_output_bSeqTxt.signal_range[ _i_ ];
   }
   if ( min_signal==MoneDefault || min_signal > t_output_bSeqTxt.min_signal ) {
      min_signal = t_output_bSeqTxt.min_signal;
   }
   if ( max_signal < t_output_bSeqTxt.max_signal ){
       max_signal = t_output_bSeqTxt.max_signal;
   }
}

void Basic_SeqTxt_Statistics::global_sum(){
   long_read_info.global_sum_no_gc();
   seq_quality_info.global_sum();

   if ( min_signal==MoneDefault){ min_signal = ZeroDefault; }
   if ( max_signal==MoneDefault){ max_signal = ZeroDefault; }
}

void Output_SeqTxt::reset(){
   all_long_read_info.reset();
   passed_long_read_info.reset();
   failed_long_read_info.reset();
}

void Output_SeqTxt::add(Output_SeqTxt& t_output_SeqTxt){
   all_long_read_info.add(t_output_SeqTxt.passed_long_read_info);
   all_long_read_info.add(t_output_SeqTxt.failed_long_read_info);
   passed_long_read_info.add(t_output_SeqTxt.passed_long_read_info);
   failed_long_read_info.add(t_output_SeqTxt.failed_long_read_info);
}

void Output_SeqTxt::global_sum(){
   all_long_read_info.global_sum();
   passed_long_read_info.global_sum();
   failed_long_read_info.global_sum();
}

// Base class for storing a read's base signal data
Base_Signals::Base_Signals(std::string read_name, std::string sequence_data_str, std::vector<std::vector<int>> basecall_signals) {
    this->read_name = read_name;  // Update the read name
    this->sequence_data_str = sequence_data_str;  // Update the sequence string
    this->basecall_signals = basecall_signals;  // Update values
    this->base_count = basecall_signals.size();  // Update read length
}

std::vector<std::vector<int>> Base_Signals::getDataVector() {
    return this->basecall_signals;
}

int Base_Signals::getBaseCount() {
    return this->base_count;
}

std::string Base_Signals::getReadName() {
    return this->read_name;
}

std::string Base_Signals::getSequenceString() {
    return this->sequence_data_str;
}


// FAST5
Output_FAST5::Output_FAST5(){
    this->read_count = 0;
    this->base_count = 0;
}

// Add read base signals
void Output_FAST5::addReadBaseSignals(Base_Signals values){
    this->read_base_signals.push_back(values);  // Update values
    this->read_count++;  // Update read count
    int base_count = values.getBaseCount();
    this->base_count += base_count;  // Update base count
}

// Get the read count
int Output_FAST5::getReadCount(){
    return this->read_count;
}

// Get the total base count across reads
int Output_FAST5::getTotalBaseCount(){
    return this->base_count;
}

// Get the Nth read's name
std::string Output_FAST5::getNthReadName(int read_index){
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::string read_name(signal_data.getReadName());
    return read_name;
}

// Get the Nth read's sequence string
std::string Output_FAST5::getNthReadSequence(int read_index){
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::string sequence_str(signal_data.getSequenceString());
    return sequence_str;
}

// Get the Nth read's base signal data
std::vector<std::vector<int>> Output_FAST5::getNthReadBaseSignals(int read_index){
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::vector<std::vector<int>> data_vector;
    data_vector = signal_data.getDataVector();

    return data_vector;
}

// Get the Nth read's base signal means
std::vector<double> Output_FAST5::getNthReadBaseMeans(int read_index){
    // Get the data vector
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::vector<std::vector<int>> data_vector;
    data_vector = signal_data.getDataVector();

    // Calculate means
    std::vector<double> output;
    output.resize( data_vector.size() );
    std::transform( data_vector.begin(), data_vector.end(), output.begin(), computeMean );

    return output;
}

// Get the Nth read's base signal standard deviations
std::vector<double> Output_FAST5::getNthReadBaseStds(int read_index){
    // Get the data vector
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::vector<std::vector<int>> data_vector;
    data_vector = signal_data.getDataVector();

    // Calculate stds
    std::vector<double> output;
    output.resize( data_vector.size() );
    std::transform( data_vector.begin(), data_vector.end(), output.begin(), computeStd );

    return output;
}

// Get the Nth read's base signal medians
std::vector<double> Output_FAST5::getNthReadBaseMedians(int read_index){
    // Get the data vector
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::vector<std::vector<int>> data_vector;
    data_vector = signal_data.getDataVector();

    // Calculate medians
    std::vector<double> output;
    output.resize( data_vector.size() );
    std::transform( data_vector.begin(), data_vector.end(), output.begin(), computeMedian );

    return output;
}

// Get the Nth read's skewness
std::vector<double> Output_FAST5::getNthReadPearsonSkewnessCoeff(int read_index){
    // Get the data vector
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::vector<std::vector<int>> data_vector;
    data_vector = signal_data.getDataVector();

    // Calculate skewness
    std::vector<double> output;
    output.resize( data_vector.size() );
    std::transform( data_vector.begin(), data_vector.end(), output.begin(), computePearsonsSkewnessCoeff );

    return output;
}

// Get the Nth read's sample kurtosis
std::vector<double> Output_FAST5::getNthReadKurtosis(int read_index){
    // Get the data vector
    Base_Signals signal_data(this->read_base_signals[read_index]);
    std::vector<std::vector<int>> data_vector;
    data_vector = signal_data.getDataVector();

    // Calculate kurtosis
    std::vector<double> output;
    output.resize( data_vector.size() );
    std::transform( data_vector.begin(), data_vector.end(), output.begin(), computeKurtosis );

    return output;
}
