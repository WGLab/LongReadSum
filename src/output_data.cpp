#include <numeric>  // std::accumulate
#include<algorithm>  // std::foreach
#include <math.h>  // sqrt
#include <iostream>
#include <sstream>

#include "output_data.h"
#include "utils.h"
#include "basic_statistics.h"

// Base class for storing error output.
Output_Info::Output_Info(){
   error_flag = 0;
   error_str = "";
}

// Base class for storing basic QC data
Basic_Seq_Statistics::Basic_Seq_Statistics(){
    read_length_count.resize(MAX_READ_LENGTH + 1, ZeroDefault);
    read_gc_content_count.resize(101, ZeroDefault);
    NXX_read_length.resize(101, ZeroDefault);
}

Basic_Seq_Statistics::~Basic_Seq_Statistics(){
}

void Basic_Seq_Statistics::resize(){
   read_gc_content_count.resize(PERCENTAGE_ARRAY_SIZE);
   for(int _i_=0; _i_<PERCENTAGE_ARRAY_SIZE; _i_++){
      read_gc_content_count[ _i_ ] = ZeroDefault;
   }
}

void Basic_Seq_Statistics::addBases(int count){
   this->total_num_bases += count;
}

// Base class for storing base quality data
Basic_Seq_Statistics::Basic_Seq_Statistics( const Basic_Seq_Statistics& _bss){
   read_gc_content_count = _bss.read_gc_content_count;
   read_length_count.resize(MAX_READ_LENGTH);
   for(int _i_=0; _i_<MAX_READ_LENGTH; _i_++){
      read_length_count[ _i_ ] = _bss.read_length_count[ _i_ ];
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

void Basic_Seq_Statistics::add(Basic_Seq_Statistics& basic_qc){
    // Update the longest read length
   if ( this->longest_read_length < basic_qc.longest_read_length){
      this->longest_read_length = basic_qc.longest_read_length;
   }

    // Update base counts
    this->total_a_cnt += basic_qc.total_a_cnt;
    this->total_c_cnt += basic_qc.total_c_cnt;
    this->total_g_cnt += basic_qc.total_g_cnt;
    this->total_tu_cnt += basic_qc.total_tu_cnt;
    this->total_n_cnt += basic_qc.total_n_cnt;

    // Update total number of reads
    this->total_num_reads += basic_qc.total_num_reads;

    // Update total number of bases
    this->total_num_bases += basic_qc.total_num_bases;

    // Add read lengths if not empty
    if (!basic_qc.read_lengths.empty()) {
        this->read_lengths.insert(this->read_lengths.end(), basic_qc.read_lengths.begin(), basic_qc.read_lengths.end());
    }

    // Add GC content if not empty
    if (!basic_qc.read_gc_content_count.empty()) {
        this->read_gc_content_count.insert(this->read_gc_content_count.end(), basic_qc.read_gc_content_count.begin(), basic_qc.read_gc_content_count.end());
    }
}

// Calculates NXX scores and GC content for BAM files
void Basic_Seq_Statistics::global_sum(){
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
        uint64_t gc_total = this->total_g_cnt + this->total_c_cnt;

        // Get total base counts
        uint64_t base_total = this->total_a_cnt + this->total_c_cnt + this->total_g_cnt + this->total_tu_cnt + this->total_n_cnt;

        // Calculate read length statistics if base counts are not zero
        if (base_total == 0) {
            std::cerr << "No bases found in input files." << std::endl;
        } else {
            // Calculate GC-content
            double percent_gc = (double) gc_total / base_total;
            this->gc_cnt = percent_gc;

            // Calculate N50 and other N-scores
            this->calculate_NXX_scores();
        }
    }
}

// Calculates NXX scores and other basic read length statistics
void Basic_Seq_Statistics::calculate_NXX_scores(){

    // Sort the read lengths in descending order
    std::sort(this->read_lengths.begin(), this->read_lengths.end(), std::greater<int>());

    // Get total base counts
    uint64_t base_total = this->total_num_bases;

    // Get the max read length
    int64_t max_length = this->read_lengths.at(0);
    this->longest_read_length = max_length;

    // Get the median read length
    double median_length = this->read_lengths[this->read_lengths.size() / 2];
    this->median_read_length = median_length;

    // Get the mean read length
    float mean_length = (double)base_total / (double)this->read_lengths.size();
    this->mean_read_length = mean_length;

    // Initialize the NXX scores
    this->NXX_read_length.resize(100, 0);
    for (int percent_value = 1; percent_value <= 100; percent_value++)
    {
        // Get the base percentage threshold for this N-score
        double base_threshold = (double) (base_total * (double) (percent_value / 100.0));


        // Calculate the NXX score
        double current_base_count = 0;
        int current_read_index = 0;
        while (current_base_count < base_threshold) {
            current_base_count += this->read_lengths.at(current_read_index);
            current_read_index++;
        }
        int nxx_read_index = current_read_index-1;
        int nxx_read_length = this->read_lengths.at(nxx_read_index);
        this->NXX_read_length[percent_value] = nxx_read_length;
    }

    // Set common score variables
    this->n50_read_length = this->NXX_read_length[50];
    this->n95_read_length = this->NXX_read_length[95];
    this->n05_read_length = this->NXX_read_length[5];
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
        // Calculate N50 and other N-scores
        this->calculate_NXX_scores();
    }
}

// Constructor
Basic_Seq_Quality_Statistics::Basic_Seq_Quality_Statistics(){
    pos_quality_distribution.resize(MAX_READ_LENGTH, ZeroDefault);
    pos_quality_distribution_dev.resize(MAX_READ_LENGTH, ZeroDefault);
    pos_quality_distribution_count.resize(MAX_READ_LENGTH, ZeroDefault);

    read_average_base_quality_distribution.resize(MAX_READ_QUALITY, ZeroDefault);
    read_quality_distribution.resize(MAX_READ_QUALITY, ZeroDefault);
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
}

Output_BAM::~Output_BAM(){
}

void Output_BAM::add_modification(int32_t ref_pos, char mod_type, char canonical_base, double likelihood, bool is_cpg)
{
    try {
        this->base_modifications.at(ref_pos);

        // If the reference position is already in the map, use the modification
        // type with the highest likelihood
        double previous_likelihood = std::get<2>(this->base_modifications[ref_pos]);
        if (likelihood > previous_likelihood){
            this->base_modifications[ref_pos] = std::make_tuple(mod_type, canonical_base, likelihood);
        }
    } catch (const std::out_of_range& oor) {

        // If the reference position is not in the map, add the modification
        this->base_modifications[ref_pos] = std::make_tuple(mod_type, canonical_base, likelihood);
    }
}

std::map<int32_t, std::tuple<char, char, double>> Output_BAM::get_modifications()
{
    return this->base_modifications;
}

void Output_BAM::add(Output_BAM &output_data)
{
    this->num_primary_alignment += output_data.num_primary_alignment;
    this->num_secondary_alignment += output_data.num_secondary_alignment;
    this->num_supplementary_alignment += output_data.num_supplementary_alignment;

    // Update the secondary alignment information
    this->reads_with_secondary.insert( output_data.reads_with_secondary.begin(), output_data.reads_with_secondary.end() );
    this->num_reads_with_secondary_alignment = this->reads_with_secondary.size();

    // Update the supplementary alignment information
    this->reads_with_supplementary.insert( output_data.reads_with_supplementary.begin(), output_data.reads_with_supplementary.end() );
    this->num_reads_with_supplementary_alignment = this->reads_with_supplementary.size();

    // Update the forward and reverse alignment information
    this->forward_alignment += output_data.forward_alignment;
    this->reverse_alignment += output_data.reverse_alignment;

//    // Resize the base quality vector if it is empty
//    if ( this->seq_quality_info.base_quality_distribution.empty() ){
//        this->seq_quality_info.base_quality_distribution.resize( MAX_READ_QUALITY );
//    }

    // Update the base quality vector if it is not empty
//    if ( !output_data.seq_quality_info.base_quality_distribution.empty() ){
    for (int i=0; i<MAX_READ_QUALITY; i++){
        this->seq_quality_info.base_quality_distribution[i] += output_data.seq_quality_info.base_quality_distribution[i];
    }

    this->num_matched_bases += output_data.num_matched_bases;
    this->num_mismatched_bases += output_data.num_mismatched_bases;
    this->num_ins_bases += output_data.num_ins_bases;
    this->num_del_bases += output_data.num_del_bases;
    this->num_clip_bases += output_data.num_clip_bases;

    this->mapped_long_read_info.add(output_data.mapped_long_read_info);
    this->unmapped_long_read_info.add(output_data.unmapped_long_read_info);

    this->long_read_info.add(output_data.mapped_long_read_info);
    this->long_read_info.add(output_data.unmapped_long_read_info);

    // Update the base modification information
    printMessage("Adding base modification information to the output data, size: " + std::to_string(output_data.base_modifications.size()));
    for (auto const &it : output_data.base_modifications) {
        int32_t ref_pos = it.first;
        char mod_type = std::get<0>(it.second);
        char canonical_base = std::get<1>(it.second);
        double likelihood = std::get<2>(it.second);
        this->add_modification(ref_pos, mod_type, canonical_base, likelihood, false);
    }
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
}

// Save the output to a file
void Output_BAM::save_summary(std::string &output_file, Input_Para &params, Output_BAM &output_data){
    FILE *fp = fopen(output_file.c_str(), "w");
    if (fp == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", output_file.c_str());
    } else {
        // Save basic statistics
        fprintf(fp, "Total number of reads\t%d\n", output_data.long_read_info.total_num_reads);
        fprintf(fp, "Total number of bases\t%ld\n", output_data.long_read_info.total_num_bases);
        fprintf(fp, "Longest read length\t%d\n", output_data.long_read_info.longest_read_length);
        fprintf(fp, "N50 read length\t%d\n", output_data.long_read_info.n50_read_length);
        fprintf(fp, "Mean read length\t%.2f\n", output_data.long_read_info.mean_read_length);
        fprintf(fp, "Median read length\t%d\n", output_data.long_read_info.median_read_length);
        fprintf(fp, "GC%%\t%.2f\n", output_data.long_read_info.gc_cnt * 100);
        fprintf(fp, "\n");

        // Save the mapping statistics
        fprintf(fp, "Total number of mapped reads\t%d\n", output_data.mapped_long_read_info.total_num_reads);
        fprintf(fp, "Total number of mapped bases\t%ld\n", output_data.mapped_long_read_info.total_num_bases);
        fprintf(fp, "Longest mapped read length\t%d\n", output_data.mapped_long_read_info.longest_read_length);
        fprintf(fp, "N50 mapped read length\t%d\n", output_data.mapped_long_read_info.n50_read_length);
        fprintf(fp, "Mean mapped read length\t%.2f\n", output_data.mapped_long_read_info.mean_read_length);
        fprintf(fp, "Median mapped read length\t%d\n", output_data.mapped_long_read_info.median_read_length);
        fprintf(fp, "GC%%\t%.2f\n", output_data.mapped_long_read_info.gc_cnt * 100);
        fprintf(fp, "\n");

        // Save the read alignment statistics
        fprintf(fp, "Total number of primary alignments\t%ld\n", output_data.num_primary_alignment);
        fprintf(fp, "Total number of secondary alignments\t%ld\n", output_data.num_secondary_alignment);
        fprintf(fp, "Total number of supplementary alignments\t%ld\n", output_data.num_supplementary_alignment);
        fprintf(fp, "Total number of reads with secondary alignments\t%ld\n", output_data.num_reads_with_secondary_alignment);
        fprintf(fp, "Total number of reads with supplementary alignments\t%ld\n", output_data.num_reads_with_supplementary_alignment);
        fprintf(fp, "Total number of reads with both secondary and supplementary alignments\t%ld\n", output_data.num_reads_with_both_secondary_supplementary_alignment);
        fprintf(fp, "Total number of reads with forward alignments\t%ld\n", output_data.forward_alignment);
        fprintf(fp, "Total number of reads with reverse alignments\t%ld\n", output_data.reverse_alignment);
        fprintf(fp, "Total number of reverse alignment\t%ld\n", output_data.reverse_alignment);
        fprintf(fp, "\n");

        // Save the base alignment statistics
        fprintf(fp, "Total number of matched bases\t%ld\n", output_data.num_matched_bases);
        fprintf(fp, "Total number of mismatched bases\t%ld\n", output_data.num_mismatched_bases);
        fprintf(fp, "Total number of insertions\t%ld\n", output_data.num_ins_bases);
        fprintf(fp, "Total number of deletions\t%ld\n", output_data.num_del_bases);
        fprintf(fp, "Total number of primary alignment clipped bases (soft + hard)\t%ld\n", output_data.num_clip_bases);

        // Close the file
        fclose(fp);
    }
}

// sequencing_summary.txt output
Basic_SeqTxt_Statistics::Basic_SeqTxt_Statistics(){
   signal_range.resize( MAX_SIGNAL_VALUE );
   for(int _i_=0; _i_<MAX_SIGNAL_VALUE; _i_++){
      signal_range[ _i_ ] = ZeroDefault;
   }
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

// Add read fastq data
void Output_FAST5::addReadFastq(std::vector<std::string> fq, FILE *read_details_fp)
{
    const char * read_name;
    double gc_content_pct;

    // Access the read name
    std::string header_str = fq[0];
    std::istringstream iss_header( header_str );
    std::string read_name_str;
    std::getline( iss_header, read_name_str, ' ' );
    read_name = read_name_str.c_str();

    // Access the sequence data
    std::string sequence_data_str = fq[1];

    // Update the total number of bases
    int base_count = sequence_data_str.length();
    long_read_info.total_num_bases += base_count;

    // Store the read length
    long_read_info.read_lengths.push_back(base_count);

    // Access base quality data
    char value;
    std::vector<int> base_quality_values;
    std::string base_quality_str = fq[3];
    std::istringstream iss( base_quality_str );
    while (iss >> value) {
        int base_quality_value = value - '!';  // '!' symbol represent 0-quality score
        base_quality_values.push_back(base_quality_value);
    }

    // Update the base quality and GC content information
    int gc_count = 0;
    double read_mean_base_qual = 0;
    char current_base;
    uint64_t base_quality_value;
    for (int i = 0; i < base_count; i++)
    {
        current_base = sequence_data_str[i];
        if (current_base == 'A' || current_base == 'a')
        {
            long_read_info.total_a_cnt += 1;
        }
        else if (current_base == 'G' || current_base == 'g')
        {
            long_read_info.total_g_cnt += 1;
            gc_count += 1;
        }
        else if (current_base == 'C' || current_base == 'c')
        {
            long_read_info.total_c_cnt += 1;
            gc_count += 1;
        }
        else if (current_base == 'T' || current_base == 't' || current_base == 'U' || current_base == 'u')
        {
            long_read_info.total_tu_cnt += 1;
        }
        // Get the base quality
        base_quality_value = (uint64_t)base_quality_values[i];
        seq_quality_info.base_quality_distribution[base_quality_value] += 1;
        read_mean_base_qual += (double)base_quality_value;
    }

    // Calculate percent guanine & cytosine
    gc_content_pct = 100.0 *( (double)gc_count / (double)base_count );

    // Look into this section
    long_read_info.read_gc_content_count[(int)(gc_content_pct + 0.5)] += 1;
    read_mean_base_qual /= (double) base_count;
    seq_quality_info.read_average_base_quality_distribution[(uint)(read_mean_base_qual + 0.5)] += 1;
    fprintf(read_details_fp, "%s\t%d\t%.2f\t%.2f\n", read_name, base_count, gc_content_pct, read_mean_base_qual);

    // Update the total number of reads
    long_read_info.total_num_reads += 1;
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
