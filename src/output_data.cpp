#include <numeric>  // std::accumulate
#include <algorithm>  // std::foreach
#include <math.h>  // sqrt
#include <iostream>
#include <sstream>
#include <cmath>  // std::round

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

    // Update the per-read GC content distribution
    for (int i = 0; i < 101; i++) {
        this->read_gc_content_count[i] += basic_qc.read_gc_content_count[i];
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
    this->num_primary_alignment = 0;
    this->num_secondary_alignment = 0;
    this->num_supplementary_alignment = 0;
    this->num_clip_bases = 0;
    this->sample_modified_base_count = 0;
    this->sample_modified_base_count_forward = 0;
    this->sample_modified_base_count_reverse = 0;
    this->forward_alignment = 0;
    this->reverse_alignment = 0;
    this->base_mod_counts = std::unordered_map<char, uint64_t>();
    this->base_mod_counts_forward = std::unordered_map<char, uint64_t>();
    this->base_mod_counts_reverse = std::unordered_map<char, uint64_t>();
}

Output_BAM::~Output_BAM(){
}

void Output_BAM::updateBaseModCounts(char mod_type, int strand)
{
    // Update the sample modified base count for predictions exceeding the threshold
    this->sample_modified_base_count++;
    this->base_mod_counts[mod_type]++;  // Update the type-specific modified base count

    // Update the modified base count for the strand from primary alignments
    if (strand == 0) {
        this->sample_modified_base_count_forward++;
        this->base_mod_counts_forward[mod_type]++;  // Update the type-specific modified base count
    } else if (strand == 1) {
        this->sample_modified_base_count_reverse++;
        this->base_mod_counts_reverse[mod_type]++;  // Update the type-specific modified base count
    }
}

void Output_BAM::updateBaseModProbabilities(char mod_type, double pct_len, double probability)
{
    // Update the base modification probabilities
    this->read_pct_len_vs_mod_prob[mod_type].push_back(std::make_pair(pct_len, probability));
}

void Output_BAM::updateReadModRate(int read_length, const std::unordered_map<char, double>& base_mod_rates) {
    ReadModData read_mod_data;
    read_mod_data.read_length = read_length;
    read_mod_data.base_mod_rates = base_mod_rates;
    this->read_mod_data.push_back(read_mod_data);
}

std::vector<char> Output_BAM::getBaseModTypes()
{
    std::vector<char> base_mod_types;
    if (this->base_mod_counts.empty()) {
        printError("No base modification counts found.");
        return base_mod_types;
    }

    for (const auto& it : this->base_mod_counts) {
        base_mod_types.push_back(it.first);
    }
    
    return base_mod_types;
}

int Output_BAM::getReadModDataSize()
{
    return this->read_mod_data.size();
}

int Output_BAM::getNthReadModLength(int read_index)
{
    return this->read_mod_data[read_index].read_length;
}

double Output_BAM::getNthReadModRate(int read_index, char mod_type)
{
    double mod_rate = 0.0;
    try {
        this->read_mod_data.at(read_index);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Read index " << read_index << " is out of range." << std::endl;
    }
    try {
        mod_rate = this->read_mod_data[read_index].base_mod_rates.at(mod_type);
    } catch (const std::out_of_range& oor) {
        // No modification rate found for the specified type in the read
        mod_rate = 0.0;
    }
    return mod_rate;
}

uint64_t Output_BAM::getModTypeCount(char mod_type)
{
    return this->base_mod_counts[mod_type];
}

uint64_t Output_BAM::getModTypeCount(char mod_type, int strand)
{
    if (strand == 0) {
        return this->base_mod_counts_forward[mod_type];
    } else {
        return this->base_mod_counts_reverse[mod_type];
    }
}

double Output_BAM::getNthReadLenPct(int read_index, char mod_type)
{
    double read_len_pct = 0.0;
    try {
        this->read_pct_len_vs_mod_prob.at(mod_type);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Read length percentage not found for type " << mod_type << std::endl;
    }
    try {
        read_len_pct = this->read_pct_len_vs_mod_prob[mod_type].at(read_index).first;
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Read length percentage not found for read index " << read_index << " and type " << mod_type << std::endl;
        return 0.0;
    }
    return read_len_pct;
}

double Output_BAM::getNthReadModProb(int read_index, char mod_type)
{
    double mod_prob = -1.0;
    try {
        this->read_pct_len_vs_mod_prob.at(mod_type);
    } catch (const std::out_of_range& oor) {
        return mod_prob;
    }
    try {
        mod_prob = this->read_pct_len_vs_mod_prob[mod_type].at(read_index).second;
    } catch (const std::out_of_range& oor) {
        // std::cerr << "Error: Modification probability not found for read index " << read_index << " and type " << mod_type << std::endl;
        return -1.0;
    }
    return mod_prob;
}

int Output_BAM::getReadCount()
{
    return this->read_move_table.size();
}

void Output_BAM::addReadMoveTable(std::string read_name, std::string sequence_data_str, std::vector<int> signal_index, int start, int end)
{
    Base_Move_Table values(sequence_data_str, signal_index, start, end);
    this->read_move_table[read_name] = values;
}

std::vector<int> Output_BAM::getReadMoveTable(std::string read_id)
{
    try {
        this->read_move_table.at(read_id);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Read name " << read_id << " is not in the move table." << std::endl;
    }
    return this->read_move_table[read_id].getBaseSignalIndex();
}

// Get the read's sequence string
std::string Output_BAM::getReadSequence(std::string read_id)
{
    try {
        this->read_move_table.at(read_id);
    } catch (const std::out_of_range& oor) {
        std::cerr << "Error: Read name " << read_id << " is not in the move table." << std::endl;
    }

    Base_Move_Table signal_data = this->read_move_table[read_id];
    std::string sequence_str(signal_data.getSequenceString());
    return sequence_str;
}

int Output_BAM::getReadSequenceStart(std::string read_id)
{
    return this->read_move_table[read_id].getSequenceStart();
}

int Output_BAM::getReadSequenceEnd(std::string read_id)
{
    return this->read_move_table[read_id].getSequenceEnd();
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

    // Update the base quality vector if it is not empty
    for (int i=0; i<MAX_READ_QUALITY; i++){
        this->seq_quality_info.base_quality_distribution[i] += output_data.seq_quality_info.base_quality_distribution[i];
    }

    // Update the read average base quality vector if it is not empty
    for (int i=0; i<MAX_READ_QUALITY; i++){
        this->seq_quality_info.read_average_base_quality_distribution[i] += output_data.seq_quality_info.read_average_base_quality_distribution[i];
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

    // Update base modification counts
    this->modified_prediction_count += output_data.modified_prediction_count;

    // Update the map of read IDs to base signal data
    for ( auto it = output_data.read_move_table.begin(); it != output_data.read_move_table.end(); ++it ){
        std::string read_id = it->first;
        std::vector<int> signal_index = it->second.getBaseSignalIndex();
        std::string sequence_data_str = it->second.getSequenceString();
        int start = it->second.getSequenceStart();
        int end = it->second.getSequenceEnd();
        this->addReadMoveTable(read_id, sequence_data_str, signal_index, start, end);
    }

    // Preprint revisions: Update base modification counts
    this->sample_modified_base_count += output_data.sample_modified_base_count;
    this->sample_modified_base_count_forward += output_data.sample_modified_base_count_forward;
    this->sample_modified_base_count_reverse += output_data.sample_modified_base_count_reverse;
    for ( auto it = output_data.sample_c_modified_positions.begin(); it != output_data.sample_c_modified_positions.end(); ++it ){
        std::string chr = it->first;
        std::vector<std::pair<int32_t, int>> positions = it->second;
        if (this->sample_c_modified_positions.find(chr) == this->sample_c_modified_positions.end()){
            this->sample_c_modified_positions[chr] = positions;
        } else {
            this->sample_c_modified_positions[chr].insert(this->sample_c_modified_positions[chr].end(), positions.begin(), positions.end());
        }
    }
}

void Output_BAM::addTINData(std::string &bam_file, TINStats &tin_data) {
    this->tin_data[bam_file] = tin_data;
}

double Output_BAM::getTINMean(std::string bam_file) {
    return this->tin_data[bam_file].mean;
}

double Output_BAM::getTINMedian(std::string bam_file) {
    return this->tin_data[bam_file].median;
}

double Output_BAM::getTINStdDev(std::string bam_file) {
    return this->tin_data[bam_file].stddev;
}

int Output_BAM::getTINCount(std::string bam_file) {
    return this->tin_data[bam_file].num_transcripts;
}

void Output_BAM::global_sum(){
    // Calculate the global sums for the basic statistics
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

void Output_SeqTxt::save_summary(std::string & output_file, Input_Para & params)
{
    
    FILE *fp = fopen(output_file.c_str(), "w");
    if (fp == NULL){
        fprintf(stderr, "Error: cannot open file %s\n", output_file.c_str());
    } else {
        // Define the types explicitly
        using ReadInfo = std::tuple<const char*, Basic_Seq_Statistics&>;

        // Save basic statistics for total, passed, and failed reads
        for (const ReadInfo& read_type : {
            ReadInfo("All", all_long_read_info.long_read_info),
            ReadInfo("Passed", passed_long_read_info.long_read_info),
            ReadInfo("Failed", failed_long_read_info.long_read_info)
        }) {
            std::string read_filter = std::get<0>(read_type);
            Basic_Seq_Statistics& long_read_info = std::get<1>(read_type);

            fprintf(fp, "%s reads:\n", read_filter.c_str());
            fprintf(fp, "Total number of reads\t%d\n", long_read_info.total_num_reads);
            fprintf(fp, "Total number of bases\t%ld\n", long_read_info.total_num_bases);
            fprintf(fp, "Longest read length\t%d\n", long_read_info.longest_read_length);
            fprintf(fp, "N50 read length\t%d\n", long_read_info.n50_read_length);
            fprintf(fp, "Mean read length\t%.2f\n", long_read_info.mean_read_length);
            fprintf(fp, "Median read length\t%d\n", long_read_info.median_read_length);
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

// Base class for storing a read's base signal data
Base_Signals::Base_Signals(std::string read_name, std::string sequence_data_str, std::vector<std::vector<int>> basecall_signals) {
    this->read_name = read_name;
    this->sequence_data_str = sequence_data_str;
    this->basecall_signals = basecall_signals;
    this->base_count = basecall_signals.size();
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

    // Access the read name
    std::string header_str = fq[0];
    std::istringstream iss_header( header_str );
    std::string read_name_str;
    std::getline( iss_header, read_name_str, ' ' );
    read_name = read_name_str.c_str();
    std::string sequence_data_str = fq[1];  // Access the sequence data

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

    // Ensure the base quality values match the sequence length
    if (base_quality_values.size() != base_count) {
        printError("Warning: Base quality values do not match the sequence length for read ID " + std::string(read_name));
    }

    // Update the base quality and GC content information
    int gc_count = 0;
    double cumulative_base_prob = 0;  // Read cumulative base quality probability
    char current_base;
    int base_quality_value;
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
        // Get the base quality (Phred) value
        base_quality_value = base_quality_values[i];

        // Update the per-base quality distribution
        try {
            seq_quality_info.base_quality_distribution[base_quality_value] += 1;
        } catch (const std::out_of_range& oor) {
            printError("Warning: Base quality value " + std::to_string(base_quality_value) + " exceeds maximum value");
        }

        // Convert the Phred quality value to a probability
        double base_quality_prob = pow(10, -base_quality_value / 10.0);
        cumulative_base_prob += base_quality_prob;
    }

    // Calculate the mean base quality probability
    cumulative_base_prob /= (double)base_count;

    // Convert the mean base quality probability to a Phred quality value
    double read_mean_base_qual = -10.0 * log10(cumulative_base_prob);

    // Update the per-read GC content distribution
    double gc_content_pct = (100.0 * gc_count) / static_cast<double>(base_count);
    int gc_content_int = static_cast<int>(std::round(gc_content_pct));
    try {
        long_read_info.read_gc_content_count[gc_content_int] += 1;
    } catch (const std::out_of_range& oor) {
        printError("Warning: Invalid GC content value " + std::to_string(gc_content_int));
    }

    // Update the per-read base quality distribution
    int read_mean_base_qual_int = static_cast<int>(std::round(read_mean_base_qual));

    try {
        seq_quality_info.read_quality_distribution[read_mean_base_qual_int] += 1;
    } catch (const std::out_of_range& oor) {
        printError("Warning: Base quality value " + std::to_string(read_mean_base_qual_int) + " exceeds maximum value");
    }

    fprintf(read_details_fp, "%s\t%d\t%.2f\t%.2f\n", read_name, base_count, gc_content_pct, read_mean_base_qual);  // Write to file

    long_read_info.total_num_reads += 1;  // Update read count
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

std::string Base_Move_Table::getSequenceString()
{
    return this->sequence_data_str;
}

std::vector<int> Base_Move_Table::getBaseSignalIndex()
{
    return this->base_signal_index;
}

int Base_Move_Table::getSequenceStart()
{
    return this->sequence_start;
}

int Base_Move_Table::getSequenceEnd()
{
    return this->sequence_end;
}

Base_Move_Table::Base_Move_Table(std::string sequence_data_str, std::vector<int> base_signal_index, int start, int end)
{
    this->sequence_data_str = sequence_data_str;
    this->base_signal_index = base_signal_index;
    this->sequence_start = start;
    this->sequence_end = end;
}

Base_Move_Table::Base_Move_Table()
{
}
