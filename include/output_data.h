#ifndef OUTPUTSTRUCTURES_H
#define OUTPUTSTRUCTURES_H

/*
OutputStructures.h
Define the output structures for each module.
*/

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

#include "input_parameters.h"

#define MAX_READ_LENGTH 10485760
#define MAX_BASE_QUALITY 100
#define MAX_READ_QUALITY 100
#define MAX_SIGNAL_VALUE 5000
#define PERCENTAGE_ARRAY_SIZE 101
#define ZeroDefault 0
#define MoneDefault -1


// Base class for storing error output.
class Output_Info
{
public:
   int error_flag;
   std::string error_str;
   Output_Info();
};


// Base class for storing basic QC data
class Basic_Seq_Statistics
{
public:
    int total_num_reads = ZeroDefault; // total number of long reads
    uint64_t total_num_bases = ZeroDefault; // total number of bases
    int longest_read_length = ZeroDefault; // the length of longest reads
    int n50_read_length = MoneDefault;      // N50
    int n95_read_length = MoneDefault;      // N95
    int n05_read_length = MoneDefault;      // N05;
    double mean_read_length = MoneDefault;      // mean of read length
    std::vector<int> NXX_read_length;       // NXX_read_length[50] means N50 read length
    int median_read_length = MoneDefault;   // median of read length

    uint64_t total_a_cnt = ZeroDefault;  // A content
    uint64_t total_c_cnt = ZeroDefault;  // C content
    uint64_t total_g_cnt = ZeroDefault;  // G content
    uint64_t total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
    uint64_t total_n_cnt = ZeroDefault;  // N content
    double gc_cnt = ZeroDefault;        // GC ratio
    std::vector<int> read_gc_content_count;  // GC content distribution
    std::vector<int> read_length_count;  // Read length distribution (TODO: Replace usages with read_lengths)
    std::vector<int> read_lengths;  // Length of reads

    void add(Basic_Seq_Statistics &t_seq_st);
    void global_sum();
    void global_sum_no_gc();
    void calculate_NXX_scores();  // Calculate N50, N95, N05, etc.
    void resize();
    void addBases(int count);

    Basic_Seq_Statistics();
    Basic_Seq_Statistics(const Basic_Seq_Statistics &_bss);
    ~Basic_Seq_Statistics();
};


// Base class for storing base quality data
class Basic_Seq_Quality_Statistics
{
public:
   //std::vector<uint64_t> base_quality_distribution;
   // Array of base quality distribution initialized to 0
   uint64_t base_quality_distribution[MAX_BASE_QUALITY] = {ZeroDefault};
   std::vector<int> read_average_base_quality_distribution;
   int min_base_quality = MoneDefault; // minimum base quality;
   int max_base_quality = MoneDefault; // maximum base quality;
   std::vector<int> pos_quality_distribution;
   std::vector<double> pos_quality_distribution_dev;
   std::vector<int> pos_quality_distribution_count;
   int max_length = ZeroDefault;

   std::vector<int> read_quality_distribution;
   int min_read_quality = MoneDefault; // the minimum mapping quality
   int max_read_quality = MoneDefault; // the maximum mapping quality

   void add(Basic_Seq_Quality_Statistics &t_qual_st);
   void global_sum();

   Basic_Seq_Quality_Statistics();
   Basic_Seq_Quality_Statistics(const Basic_Seq_Quality_Statistics &_bsqs);
};

// FASTA output
class Output_FA : public Output_Info
{
public:
   Basic_Seq_Statistics long_read_info;
};

// FASTQ output
class Output_FQ : public Output_FA
{
public:
   Basic_Seq_Quality_Statistics seq_quality_info;
};


// Define the base modification data structure (modification type, canonical
// base, likelihood, strand: 0 for forward, 1 for reverse, and CpG flag: T/F)
using Base_Modification = std::tuple<char, char, double, int, bool>;

// BAM output
class Output_BAM : public Output_FQ
{
public:
    uint64_t num_primary_alignment = ZeroDefault;                                 // the number of primary alignment/
    uint64_t num_secondary_alignment = ZeroDefault;                               // the number of secondary alignment
    uint64_t num_reads_with_secondary_alignment = ZeroDefault;                    // the number of long reads with the secondary alignment: one read might have multiple seconard alignment
    uint64_t num_supplementary_alignment = ZeroDefault;                           // the number of supplementary alignment
    uint64_t num_reads_with_supplementary_alignment = ZeroDefault;                // the number of long reads with secondary alignment;
    uint64_t num_reads_with_both_secondary_supplementary_alignment = ZeroDefault; // the number of long reads with both secondary and supplementary alignment.
    uint64_t forward_alignment = ZeroDefault;  // Total number of forward alignments
    uint64_t reverse_alignment = ZeroDefault;  // Total number of reverse alignments
    std::map<std::string, bool> reads_with_supplementary;  // Map of reads with supplementary alignments
    std::map<std::string, bool> reads_with_secondary;  // Map of reads with secondary alignments

    // Similar to Output_FA: below are for mapped.
    uint64_t num_matched_bases = ZeroDefault;    // the number of matched bases with =
    uint64_t num_mismatched_bases = ZeroDefault; // the number of mismatched bases X
    uint64_t num_ins_bases = ZeroDefault;        // the number of inserted bases;
    uint64_t num_del_bases = ZeroDefault;        // the number of deleted bases;
    uint64_t num_clip_bases = ZeroDefault;       // the number of soft-clipped bases;

    // The number of columns can be calculated by summing over the lengths of M/I/D CIGAR operators
    int num_columns = ZeroDefault; // the number of columns
    double percent_identity = ZeroDefault;  // Percent identity = (num columns - NM) / num columns
    std::vector<int> accuracy_per_read;

    // Modified base statistics
    // Number of modified bases by position in the reference:
    // chr -> reference position -> (modification type, canonical base, maximum
    // likelihood, strand)
    std::map<std::string, std::map<int32_t, Base_Modification>> base_modifications;
   //  std::map<int32_t, Base_Modification> base_modifications;
    uint64_t modified_base_count = ZeroDefault;  // Total number of modified bases in the genome
    uint64_t cpg_modified_base_count = ZeroDefault;  // Total number of CpG modified bases in the genome

    Basic_Seq_Statistics mapped_long_read_info;
    Basic_Seq_Statistics unmapped_long_read_info;

    Basic_Seq_Quality_Statistics mapped_seq_quality_info;
    Basic_Seq_Quality_Statistics unmapped_seq_quality_info;

    // Add modified base data
    void add_modification(std::string chr, int32_t ref_pos, char mod_type, char canonical_base, double likelihood, int strand);

    // Return the modification information
    std::map<std::string, std::map<int32_t, Base_Modification>> get_modifications();

    // Add a batch of records to the output
    void add(Output_BAM &t_output_bam);

    // Calculate QC across all records
    void global_sum();

    // Save the output to a summary text file
    void save_summary(std::string &output_file, Input_Para &params, Output_BAM &output_data);

    Output_BAM();
    ~Output_BAM();
};


// sequencing_summary.txt output
class Basic_SeqTxt_Statistics
{
public:
   Basic_Seq_Statistics long_read_info;
   Basic_Seq_Quality_Statistics seq_quality_info;

   std::vector<int> signal_range;
   int min_signal = MoneDefault; // minimum signals;
   int max_signal = MoneDefault; // maximum signals;

   void add(Basic_SeqTxt_Statistics &output_data);
   void global_sum();

   Basic_SeqTxt_Statistics();
};

class Output_SeqTxt : public Output_Info
{
public:
   Basic_SeqTxt_Statistics all_long_read_info;
   Basic_SeqTxt_Statistics passed_long_read_info;
   Basic_SeqTxt_Statistics failed_long_read_info;

   void add(Output_SeqTxt &output_data);
   void global_sum();
};

// Base class for storing a read's base signal data
class Base_Signals
{
public:
    std::string read_name;
    int base_count;
    std::string sequence_data_str;  // Sequence of bases
    std::vector<std::vector<int>> basecall_signals;  // 2D vector of base signals

    // Methods
    int getBaseCount();
    std::string getReadName();
    std::string getSequenceString();
    std::vector<std::vector<int>> getDataVector();
    Base_Signals(std::string read_name, std::string sequence_data_str, std::vector<std::vector<int>> basecall_signals);
};

// FAST5 output
class Output_FAST5 : public Output_FA
{
public:
    // Base quality section
    Basic_Seq_Quality_Statistics seq_quality_info;

    // Signal data section
    int read_count;
    int base_count;
    std::vector<Base_Signals> read_base_signals;

    // Methods
    int getReadCount();
    int getTotalBaseCount();
    std::string getNthReadName(int read_index);
    std::string getNthReadSequence(int read_index);
    void addReadBaseSignals(Base_Signals values);
    void addReadFastq(std::vector<std::string> fq, FILE *read_details_fp);
    std::vector<std::vector<int>> getNthReadBaseSignals(int read_index);
    std::vector<double> getNthReadBaseMeans(int read_index);
    std::vector<double> getNthReadBaseStds(int read_index);
    std::vector<double> getNthReadBaseMedians(int read_index);
    std::vector<double> getNthReadPearsonSkewnessCoeff(int read_index);
    std::vector<double> getNthReadKurtosis(int read_index);

    Output_FAST5();
};

#endif
