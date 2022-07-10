#ifndef CXX_TO_PYTHON_CLASS_H
#define CXX_TO_PYTHON_CLASS_H

/*
CXX_to_python_class.h:
Define the C++ bindings from our Python modules
*/

#include <string>
#include <vector>

//#include <stdint.h>
//#include <cstddef>

//                        ,  ,
#define MAX_READ_LENGTH 10485760
#define MAX_MAP_QUALITY 256
#define MAX_BASE_QUALITY 256
#define MAX_READ_QUALITY 256
#define MAX_SIGNAL_VALUE 5000

#define PERCENTAGE_ARRAY_SIZE 101

#define ZeroDefault 0
#define MoneDefault -1


/*
Output_FQ: Base class which contains error output.
*/
class Output_Info
{
public:
   int error_flag;
   std::string error_str;
   Output_Info();
};



/*
Basic_Seq_Statistics: Python-invoked C++ class which contains basic statistics from the C++ module
*/
class Basic_Seq_Statistics
{
public:
   int64_t total_num_reads = ZeroDefault; // total number of long reads
   int64_t total_num_bases = ZeroDefault; // total number of bases

   uint64_t longest_read_length = ZeroDefault; // the length of longest reads
   int64_t n50_read_length = MoneDefault;      // N50
   int64_t n95_read_length = MoneDefault;      // N95
   int64_t n05_read_length = MoneDefault;      // N05;
   double mean_read_length = MoneDefault;      // mean of read length
   std::vector<int64_t> nx_read_length;        // deprecated
   std::vector<int64_t> NXX_read_length;       // NXX_read_length[50] means N50 read length
   int64_t median_read_length = MoneDefault;   // median of read length

   int64_t total_a_cnt = ZeroDefault;  // A content
   int64_t total_c_cnt = ZeroDefault;  // C content
   int64_t total_g_cnt = ZeroDefault;  // G content
   int64_t total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
   int64_t total_n_cnt = ZeroDefault;  // N content
   double gc_cnt = ZeroDefault;        // GC ratio

   std::vector<int64_t> read_length_count;
   // read_length_count[x] is the number of reads that length is equal to x. read_length_count[MAX_READ_LENGTH] includes the number of reads that are longer than MAX_READ_LENGTH. size of read_length_count is MAX_READ_LENGTH+1

   std::vector<int64_t> read_gc_content_count;

   void reset();
   void add(Basic_Seq_Statistics &t_seq_st);
   void global_sum();
   void resize();

   Basic_Seq_Statistics();
   Basic_Seq_Statistics(const Basic_Seq_Statistics &_bss);
   ~Basic_Seq_Statistics();
};

class Output_FA : public Output_Info
{
public:
   Basic_Seq_Statistics long_read_info;
};



class Basic_Seq_Quality_Statistics
{
public:
   std::vector<int64_t> base_quality_distribution;
   std::vector<int64_t> read_average_base_quality_distribution;
   int min_base_quality = MoneDefault; // minimum base quality;
   int max_base_quality = MoneDefault; // maximum base quality;
   std::vector<int> pos_quality_distribution;
   std::vector<double> pos_quality_distribution_dev;
   std::vector<int64_t> pos_quality_distribution_count;
   int64_t max_length = ZeroDefault;

   std::vector<int64_t> read_quality_distribution;
   int min_read_quality = MoneDefault; // the minimum mapping quality
   int max_read_quality = MoneDefault; // the maximum mapping quality

   void reset();
   void add(Basic_Seq_Quality_Statistics &t_qual_st);
   void global_sum();

   Basic_Seq_Quality_Statistics();
   Basic_Seq_Quality_Statistics(const Basic_Seq_Quality_Statistics &_bsqs);
   ~Basic_Seq_Quality_Statistics();
};


class Output_FQ : public Output_FA
{
public:
   Basic_Seq_Quality_Statistics seq_quality_info;
};


class Output_BAM : public Output_FQ
{
public:
   int64_t num_primary_alignment = ZeroDefault;                                 // the number of primary alignment/
   int64_t num_secondary_alignment = ZeroDefault;                               // the number of secondary alignment
   int64_t num_reads_with_secondary_alignment = ZeroDefault;                    // the number of long reads with the secondary alignment: one read might have multiple seconard alignment
   int64_t num_supplementary_alignment = ZeroDefault;                           // the number of supplementary alignment
   int64_t num_reads_with_supplementary_alignment = ZeroDefault;                // the number of long reads with secondary alignment;
   int64_t num_reads_with_both_secondary_supplementary_alignment = ZeroDefault; // the number of long reads with both secondary and supplementary alignment.

   int64_t forward_alignment = ZeroDefault;
   int64_t reverse_alignment = ZeroDefault;

   //int64_t map_quality_distribution[MAX_MAP_QUALITY]; // The mapping quality distribution;
   std::vector<int64_t> map_quality_distribution;
   int min_map_quality = MoneDefault; // the minimum mapping quality
   int max_map_quality = MoneDefault; // the maximum mapping quality

   // Similar to Output_FA: below are for mapped.
   //int64_t num_mapped_reads = ZeroDefault; // the number of mapped long reads
   //int64_t num_mapped_bases = ZeroDefault; // the number of mapped bases with X/=/M
   int64_t num_matched_bases = ZeroDefault;    // the number of matched bases with =
   int64_t num_mismatched_bases = ZeroDefault; // the number of mismatched bases X
   int64_t num_ins_bases = ZeroDefault;        // the number of inserted bases;
   int64_t num_del_bases = ZeroDefault;        // the number of deleted bases;
   int64_t num_clip_bases = ZeroDefault;       // the number of soft-clipped bases;

   //int64_t accuracy_per_read[PERCENTAGE_ARRAY_SIZE]; //
   std::vector<int64_t> accuracy_per_read;

   Basic_Seq_Statistics mapped_long_read_info;
   Basic_Seq_Statistics unmapped_long_read_info;

   Basic_Seq_Quality_Statistics mapped_seq_quality_info;
   Basic_Seq_Quality_Statistics unmapped_seq_quality_info;

   void reset();
   void add(Output_BAM &t_output_bam);
   void global_sum();

   Output_BAM();
   ~Output_BAM();
};

class Basic_SeqTxt_Statistics
{
public:
   Basic_Seq_Statistics long_read_info;
   Basic_Seq_Quality_Statistics seq_quality_info;

   //int signal_range[MAX_SIGNAL_VALUE]; // statistics of all signals
   std::vector<int> signal_range;
   int min_signal = MoneDefault; // minimum signals;
   int max_signal = MoneDefault; // maximum signals;

   //int64_t read_quality_distribution[MAX_READ_QUALITY]; // The mapping quality distribution;
   /*std::vector<int64_t> read_quality_distribution;
   int min_read_quality = MoneDefault; // the minimum mapping quality
   int max_read_quality = MoneDefault; // the maximum mapping quality
   */

   // int64_t num_reads = ZeroDefault; // The number of long reads
   // int64_t *read_length_list; // statistics of read length for long reads: each position is the number of reads with the length of the index;

   void reset();
   void add(Basic_SeqTxt_Statistics &t_output_bSeqTxt);
   void global_sum();

   Basic_SeqTxt_Statistics();
   ~Basic_SeqTxt_Statistics();
};

class Output_SeqTxt : public Output_Info
{
public:
   Basic_SeqTxt_Statistics all_long_read_info;
   Basic_SeqTxt_Statistics passed_long_read_info;
   Basic_SeqTxt_Statistics failed_long_read_info;

   void reset();
   void add(Output_SeqTxt &t_output_SeqTxt);
   void global_sum();
};

#endif
