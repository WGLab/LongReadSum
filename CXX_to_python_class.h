#ifndef CXX_TO_PYTHON_CLASS_H
#define CXX_TO_PYTHON_CLASS_H
// CXX_to_python_class

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

class Output_Info{
public:
    int error_flag;
    std::string error_str;
    Output_Info();
};

class Basic_Seq_Statistics{
public:
   int64_t total_num_reads = ZeroDefault; // total number of long reads
   int64_t total_num_bases = ZeroDefault; // total number of bases

   uint64_t longest_read_length = ZeroDefault; // the length of longest reads
   int64_t n50_read_length = MoneDefault; // N50
   int64_t n95_read_length = MoneDefault; // N95
   int64_t n05_read_length = MoneDefault; // N05;
   int64_t mean_read_length = MoneDefault; // mean of read length
   //int64_t nx_read_length[10]; // can provide n10, n20, n30, ... n90; The last indicates whether this info is generated.
   std::vector<int64_t> nx_read_length;
   int64_t median_read_length = MoneDefault; // median of read length

   int64_t total_a_cnt = ZeroDefault; // A content
   int64_t total_c_cnt = ZeroDefault; // C content
   int64_t total_g_cnt = ZeroDefault; // G content
   int64_t total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
   int64_t total_n_cnt = ZeroDefault; // N content
   double gc_cnt = ZeroDefault; // GC ratio

   //int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;
   std::vector<int64_t> read_length_count; 
   //int64_t read_length_count[MAX_READ_LENGTH];

   std::vector<int64_t> read_quality_distribution;
   int min_read_quality = MoneDefault; // the minimum mapping quality
   int max_read_quality = MoneDefault; // the maximum mapping quality

   void reset();
   void add(Basic_Seq_Statistics& t_seq_st);
   void global_sum();

   Basic_Seq_Statistics();
   Basic_Seq_Statistics( const Basic_Seq_Statistics& _bss);
   ~Basic_Seq_Statistics();
};

class Output_FA : public Output_Info {
public:
   Basic_Seq_Statistics long_read_info;
};

class Basic_Seq_Quality_Statistics{
public:
   //int64_t base_quality_distribution[MAX_BASE_QUALITY];  // base quality distribution: each position is the number of bases with the quality of the index.
   std::vector<int64_t> base_quality_distribution;
   int min_base_quality = MoneDefault; // minimum base quality;
   int max_base_quality = MoneDefault; // maximum base quality;

   //int *pos_quality_distribution; // the base quality according to positions of long reads: each position is the mean quality at the position of long reads
   //double *pos_quality_distribution_dev; // similar to above but for the standard deviation
   //int64_t *pos_quality_distribution_count; // similar to above but for the number of bases
   std::vector<int> pos_quality_distribution; 
   std::vector<double> pos_quality_distribution_dev;
   std::vector<int64_t> pos_quality_distribution_count;
   //int pos_quality_distribution[MAX_READ_LENGTH]; 
   //double pos_quality_distribution_dev[MAX_READ_LENGTH]; 
   //int64_t pos_quality_distribution_count[MAX_READ_LENGTH];
   int64_t max_length = ZeroDefault;

   void reset();
   void add(Basic_Seq_Quality_Statistics& t_qual_st);
   void global_sum();

   Basic_Seq_Quality_Statistics();
   Basic_Seq_Quality_Statistics( const Basic_Seq_Quality_Statistics& _bsqs);
   ~Basic_Seq_Quality_Statistics();
};

class Output_FQ: public Output_FA {
public:
   Basic_Seq_Quality_Statistics seq_quality_info;
};

class Output_BAM: public Output_FQ {
public:
   int64_t num_primary_alignment = ZeroDefault; // the number of primary alignment/
   int64_t num_secondary_alignment = ZeroDefault; // the number of secondary alignment
   int64_t num_reads_with_secondary_alignment = ZeroDefault; // the number of long reads with the secondary alignment: one read might have multiple seconard alignment
   int64_t num_supplementary_alignment = ZeroDefault; // the number of supplementary alignment
   int64_t num_reads_with_supplementary_alignment = ZeroDefault; // the number of long reads with secondary alignment;
   int64_t num_reads_with_both_secondary_supplementary_alignment = ZeroDefault; // the number of long reads with both secondary and supplementary alignment.

   int64_t forward_alignment = ZeroDefault;
   int64_t reverse_alignment = ZeroDefault;

   //int64_t map_quality_distribution[MAX_MAP_QUALITY]; // The mapping quality distribution;
   std::vector<int64_t>  map_quality_distribution;
   int min_map_quality = MoneDefault; // the minimum mapping quality
   int max_map_quality = MoneDefault; // the maximum mapping quality

   // Similar to Output_FA: below are for mapped.
   //int64_t num_mapped_reads = ZeroDefault; // the number of mapped long reads
   //int64_t num_mapped_bases = ZeroDefault; // the number of mapped bases with X/=/M
   int64_t num_matched_bases = ZeroDefault; // the number of matched bases with =
   int64_t num_mismatched_bases = ZeroDefault; // the number of mismatched bases X
   int64_t num_ins_bases = ZeroDefault; // the number of inserted bases;
   int64_t num_del_bases = ZeroDefault; // the number of deleted bases;
   int64_t num_clip_bases = ZeroDefault; // the number of soft-clipped bases;

   //int64_t accuracy_per_read[PERCENTAGE_ARRAY_SIZE]; // 
   std::vector<int64_t> accuracy_per_read;

   Basic_Seq_Statistics mapped_long_read_info;
   Basic_Seq_Statistics unmapped_long_read_info;
   
   Basic_Seq_Quality_Statistics mapped_seq_quality_info;
   Basic_Seq_Quality_Statistics unmapped_seq_quality_info;

   void reset();
   void add(Output_BAM& t_output_bam); 
   void global_sum();

   Output_BAM();
   ~Output_BAM();
};

class Basic_F5_Statistics{
public: 
   Basic_Seq_Statistics long_read_info;
   Basic_Seq_Quality_Statistics seq_quality_info;
   
   //int signal_range[MAX_SIGNAL_VALUE]; // statistics of all signals
   std::vector<int> signal_range;
   int min_signal = MoneDefault;  // minimum signals;
   int max_signal = MoneDefault;  // maximum signals;

   //int64_t read_quality_distribution[MAX_READ_QUALITY]; // The mapping quality distribution;
   /*std::vector<int64_t> read_quality_distribution;
   int min_read_quality = MoneDefault; // the minimum mapping quality
   int max_read_quality = MoneDefault; // the maximum mapping quality
   */

   // int64_t num_reads = ZeroDefault; // The number of long reads
   // int64_t *read_length_list; // statistics of read length for long reads: each position is the number of reads with the length of the index;

   void reset();
   void add(Basic_F5_Statistics& t_output_bf5);
   void global_sum();

   Basic_F5_Statistics();
   ~Basic_F5_Statistics();
};

class Output_F5 : public Output_Info {
public:
   Basic_F5_Statistics f5_long_read_info;
   Basic_F5_Statistics f5_passed_long_read_info;
   Basic_F5_Statistics f5_failed_long_read_info;

   void reset();
   void add(Output_F5& t_output_f5);
   void global_sum();
};


#endif
