#ifndef CXX_TO_PYTHON_CLASS_H
#define CXX_TO_PYTHON_CLASS_H
// CXX_to_python_class

#include <string>

//                        ,  ,
#define MAX_READ_LENGTH 10485760 
#define MAX_MAP_QUALITY 256
#define MAX_BASE_QUALITY 256
#define MAX_SIGNAL_VALUE 5000

#define ZeroDefault 0
#define MoneDefault -1

class Output_FA {
public:
   int64_t total_num_reads = ZeroDefault; // total number of long reads
   int64_t total_num_bases = ZeroDefault; // total number of bases

   int64_t longest_read_length = MoneDefault; // the length of longest reads
   int64_t n50_read_length = MoneDefault; // N50
   int64_t n95_read_length = MoneDefault; // N95
   int64_t n05_read_length = MoneDefault; // N05;
   int64_t mean_read_length = MoneDefault; // mean of read length
   int64_t nx_read_length[10]; // can provide n10, n20, n30, ... n90; The last indicates whether this info is generated.
   int64_t median_read_length = MoneDefault; // median of read length

   int64_t total_a_cnt = ZeroDefault; // A content
   int64_t total_c_cnt = ZeroDefault; // C content
   int64_t total_g_cnt = ZeroDefault; // G content
   int64_t total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
   int64_t total_n_cnt = ZeroDefault; // N content
   double gc_cnt = ZeroDefault; // GC ratio

   int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

   Output_FA();
   ~Output_FA();
};

class Output_FQ: public Output_FA {
public:
   int64_t base_quality_distribution[MAX_BASE_QUALITY];  // base quality distribution: each position is the number of bases with the quality of the index.
   int min_base_quality = MoneDefault; // minimum base quality;
   int max_base_quality = MoneDefault; // maximum base quality;

   int *pos_quality_distribution; // the base quality according to positions of long reads: each position is the mean quality at the position of long reads
   float *pos_quality_distribution_dev; // similar to above but for the standard deviation
   int64_t *pos_quality_distribution_count; // similar to above but for the number of bases

   Output_FQ();
   ~Output_FQ();
};

class Output_BAM: public Output_FQ {
public:
   int64_t num_primary_alignment = ZeroDefault; // the number of primary alignment/
   int64_t num_secondary_alignment = ZeroDefault; // the number of secondary alignment
   int64_t num_reads_with_secondary_alignment = ZeroDefault; // the number of long reads with the secondary alignment: one read might have multiple seconard alignment
   int64_t num_supplementary_alignment = ZeroDefault; // the number of supplementary alignment
   int64_t num_reads_with_supplementary_alignment = ZeroDefault; // the number of long reads with secondary alignment;
   int64_t num_reads_with_both_secondary_supplementary_alignment = ZeroDefault; // the number of long reads with both secondary and supplementary alignment.

   int64_t map_quality_distribution[MAX_MAP_QUALITY]; // The mapping quality distribution;
   int min_map_quality = MoneDefault; // the minimum mapping quality
   int max_map_quality = MoneDefault; // the maximum mapping quality

   // Similar to Output_FA: below are for mapped.
   int64_t num_mapped_reads = ZeroDefault; // the number of mapped long reads
   int64_t num_mapped_bases = ZeroDefault; // the number of mapped bases with X/=/M
   int64_t num_matched_bases = ZeroDefault; // the number of matched bases with =
   int64_t num_mismatched_bases = ZeroDefault; // the number of mismatched bases X
   int64_t num_ins_bases = ZeroDefault; // the number of inserted bases;
   int64_t num_del_bases = ZeroDefault; // the number of deleted bases;
   int64_t num_clip_bases = ZeroDefault; // the number of soft-clipped bases;
   
   // similar to reads; below for mapped reads
   int64_t mapped_longest_read_length = MoneDefault; // the length of longest reads
   int64_t mapped_n50_read_length = MoneDefault; // N50
   int64_t mapped_n95_read_length = MoneDefault; // N95;
   int64_t mapped_n05_read_length = MoneDefault; // N05;
   int64_t nmapped_x_read_length[10]; // can provide n10, n20, n30, ... n90; The last indicates whether this info is generated.
   int64_t mapped_mean_read_length = MoneDefault; // mean of read length
   int64_t mapped_median_read_length = MoneDefault; // median of read length

   int64_t mapped__a_cnt = ZeroDefault; // Mapped A content;
   int64_t mapped__c_cnt = ZeroDefault; // Napped C content;
   int64_t mapped__g_cnt = ZeroDefault; // Mapped G content;
   int64_t mapped__tu_cnt = ZeroDefault; // Mapped T content or U content;
   int64_t mapped__n_cnt = ZeroDefault; // Mapped N content;
   double mapped_gc_cnt = ZeroDefault;  // GC content for mapped bases

   int64_t *mapped_read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

   Output_BAM();
   ~Output_BAM();
};

class Output_F5 {
public:
   Output_BAM bam_output;

   int64_t num_passed_reads = ZeroDefault; // The number of passed long reads
   int64_t num_failed_reads = ZeroDefault; // The number of failed long reads
   
   int64_t *passed_read_length_list; // statistics of read length for passed long reads: each position is the number of reads with the length of the index;
   int64_t *failed_read_length_list; // statistics of read length for failed long reads: each position is the number of reads with the length of the index;

   int signal_range[MAX_SIGNAL_VALUE]; // statistics of all signals
   int min_signal = MoneDefault;  // minimum signals;
   int max_signal = MoneDefault;  // maximum signals;

   Output_F5();
   ~Output_F5(); 
};


#endif
