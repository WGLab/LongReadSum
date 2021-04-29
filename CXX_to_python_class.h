#ifndef CXX_TO_PYTHON_CLASS_H
#define CXX_TO_PYTHON_CLASS_H
// CXX_to_python_class

#include <string>
#include <vector>

//                        ,  ,
#define MAX_READ_LENGTH 10485760 
#define MAX_MAP_QUALITY 256

#define MAX_SIGNAL_VALUE 5000

class Output_FA {
public:
   int64_t total_num_reads;
   int64_t total_num_bases;

   int64_t longest_read_length;
   int64_t n50_read_length;
   int64_t n95_read_length;
   int64_t n05_read_length;
   int64_t mean_read_length;
   int64_t median_read_length;

   int64_t total_a_cnt;
   int64_t total_c_cnt;
   int64_t total_g_cnt;
   int64_t total_tu_cnt;
   int64_t total_n_cnt;
   double gc_cnt;

   int64_t *read_length_count;

   Output_FA();
   ~Output_FA();
};

class Output_FQ: public Output_FA {
public:
   
};

class Output_BAM: public Output_FQ {
public:
   int64_t num_primary_alignment;
   int64_t num_secondary_alignment;
   int64_t num_reads_with_secondary_alignment;
   int64_t num_supplementary_alignment;
   int64_t num_reads_with_supplementary_alignment;
   int64_t num_reads_with_both_secondary_supplementary_alignment;

   int64_t map_quality_distribution[MAX_MAP_QUALITY];

   int64_t num_mapped_reads;
   int64_t num_mapped_bases;
   int64_t num_matched_bases;
   int64_t num_mismatched_bases;
   int64_t num_ins_bases;
   int64_t num_del_bases;

   // similar to reads; below for mapped reads
   int64_t mapped_longest_read_length;
   int64_t mapped_n50_read_length;
   int64_t mapped_n95_read_length;
   int64_t mapped_n05_read_length;
   int64_t mapped_mean_read_length;
   int64_t mapped_median_read_length;

   int64_t mapped__a_cnt;
   int64_t mapped__c_cnt;
   int64_t mapped__g_cnt;
   int64_t mapped__tu_cnt;
   int64_t mapped__n_cnt;
   double mapped_gc_cnt;

   int64_t *mapped_read_length_count;

   Output_BAM();
   ~Output_BAM();
};

class Output_F5 {
public:
   Output_BAM bam_output;

   int64_t num_passed_reads;
   int64_t num_failed_reads;
   
   int64_t *passed_reads_list;
   int64_t *failed_reads_list;

   int signal_range[MAX_SIGNAL_VALUE];
   int min_signal;
   int max_signal;

   Output_F5();
   ~Output_F5(); 
};


#endif
