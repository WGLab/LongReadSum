#ifndef HTSREADER_H_
#define HTSREADER_H_

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <mutex>
#include <stdint.h>
#include <atomic>

#include "output_data.h"

#define BAM_UN_OPEN 1
#define BAM_FAILED 2
#define BAM_OPEN 3
#define BAM_CLOSE 4

#define DNASEQ 1
#define RNASEQ 2

// Class representing an HTSlib BAM file
class HTSReader {
    public:
        htsFile* bam_file; // open the BAM file for reading
        bam_hdr_t* header; // read the BAM header
        bam1_t* record;
        int record_count = 0;
        
        // Atomic flags for whether certain BAM flags are present
        std::atomic_flag has_nm_tag = ATOMIC_FLAG_INIT;  // NM tag for number of mismatches using edit distance
        std::atomic_flag has_mm_ml_tags = ATOMIC_FLAG_INIT;  // MM and ML tags for modified base information
        std::atomic_flag has_pod5_tags = ATOMIC_FLAG_INIT;  // POD5 tags for signal-level information (ts, ns)

        // Bool for whether the reading is complete
        bool reading_complete = false;

        // Update read and base counts
        int updateReadAndBaseCounts(bam1_t* record, Basic_Seq_Statistics& basic_qc, Basic_Seq_Quality_Statistics& seq_quality_info, bool is_primary);

        // Read the next batch of records from the BAM file
        int readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex, std::unordered_set<std::string>& read_ids, double base_mod_threshold);

        // Return if the reader has finished reading the BAM file
        bool hasNextRecord();

        // Return the number of records in the BAM file using the BAM index
        int getNumRecords(const std::string &bam_file_name, int thread_count);

        // Run base modification analysis
        void runBaseModificationAnalysis(const std::string &bam_filename, Output_BAM& final_output, double base_mod_threshold, int read_count, int sample_count, int thread_count);

        std::map<int, int> getQueryToRefMap(bam1_t* record);

        // Add a modification to the base modification map
        void addModificationToQueryMap(std::map<int32_t, std::tuple<char, char, double, int>> &base_modifications, int32_t pos, char mod_type, char canonical_base, double likelihood, int strand);

        HTSReader(const std::string &bam_file_name);
        ~HTSReader();
};

#endif


