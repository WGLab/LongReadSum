#ifndef HTSREADER_H_
#define HTSREADER_H_

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <mutex>
#include <stdint.h>

#include "output_data.h"

#define BAM_UN_OPEN 1
#define BAM_FAILED 2
#define BAM_OPEN 3
#define BAM_CLOSE 4

#define DNASEQ 1
#define RNASEQ 2

// Define the modified base map type (position -> mod_type, canonical_base, likelihood)
typedef std::map<int32_t, std::tuple<char, char, double>> Base_Modification_Map;

// Class representing an HTSlib BAM file
class HTSReader {
    public:
        htsFile* bam_file; // open the BAM file for reading
        bam_hdr_t* header; // read the BAM header
        bam1_t* record;
        int record_count = 0;

        // Bool for whether the reading is complete
        bool reading_complete = false;

        // Update read and base counts
        int updateReadAndBaseCounts(bam1_t* record, Basic_Seq_Statistics* basic_qc, uint64_t *base_quality_distribution);

        // Read the next batch of records from the BAM file
        int readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex, std::unordered_set<std::string>& read_ids);

        // Return if the reader has finished reading the BAM file
        bool hasNextRecord();

        // Return the number of records in the BAM file using the BAM index
        int64_t getNumRecords(const std::string &bam_file_name);

        // Input a query position and return the reference position by querying
        // the BAM index
        int getRefPos(bam1_t* record, std::vector<int> query_pos, std::vector<int> &ref_pos);
        // int getRefPos(const std::string &bam_file_name, const std::string &chromosome, std::vector<int> query_pos, std::vector<int> &ref_pos);

        // Add a modification to the base modification map
        void addModificationToQueryMap(Base_Modification_Map &base_modifications, int32_t pos, char mod_type, char canonical_base, double likelihood, bool is_cpg);

        HTSReader(const std::string &bam_file_name);
        ~HTSReader();
};

#endif


