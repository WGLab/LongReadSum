#ifndef BAMREADER_H_
#define BAMREADER_H_

#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include <iostream>
#include <sstream>
#include <fstream>

#include <mutex>

#include "ComStruct.h"
#include "OutputStructures.h"

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

        // Bool for whether the reading is complete
        bool reading_complete = false;

        // Read the next batch of records from the BAM file
        int readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex);

        // Return if the reader has finished reading the BAM file
        bool hasNextRecord();

        HTSReader(const std::string &bam_file_name);
        ~HTSReader();
};

#endif


