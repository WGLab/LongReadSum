/*

BamReader.cpp:
Class for reading a set number of records from a BAM file. Used for multi-threading.

*/


#include <iostream>
#include <sstream>
#include <fstream>
#include <htslib/sam.h>
#include <math.h>

#include "hts_reader.h"
#include "ComFunction.h"

// HTSReader constructor
HTSReader::HTSReader(const std::string & bam_file_name){
    this->bam_file = hts_open(bam_file_name.c_str(), "r");
    this->header = sam_hdr_read(this->bam_file);
}

// HTSReader destructor
HTSReader::~HTSReader(){
    bam_hdr_destroy(this->header);
    hts_close(this->bam_file);
}

// Update read and base counts
int HTSReader::updateReadAndBaseCounts(bam1_t* record, Basic_Seq_Statistics *basic_qc, std::vector<int> & base_quality_distribution){
    int exit_code = 0;

    // Update the total number of reads
    basic_qc->total_num_reads++;

    // Update read length statistics
    int64_t read_length = (int64_t) record->core.l_qseq;
    basic_qc->total_num_bases += read_length;  // Update the total number of bases
    basic_qc->read_lengths.push_back(read_length);

    // Loop and count the number of each base
    uint8_t *seq = bam_get_seq(record);
    for (int i = 0; i < read_length; i++) {
        // Get the base quality and update the base quality histogram
        int base_quality = (int)bam_get_qual(record)[i];
        base_quality_distribution[base_quality]++;

        // Get the base and update the base count
        char base = seq_nt16_str[bam_seqi(seq, i)];
        switch (base) {
            case 'A':
                basic_qc->total_a_cnt++;
                break;
            case 'C':
                basic_qc->total_c_cnt++;
                break;
            case 'G':
                basic_qc->total_g_cnt++;
                break;
            case 'T':
                basic_qc->total_tu_cnt++;
                break;
            case 'N':
                basic_qc->total_n_cnt++;
                break;
            default:
                std::cerr << "Error reading nucleotide: " << base << std::endl;
                exit_code = 1;
                break;
        }
    }

    return exit_code;
}

// Read the next batch of records from the BAM file and store QC in the output_data object
int HTSReader::readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex){
    int record_count = 0;
    int exit_code = 0;

    // Set up the base quality histogram
    std::vector<int> base_quality_distribution(256, 0);

    // Loop through each alignment record in the BAM file
    // Do QC on each record and store the results in the output_data object
    while ((record_count < batch_size) && (exit_code >= 0)) {
        // Create a record object
        bam1_t* record = bam_init1();

        // read the next record in a thread-safe manner
        read_mutex.lock();
        exit_code = sam_read1(this->bam_file, this->header, record);
        read_mutex.unlock();

        if (exit_code < 0) {
            this->reading_complete = true;
            bam_destroy1(record);
            break; // error or EOF
        }

        // Determine if this is an unmapped read
        if (record->core.flag & BAM_FUNMAP) {
            Basic_Seq_Statistics *basic_qc = &output_data.unmapped_long_read_info;

            // Update read and base QC
            this->updateReadAndBaseCounts(record, basic_qc, base_quality_distribution);

        } else {
            // Set up the basic QC object
            Basic_Seq_Statistics *basic_qc = &output_data.mapped_long_read_info;

            // Calculate base alignment statistics on non-secondary alignments
            if (!(record->core.flag & BAM_FSECONDARY)) {

                // Update the number of mismatches
                uint8_t *nmTag = bam_aux_get(record, "NM");
                output_data.num_mismatched_bases += bam_aux2i(nmTag);

                // Loop through the cigar string and count the number of insertions, deletions, and matches
                uint32_t *cigar = bam_get_cigar(record);
                for (uint32_t i = 0; i < record->core.n_cigar; i++) {
                    int cigar_op = bam_cigar_op(cigar[i]);
                    int cigar_len = bam_cigar_oplen(cigar[i]);
                    switch (cigar_op) {
                        case BAM_CMATCH:
                            output_data.num_matched_bases += cigar_len;
                            break;
                        case BAM_CINS:
                            output_data.num_ins_bases += cigar_len;
                            break;
                        case BAM_CDEL:
                            output_data.num_del_bases += cigar_len;
                            break;
                        case BAM_CSOFT_CLIP:
                            output_data.num_clip_bases += cigar_len;
                            break;
                        case BAM_CHARD_CLIP:
                            output_data.num_clip_bases += cigar_len;
                            break;
                        default:
                            break;
                    }
                }
            }

            // Determine if this is a secondary alignment (not included in QC, only read count)
            if (record->core.flag & BAM_FSECONDARY) {
                output_data.num_secondary_alignment++;

                // Get the alignment's query name (the read name)
                std::string query_name = bam_get_qname(record);

                // Update the read's secondary alignments (count once per read)
                output_data.reads_with_secondary[query_name] = true;

            // Determine if this is a supplementary alignment (not included in QC, only read count)
            } else if (record->core.flag & BAM_FSUPPLEMENTARY) {
                output_data.num_supplementary_alignment++;

                // Get the alignment's query name (the read name)
                std::string query_name = bam_get_qname(record);

                // Update the read's supplementary alignments (count once per read)
                output_data.reads_with_supplementary[query_name] = true;

            // Determine if this is a primary alignment
            } else if (!(record->core.flag & BAM_FSECONDARY || record->core.flag & BAM_FSUPPLEMENTARY)) {
                // Determine if this is a forward or reverse read
                if (record->core.flag & BAM_FREVERSE) {
                    output_data.forward_alignment++;
                } else {
                    output_data.reverse_alignment++;
                }
                output_data.num_primary_alignment++;  // Update the number of primary alignments

                // Update read and base QC
                this->updateReadAndBaseCounts(record, basic_qc, base_quality_distribution);

                // Calculate the percent GC content
                int percent_gc = round((basic_qc->total_g_cnt + basic_qc->total_c_cnt) / (double) (basic_qc->total_a_cnt + basic_qc->total_c_cnt + basic_qc->total_g_cnt + basic_qc->total_tu_cnt) * 100);

                // Update the GC content histogram
                basic_qc->read_gc_content_count.push_back(percent_gc);
            } else {
                std::cerr << "Error: Unknown alignment type" << std::endl;
                std::cerr << "Flag: " << record->core.flag << std::endl;
            }
        }

        // Delete the record object
        bam_destroy1(record);

        record_count++;
    }

    // Update the base quality histogram
    output_data.seq_quality_info.base_quality_distribution = base_quality_distribution;

    return exit_code;
}

// Return if the file has any more records
bool HTSReader::hasNextRecord(){
    return !this->reading_complete;
}

// Return the number of records in the BAM file using the BAM index
int64_t HTSReader::getNumRecords(const std::string & bam_filename){
    samFile* bam_file = sam_open(bam_filename.c_str(), "r");
    bam_hdr_t* bam_header = sam_hdr_read(bam_file);
    bam1_t* bam_record = bam_init1();

    int64_t num_reads = 0;
    while (sam_read1(bam_file, bam_header, bam_record) >= 0) {
        num_reads++;
    }

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_file);

    std::cout << "Number of reads in " << bam_filename << ": " << num_reads << std::endl;

    return num_reads;
}
