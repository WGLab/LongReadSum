/*

BamReader.cpp:
Class for reading a set number of records from a BAM file. Used for multi-threading.

*/


#include <iostream>
#include <sstream>
#include <fstream>
#include <htslib/sam.h>
#include <math.h>
#include <algorithm>  // std::find

#include "hts_reader.h"
#include "utils.h"

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
int HTSReader::updateReadAndBaseCounts(bam1_t* record, Basic_Seq_Statistics *basic_qc, uint64_t *base_quality_distribution){
    int exit_code = 0;

    // Update the total number of reads
    basic_qc->total_num_reads++;

    // Update read length statistics
    int read_length = (int) record->core.l_qseq;
    basic_qc->total_num_bases += (uint64_t) read_length;  // Update the total number of bases
    basic_qc->read_lengths.push_back(read_length);

    // Loop and count the number of each base
    uint8_t *seq = bam_get_seq(record);
    for (int i = 0; i < read_length; i++) {
        // Get the base quality and update the base quality histogram
        uint64_t base_quality = (uint64_t)bam_get_qual(record)[i];
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
                std::cerr << "Warning: N base found in read " << bam_get_qname(record) << std::endl;
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
int HTSReader::readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex, std::unordered_set<std::string>& read_ids){
    int record_count = 0;
    int exit_code = 0;

    // Determine if filtering by read ID
    bool read_ids_present = false;
    if (read_ids.size() > 0){
        read_ids_present = true;
    }

    // Access the base quality histogram from the output_data object
    uint64_t *base_quality_distribution = output_data.seq_quality_info.base_quality_distribution;

    // Loop through each alignment record in the BAM file
    // Do QC on each record and store the results in the output_data object
    // bool nm_tag_present = false;  // Flag to determine if the NM tag is present (for mismatch counting)
    bool mod_tag_present = false;  // Flag to determine if the base modification tags (MM, ML) are present
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

        // Follow here to get base modification tags:
        // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/sam_mods.c
        // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/htslib/sam.h#L2274
        hts_base_mod_state *state = hts_base_mod_state_alloc();
        if (bam_parse_basemod(record, state) >= 0) {
            printMessage("Base modification tags found");
            // std::cout << "Base modification tags found" << std::endl;
            mod_tag_present = true;

            // Iterate over the state object to get the base modification tags
            // using bam_next_basemod
            hts_base_mod mods[10];
            int n = 0;
            int pos = 0;
            while ((n=bam_next_basemod(record, state, mods, 10, &pos)) > 0) {
                for (int i = 0; i < n; i++) {
                    // Struct definition: https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/htslib/sam.h#L2226
                    printMessage("Found base modification at position " + std::to_string(pos));
                    printMessage("Modification type: " + std::string(1, mods[i].modified_base));
                    printMessage("Canonical base: " + std::string(1, mods[i].canonical_base));
                    printMessage("Likelihood: " + std::to_string(mods[i].qual / 256.0));
                    printMessage("Strand: " + std::to_string(mods[i].strand));

                    //
                    // std::cout << "Base modification at position " << pos << std::endl;
                    // std::cout << "Base modification type: " << mods[i].modified_base << std::endl;
                    // std::cout << "Base modification likelihood: " << mods[i].qual / 256.0 << std::endl;
                    // std::cout << "Base modification strand: " << mods[i].strand << std::endl;
                }
            }

            // Iterating by position
            // hts_base_mod mods[10];
            // int n = bam_mods_at_next_pos(record, state, mods, 10);
            // for (int i = 0; i < n; i++) {
            //     std::cout << "Base modification at position " << mods[i].pos << std::endl;
            //     std::cout << "Base modification type: " << mods[i].type << std::endl;
            //     std::cout << "Base modification likelihood: " << mods[i].likelihood << std::endl;
            // }
            

            // // Get the ML tag (base modification likelihoods) from the state
            // // object (https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/sam_mods.c#L176)
            // if (state->ml) {
            //     std::cout << "ML tag found" << std::endl;
            //     // std::cout << "ML tag: " << state->ml << std::endl;
            //     // std::cout << "ML tag length: " << state->ml_len << std::endl;
            //     // std::cout << "ML tag type: " << state->ml_type << std::endl;
            //     // std::cout << "ML tag num: " << state->ml_num << std::endl;
            //     // std::cout << "ML tag num length: " << state->ml_num_len << std::endl;
            //     // std::cout << "ML tag num type: " << state->ml_num_type <<
            //     // std::endl;
            // }
        } else {
            std::cout << "No base modification tags found" << std::endl;
        }

        // Determine if this read should be skipped
        if (read_ids_present){
            // Get the alignment's query name (the read name)
            std::string query_name = bam_get_qname(record);
            //std::cout << "Query name: " << query_name << std::endl;

            // Determine if this read should be skipped
            if (read_ids.find(query_name) == read_ids.end()){
                // std::cout << "Skipping read " << query_name << std::endl;
                continue;  // Skip this read
            }
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

                // Determine if this is a forward or reverse read
                if (record->core.flag & BAM_FREVERSE) {
                    output_data.forward_alignment++;
                } else {
                    output_data.reverse_alignment++;
                }

                // Loop through the cigar string and count the number of insertions, deletions, and matches
                uint32_t *cigar = bam_get_cigar(record);
                uint64_t num_mismatches = 0;
                for (uint32_t i = 0; i < record->core.n_cigar; i++) {
                    int cigar_op = bam_cigar_op(cigar[i]);
                    uint64_t cigar_len = (uint64_t)bam_cigar_oplen(cigar[i]);
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
                        case BAM_CDIFF:
                            num_mismatches += cigar_len;
                            break;
                        default:
                            break;
                    }
                }

                // Update the number of mismatches if the NM tag is present (more accurate than CIGAR)
                uint8_t *nmTag = bam_aux_get(record, "NM");
                if (nmTag != NULL) {
                    num_mismatches = (uint64_t) bam_aux2i(nmTag);
                    // nm_tag_present = true;
                }
                output_data.num_mismatched_bases += num_mismatches;
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
                output_data.num_primary_alignment++;  // Update the number of primary alignments

                // Loop through the cigar string and count the number of clipped bases
                uint32_t *cigar = bam_get_cigar(record);
                for (uint32_t i = 0; i < record->core.n_cigar; i++) {
                    int cigar_op = bam_cigar_op(cigar[i]);
                    uint64_t cigar_len = (uint64_t)bam_cigar_oplen(cigar[i]);
                    switch (cigar_op) {
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

                // Update read and base QC
                this->updateReadAndBaseCounts(record, basic_qc, base_quality_distribution);

                // Calculate the percent GC content
                int percent_gc = round((basic_qc->total_g_cnt + basic_qc->total_c_cnt) / (double) (basic_qc->total_a_cnt + basic_qc->total_c_cnt + basic_qc->total_g_cnt + basic_qc->total_tu_cnt) * 100);

                // Update the GC content histogram
                basic_qc->read_gc_content_count.push_back(percent_gc);

                // Determine if the base modification tags are present
                uint8_t *mmTag = bam_aux_get(record, "MM:Z");
                uint8_t *mlTag = bam_aux_get(record, "ML:B:C");
                // uint8_t *mmTag = bam_aux_get(record, "mm");
                // uint8_t *mlTag = bam_aux_get(record, "ml");
                // uint8_t *mmTag = bam_aux_get(record, "MM");
                // uint8_t *mlTag = bam_aux_get(record, "ML");

                if (mmTag != NULL || mlTag != NULL) {
                    mod_tag_present = true;
                }

            } else {
                std::cerr << "Error: Unknown alignment type" << std::endl;
                std::cerr << "Flag: " << record->core.flag << std::endl;
            }
        }

        // Delete the record object
        bam_destroy1(record);

        record_count++;
    }

    // Print if the NM tag was not present
    // if (nm_tag_present)
    // {
    //     std::cout << "NM tag found, used NM tag for mismatch count" << std::endl;
    // } else {
    //     std::cout << "No NM tag found, used CIGAR for mismatch count" << std::endl;
    // }

    // Print if the base modification tags were present
    if (mod_tag_present)
    {
        std::cout << "Base modification tags found" << std::endl;
    } else {
        std::cout << "[test2] No base modification tags found" << std::endl;
    }

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

    return num_reads;
}
