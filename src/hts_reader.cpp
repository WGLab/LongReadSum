/*

BamReader.cpp:
Class for reading a set number of records from a BAM file. Used for multi-threading.

*/

#include "hts_reader.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <algorithm>  // std::find
#include <htslib/sam.h>

#include "utils.h"

void HTSReader::addModificationToQueryMap(std::map<int32_t, std::tuple<char, char, double, int>> &base_modifications, int32_t pos, char mod_type, char canonical_base, double likelihood, int strand)
{
    // Add the modification type to the map
    base_modifications[pos] = std::make_tuple(mod_type, canonical_base, likelihood, strand);
}

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
int HTSReader::readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex, std::unordered_set<std::string>& read_ids, double base_mod_threshold) {
    int record_count = 0;
    int exit_code = 0;

    // Determine if filtering by read ID
    bool read_ids_present = false;
    if (read_ids.size() > 0){
        read_ids_present = true;
        printMessage("Filtering reads by read ID");

        printMessage("Number of read IDs: " + std::to_string(read_ids.size()));
        printMessage("First read ID: " + *read_ids.begin());
        // Check if the first read ID has any newlines, carriage returns, tabs,
        // or spaces
        if (read_ids.begin()->find_first_of("\n\r\t ") != std::string::npos) {
            printError("Read IDs cannot contain newlines, carriage returns, tabs, or spaces");
            return 1;
        }
    }

    // Access the base quality histogram from the output_data object
    uint64_t *base_quality_distribution = output_data.seq_quality_info.base_quality_distribution;

    // Do QC on each record and store the results in the output_data object
    while ((record_count < batch_size) && (exit_code >= 0)) {
        // Create a record object
        bam1_t* record = bam_init1();

        // Read the next record in a thread-safe manner
        read_mutex.lock();
        exit_code = sam_read1(this->bam_file, this->header, record);
        read_mutex.unlock();

        if (exit_code < 0) {
            this->reading_complete = true;
            bam_destroy1(record);
            break; // error or EOF
        }


        // Get the read (query) name
        std::string query_name = bam_get_qname(record);

        // Check if the read name has any newlines, carriage returns, tabs, or
        // spaces
        if (query_name.find_first_of("\n\r\t ") != std::string::npos) {
            printError("BAM Read names cannot contain newlines, carriage returns, tabs, or spaces");
            return 1;
        }

        // Determine if this read should be skipped
        if (read_ids_present) {

            // Determine if this read should be skipped
            if (read_ids.find(query_name) == read_ids.end()){
                continue;  // Skip this read
            }
        }

        // For POD5 files, corresponding BAM files will have tags for
        // indexing signal data for each base (ts, ns, mv). Find the
        // tags and store in the output data object
        uint8_t *ts_tag = bam_aux_get(record, "ts");
        uint8_t *ns_tag = bam_aux_get(record, "ns");
        uint8_t *mv_tag = bam_aux_get(record, "mv");

        // Get POD5 signal tag values if they exist
        if (mv_tag != NULL && ts_tag != NULL && ns_tag != NULL) {
            // Set the atomic flag and print a message if the POD5 tags are
            // present
            if (!this->has_pod5_tags.test_and_set()) {
                printMessage("POD5 tags found (ts, ns, mv)");
            }

            // Get the ts and ns tags
            int32_t ts = bam_aux2i(ts_tag);
            int32_t ns = bam_aux2i(ns_tag);

            // Get the move table (start at 1 to skip the tag type)
            int32_t length = bam_auxB_len(mv_tag);
            std::vector<int32_t> move_table(length);
            std::vector<std::vector<int>> sequence_move_table;  // Store the sequence move table with indices
            
            // Iterate over the move table values
            int current_index = ts;
            std::vector<int> signal_index_vector;
            int move_value = 0;
            for (int32_t i = 1; i < length; i++) {
                move_value = bam_auxB2i(mv_tag, i);
                if (move_value == 1) {
                    signal_index_vector.push_back(current_index);
                }

                current_index++;
            }
            // Create a tuple and add the read's signal data to the output data
            std::string seq_str = "";
            for (int i = 0; i < record->core.l_qseq; i++) {
                seq_str += seq_nt16_str[bam_seqi(bam_get_seq(record), i)];
            }

            // Throw an error if the query name is empty
            if (query_name.empty()) {
                std::cerr << "Error: Query name is empty" << std::endl;
                exit_code = 1;
                break;
            }
            output_data.addReadMoveTable(query_name, seq_str, signal_index_vector, ts, ns);
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

                    // Set the atomic flag and print a message if the NM tag is
                    // present
                    if (!this->has_nm_tag.test_and_set()) {
                        printMessage("NM tag found, used NM tag for mismatch count");
                    }
                } else {
                    // Set the atomic flag and print a message if the NM tag is
                    // not present
                    if (!this->has_nm_tag.test_and_set()) {
                        printMessage("No NM tag found, using CIGAR for mismatch count");
                    }
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

                // Loop through the cigar string, count the number of clipped
                // bases, and also get the reference position of the read if
                // there is a modification tag
                uint32_t *cigar = bam_get_cigar(record);
                int32_t ref_pos = record->core.pos;
                int query_pos = 0;
                for (uint32_t i = 0; i < record->core.n_cigar; i++) {
                    int cigar_op = bam_cigar_op(cigar[i]);
                    uint64_t cigar_len = (uint64_t)bam_cigar_oplen(cigar[i]);
                    switch (cigar_op) {
                        case BAM_CSOFT_CLIP:
                            output_data.num_clip_bases += cigar_len;
                            query_pos += cigar_len; // Consumes query bases
                            break;
                        case BAM_CHARD_CLIP:
                            output_data.num_clip_bases += cigar_len;
                            break;
                        case BAM_CMATCH:
                        case BAM_CEQUAL:
                        case BAM_CDIFF:
                            ref_pos += cigar_len;
                            query_pos += cigar_len;
                            break;
                        case BAM_CINS:
                            query_pos += cigar_len;
                            break;
                        case BAM_CDEL:
                            ref_pos += cigar_len;
                            break;
                        case BAM_CREF_SKIP:
                            ref_pos += cigar_len;
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

            } else {
                std::cerr << "Error: Unknown alignment type" << std::endl;
                std::cerr << "Flag: " << record->core.flag << std::endl;
            }
        }

        // Delete the record object
        bam_destroy1(record);

        record_count++;
    }

    return exit_code;
}

// Return if the file has any more records
bool HTSReader::hasNextRecord(){
    return !this->reading_complete;
}

// Return the number of records in the BAM file using the BAM index
int64_t HTSReader::getNumRecords(const std::string & bam_filename, Output_BAM &final_output, bool mod_analysis, double base_mod_threshold) {
    samFile* bam_file = sam_open(bam_filename.c_str(), "r");
    bam_hdr_t* bam_header = sam_hdr_read(bam_file);
    bam1_t* bam_record = bam_init1();

    int64_t num_reads = 0;
    while (sam_read1(bam_file, bam_header, bam_record) >= 0) {
        num_reads++;

        if (mod_analysis) {

            // Base modification tag analysis
            // Follow here to get base modification tags:
            // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/sam_mods.c
            // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/htslib/sam.h#L2274
            hts_base_mod_state *state = hts_base_mod_state_alloc();

            // Preprint revisions: New data structure that does not require unique
            // positions for each base modification
            // chr -> vector of (position, strand) for C modified bases passing the threshold
            std::vector<std::pair<int32_t, int>> c_modified_positions;

            // Parse the base modification tags if a primary alignment
            int ret = bam_parse_basemod(bam_record, state);
            if (ret >= 0 && !(bam_record->core.flag & BAM_FSECONDARY) && !(bam_record->core.flag & BAM_FSUPPLEMENTARY) && !(bam_record->core.flag & BAM_FUNMAP)) {

                // Get the chromosome if alignments are present
                bool alignments_present = true;
                std::string chr;
                std::map<int, int> query_to_ref_map;
                if (bam_record->core.tid < 0) {
                    alignments_present = false;
                } else {
                    chr = bam_header->target_name[bam_record->core.tid];

                    // Get the query to reference position mapping
                    query_to_ref_map = this->getQueryToRefMap(bam_record);
                }

                // Get the strand from the alignment flag (hts_base_mod uses 0 for positive and 1 for negative,
                // but it always yields 0...)
                int strand = (bam_record->core.flag & BAM_FREVERSE) ? 1 : 0;

                // Iterate over the state object to get the base modification tags
                // using bam_next_basemod
                hts_base_mod mods[10];
                int n = 0;
                int32_t pos = 0;
                std::vector<int> query_pos;
                while ((n=bam_next_basemod(bam_record, state, mods, 10, &pos)) > 0) {
                    for (int i = 0; i < n; i++) {
                        // Update the prediction count
                        final_output.modified_prediction_count++;

                        // Note: The modified base value can be a positive char (e.g. 'm',
                        // 'h') (DNA Mods DB) or negative integer (ChEBI ID):
                        // https://github.com/samtools/hts-specs/issues/741
                        // DNA Mods: https://dnamod.hoffmanlab.org/
                        // ChEBI: https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:21839
                        // Header line:
                        // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/htslib/sam.h#L2215

                        // Determine the probability of the modification (-1 if
                        // unknown)
                        double probability = -1;
                        if (mods[i].qual != -1) {
                            probability = mods[i].qual / 256.0;

                            // If the probability is greater than the threshold,
                            // update the count
                            if (probability >= base_mod_threshold) {
                                final_output.sample_modified_base_count++;

                                // Update the modified base count for the strand
                                if (strand == 0) {
                                    final_output.sample_modified_base_count_forward++;
                                } else {
                                    final_output.sample_modified_base_count_reverse++;
                                }

                                // Preprint revisions: Store the modified positions
                                char canonical_base_char = std::toupper(mods[i].canonical_base);
                                char mod_type = mods[i].modified_base;
                                if (canonical_base_char == 'C' && mod_type == 'm') {

                                    // Convert the query position to reference position if available
                                    if (alignments_present) {
                                        if (query_to_ref_map.find(pos) != query_to_ref_map.end()) {
                                            int32_t ref_pos = query_to_ref_map[pos];
                                            c_modified_positions.push_back(std::make_pair(ref_pos, strand));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // Preprint revisions: Append the modified positions to the output data
                if (c_modified_positions.size() > 0) {
                    // Set the atomic flag and print a message if base
                    // modification tags are present in the file
                    if (!this->has_mm_ml_tags.test_and_set()) {
                        printMessage("Base modification data found (MM, ML tags)");
                    }

                    // Add the modified positions to the output data
                    if (final_output.sample_c_modified_positions.find(chr) == final_output.sample_c_modified_positions.end()) {
                        final_output.sample_c_modified_positions[chr] = c_modified_positions;
                    } else {
                        final_output.sample_c_modified_positions[chr].insert(final_output.sample_c_modified_positions[chr].end(), c_modified_positions.begin(), c_modified_positions.end());
                    }
                }
            }

            // Deallocate the state object
            hts_base_mod_state_free(state);
        }
    }

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_file);

    return num_reads;
}

// Get the mapping of query positions to reference positions for a given alignment record
std::map<int, int> HTSReader::getQueryToRefMap(bam1_t *record)
{
    std::map<int, int> query_to_ref_map;

    // Initialize the starting reference and query positions
    int32_t current_ref_pos = record->core.pos;  // Get the reference position
    int current_query_pos = 0;
    uint32_t *cigar = bam_get_cigar(record);

    // Iterate over the CIGAR operations
    int cigar_len = record->core.n_cigar;
    for (int i = 0; i < cigar_len; i++) {
        int cigar_op = bam_cigar_op(cigar[i]);  // Get the CIGAR operation
        int op_len = bam_cigar_oplen(cigar[i]);  // Get the CIGAR operation length
        
        switch (cigar_op) {
            case BAM_CDIFF:
                current_ref_pos += op_len;
                current_query_pos += op_len;
                break;
            case BAM_CMATCH:
            case BAM_CEQUAL:
                for (int j = 0; j < op_len; j++) {
                    query_to_ref_map[current_query_pos] = current_ref_pos + 1;  // Use 1-indexed positions
                    current_ref_pos++;
                    current_query_pos++;
                    // query_to_ref_map[current_query_pos] = current_ref_pos + 1;  // Use 1-indexed positions
                }
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                current_query_pos += op_len;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                current_ref_pos += op_len;
                break;
            default:
                // Handle unexpected CIGAR operations if needed
                break;
        }
    }

    return query_to_ref_map;
}
