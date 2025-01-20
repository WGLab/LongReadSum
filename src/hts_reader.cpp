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
#include <random>
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
int HTSReader::updateReadAndBaseCounts(bam1_t* record, Basic_Seq_Statistics& basic_qc, Basic_Seq_Quality_Statistics& seq_quality_info, bool is_primary) {

    // Update read QC
    basic_qc.total_num_reads++;  // Update the total number of reads
    int read_length = (int) record->core.l_qseq;
    basic_qc.total_num_bases += (uint64_t) read_length;  // Update the total number of bases
    basic_qc.read_lengths.push_back(read_length);

    // Get base counts, quality, and GC content
    double read_gc_count = 0.0;  // For GC content calculation
    double read_base_total = 0.0;  // For GC content calculation
    double cumulative_base_prob = 0.0;  // For mean base quality probability calculation
    uint8_t *seq = bam_get_seq(record);
    for (int i = 0; i < read_length; i++) {
        // Get the base quality and update the base quality histogram
        int base_quality = (int)bam_get_qual(record)[i];
        seq_quality_info.base_quality_distribution[(uint64_t)base_quality]++;

        // Convert the Phred quality value to a probability
        double base_quality_prob = pow(10, -base_quality / 10.0);
        cumulative_base_prob += base_quality_prob;

        // Get the base and update the base count
        char base = seq_nt16_str[bam_seqi(seq, i)];
        switch (base) {
            case 'A':
                basic_qc.total_a_cnt++;
                read_base_total++;
                break;
            case 'C':
                basic_qc.total_c_cnt++;
                read_gc_count++;
                read_base_total++;
                break;
            case 'G':
                basic_qc.total_g_cnt++;
                read_gc_count++;
                read_base_total++;
                break;
            case 'T':
                basic_qc.total_tu_cnt++;
                read_base_total++;
                break;
            case 'N':
                basic_qc.total_n_cnt++;
                std::cerr << "Warning: N base found in read " << bam_get_qname(record) << std::endl;
                break;
            default:
                printError("Invalid base: " + std::to_string(base));
                break;
        }
    }

    // Calculate the mean base quality probability
    cumulative_base_prob /= (double)read_length;

    // Convert the mean base quality probability to a Phred quality value
    double read_mean_base_qual = -10.0 * log10(cumulative_base_prob);

    // Update the per-read mean base quality distribution
    int read_mean_base_qual_int = static_cast<int>(std::round(read_mean_base_qual));
    try {
        seq_quality_info.read_average_base_quality_distribution[read_mean_base_qual_int]++;
    } catch (const std::out_of_range& oor) {
        printError("Warning: Base quality value " + std::to_string(read_mean_base_qual_int) + " exceeds maximum value");
    }

    // Calculate the read GC content percentage if a primary alignment
    if (is_primary) {
        double gc_content = read_gc_count / read_base_total;
        int gc_content_percent = (int) round(gc_content * 100);
        std::string query_name = bam_get_qname(record);
        // printMessage("Read name: " + query_name + ", GC content: " + std::to_string(gc_content) + ", GC count: " + std::to_string(read_gc_count) + ", Total count: " + std::to_string(read_base_total));
        basic_qc.read_gc_content_count[gc_content_percent]++;
    }

    return 0;
}

// Read the next batch of records from the BAM file and store QC in the output_data object
int HTSReader::readNextRecords(int batch_size, Output_BAM & output_data, std::mutex & read_mutex, std::unordered_set<std::string>& read_ids, double base_mod_threshold) {
    int record_count = 0;
    int exit_code = 0;

    // Determine if filtering by read ID
    bool read_ids_present = false;
    if (read_ids.size() > 0){
        read_ids_present = true;
        // printMessage("Filtering reads by read ID");

        // printMessage("Number of read IDs: " + std::to_string(read_ids.size()));
        // printMessage("First read ID: " + *read_ids.begin());
        // Check if the first read ID has any newlines, carriage returns, tabs,
        // or spaces
        if (read_ids.begin()->find_first_of("\n\r\t ") != std::string::npos) {
            printError("Read IDs cannot contain newlines, carriage returns, tabs, or spaces");
            return 1;
        }
    }

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

        // Unmapped reads
        if (record->core.flag & BAM_FUNMAP) {
            Basic_Seq_Statistics& basic_qc = output_data.unmapped_long_read_info;
            Basic_Seq_Quality_Statistics& seq_quality_info = output_data.unmapped_seq_quality_info;
            this->updateReadAndBaseCounts(record, basic_qc, seq_quality_info, false);

        } else {
            // Calculate base alignment statistics on non-secondary alignments
            Basic_Seq_Statistics& basic_qc = output_data.mapped_long_read_info;
            Basic_Seq_Quality_Statistics& seq_quality_info = output_data.seq_quality_info;
            if (!(record->core.flag & BAM_FSECONDARY)) {

                // Determine if this is a forward or reverse read
                if (record->core.flag & BAM_FREVERSE) {
                    output_data.reverse_alignment++;
                } else {
                    output_data.forward_alignment++;
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

            // Secondary alignment (not included in QC, only read count)
            if (record->core.flag & BAM_FSECONDARY) {
                output_data.num_secondary_alignment++;

                // Get the alignment's query name (the read name)
                std::string query_name = bam_get_qname(record);

                // Update the read's secondary alignments (count once per read)
                output_data.reads_with_secondary[query_name] = true;

            // Supplementary alignment (not included in QC, only read count)
            } else if (record->core.flag & BAM_FSUPPLEMENTARY) {
                output_data.num_supplementary_alignment++;

                // Get the alignment's query name (the read name)
                std::string query_name = bam_get_qname(record);

                // Update the read's supplementary alignments (count once per read)
                output_data.reads_with_supplementary[query_name] = true;

            // Primary alignment
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
                this->updateReadAndBaseCounts(record, basic_qc, seq_quality_info, true);

            } else {
                printError("Error: Unknown alignment type with flag " + std::to_string(record->core.flag));
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
int HTSReader::getNumRecords(const std::string& bam_filename, int thread_count) {
    samFile* bam_file = sam_open(bam_filename.c_str(), "r");
    hts_set_threads(bam_file, thread_count);  // Enable multi-threading
    bam_hdr_t* bam_header = sam_hdr_read(bam_file);
    bam1_t* bam_record = bam_init1();
    int num_reads = 0;
    while (sam_read1(bam_file, bam_header, bam_record) >= 0) {
        num_reads++;
    }

    // Close the BAM file
    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_file);

    return num_reads;
}

void HTSReader::runBaseModificationAnalysis(const std::string &bam_filename, Output_BAM &final_output, double base_mod_threshold, int read_count, int sample_count, int thread_count)
{
    samFile* bam_file = sam_open(bam_filename.c_str(), "r");
    hts_set_threads(bam_file, thread_count);  // Enable multi-threading
    bam_hdr_t* bam_header = sam_hdr_read(bam_file);
    bam1_t* bam_record = bam_init1();
    int64_t read_index = 0;

    // Create a random number generator and seed it with the current time
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);

    // Create a list of read indices to sample, and only keep the first
    // sample_count reads
    std::vector<int> read_indices;
    for (int i = 0; i < read_count; i++) {
        read_indices.push_back(i);
    }
    std::shuffle(read_indices.begin(), read_indices.end(), generator);
    read_indices.resize(sample_count);
    std::unordered_set<int> read_indices_set(read_indices.begin(), read_indices.end());
    printMessage("Number of sampled reads for base modification analysis = " + std::to_string(read_indices_set.size()));

    while (sam_read1(bam_file, bam_header, bam_record) >= 0) {

        // if (read_indices_set.find(read_index) == read_indices_set.end()) {
        //     read_index++;
        //     continue;
        // }

        // Base modification tag analysis
        // Follow here to get base modification tags:
        // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/sam_mods.c
        // https://github.com/samtools/htslib/blob/11205a9ba5e4fc39cc8bb9844d73db2a63fb8119/htslib/sam.h#L2274
        int read_length = bam_record->core.l_qseq;
        hts_base_mod_state *state = hts_base_mod_state_alloc();
        std::vector<std::pair<int32_t, int>> c_modified_positions;  // C-modified positions for CpG analysis (chr->(position, strand))
        std::unordered_map<char, std::unordered_map<char, int>> base_mod_counts;  // Type-specific base modification probabilities (canonical base -> modified base -> [read length %, probability])
        std::unordered_map<char, int> base_primary_count;  // Total base counts for the alignment

        // Parse the base modification tags if a primary alignment
        int read_mod_count = 0;
        int ret = bam_parse_basemod(bam_record, state);
        bool is_primary = !(bam_record->core.flag & BAM_FSECONDARY) && !(bam_record->core.flag & BAM_FSUPPLEMENTARY) && !(bam_record->core.flag & BAM_FUNMAP);

        if (ret >= 0 && is_primary) {
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

            // Get the number of each type of base for the read
            uint8_t *seq = bam_get_seq(bam_record);
            for (int i = 0; i < read_length; i++) {
                char base = seq_nt16_str[bam_seqi(seq, i)];
                base_primary_count[std::toupper(base)]++;
            }

            // Iterate over the state object to get the base modification tags
            // using bam_next_basemod
            hts_base_mod mods[10];
            int n = 0;
            int32_t pos = 0;
            std::vector<int> query_pos;
            bool first_mod_found = false;
            while ((n=bam_next_basemod(bam_record, state, mods, 10, &pos)) > 0) {

                for (int i = 0; i < n; i++) {
                    // Update the modified prediction counts
                    read_mod_count++;  // Read-specific count
                    final_output.modified_prediction_count++;  // Cumulative count
                    char canonical_base_char = std::toupper(mods[i].canonical_base);
                    char mod_type = mods[i].modified_base;
                    // base_mod_counts[mod_type]++;  // Update the type-specific count

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

                        // Update the read length % and probability for the
                        // modification
                        double read_len_pct = (double) (pos + 1) / read_length;
                        // std::cout << "Read length %: " << read_len_pct << ",
                        // probability: " << probability << std::endl;
                        
                        // Update the base modification probabilities for
                        // sampled reads only (10,000 maximum)
                        if (read_indices_set.find(read_index) != read_indices_set.end()) {
                            final_output.updateBaseModProbabilities(mod_type, read_len_pct, probability);  // Update the base modification probabilities
                        }

                        // Update counts for predictions exceeding the threshold
                        if (probability >= base_mod_threshold) {
                            final_output.updateBaseModCounts(mod_type, strand);  // Update the base modification counts

                            // Store the modified positions for later CpG
                            // analysis if a C modification on a primary alignment
                            if (canonical_base_char == 'C' && mod_type != 'C') {

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

            // Append the modified positions to the output data
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
        hts_base_mod_state_free(state);  // Deallocate the base modification state object

        read_index++;  // Update the read index
    }

    bam_destroy1(bam_record);
    bam_hdr_destroy(bam_header);
    sam_close(bam_file);
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
