#include "tin.h"

/// @cond
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
/// @endcond

#include "tin_stats.h"

std::unordered_map<int, int> getReadDepths(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end, bool transcript_strand)
{
    // Set up the region to fetch reads (1-based)
    std::string region = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str());
    if (iter == NULL) {
        std::cerr << "Error creating iterator for region " << region << std::endl;
        exit(1);
    }

    // Initialize a map to store the read depth values at each position
    std::unordered_map<int, int> C;
    for (int i = start; i <= end; i++) {
        C[i] = 0;
    }

    // Loop through the reads in the region and calculate the read depth
    // values
    bam1_t* record = bam_init1();
    int read_count = 0;
    int skip_count = 0;

    // Vector for keeping track of each read's positions and avoid
    // counting overlapping reads
    std::vector<int> read_positions;
    while (sam_itr_next(bam_file, iter, record) >= 0) {

        // Skip unmapped and secondary alignments, and QC failures
        if (record->core.flag & BAM_FUNMAP || record->core.flag & BAM_FSECONDARY || record->core.flag & BAM_FQCFAIL || record->core.flag & BAM_FDUP) {
            skip_count++;
            continue;
        }

        // Skip if on a different strand
        bool read_strand = (record->core.flag & BAM_FREVERSE) ? false : true;
        if (read_strand != transcript_strand) {
            skip_count++;
            continue;
        }
        
        read_count++;

        // Clear positions for each read
        read_positions.clear();

        // Get the alignment position (0-based)
        int pos = record->core.pos + 1;  // 1-based
        int query_pos = 0;

        // Loop through the CIGAR string and update depth at matches
        uint32_t* cigar = bam_get_cigar(record);
        int base_skip = 0;
        for (uint32_t i = 0; i < record->core.n_cigar; i++) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                for (int j = 0; j < len; j++) {

                    // Check if the position has already been counted
                    if (std::find(read_positions.begin(), read_positions.end(), pos + j) != read_positions.end()) {
                        // std::cout << "Read position " << pos + j << " already counted" << std::endl;
                        continue;
                    }

                    // Check if the base is N
                    // if (bam_seqi(seq, query_pos + j) == 15) {
                    if (bam_seqi(bam_get_seq(record), query_pos + j) == 15) {
                        base_skip++;
                        // std::cout << "Base is N at position " << pos + j << std::endl;
                        continue;
                    }

                    // Skip if base quality is less than 13
                    if (bam_get_qual(record)[query_pos + j] < 13) {
                        base_skip++;
                        // std::cout << "Base quality is less than 13 at position " << pos + j << std::endl;
                        continue;
                    }

                    // Check if the position is within the exon
                    if (pos + j >= start && pos + j <= end) {
                        C[pos + j]++;

                        // Add the position to the read positions
                        read_positions.push_back(pos + j);
                    }
                }
            }

            // Update the reference position
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
                pos += len;
            }

            // Update the query position
            if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF) {
                query_pos += len;
            }
        }
    }

    // Destroy the iterator
    hts_itr_destroy(iter);

    // Destroy the record
    bam_destroy1(record);

    return C;
}

bool checkMinReads(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end, int min_reads)
{
    // Set up the region to fetch reads (1-based)
    std::string region = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str());
    if (iter == NULL) {
        std::cerr << "Error creating iterator for region " << region << std::endl;
        exit(1);
    }

    // Create a set for storing the read positions
    std::unordered_set<int> read_positions;

    // Calculate the read depth values for the region
    int read_count = 0;
    bool min_reads_met = false;
    bam1_t* record = bam_init1();
    while (sam_itr_next(bam_file, iter, record) >= 0) {
        // Skip unmapped and secondary alignments, and QC failures
        if (record->core.flag & BAM_FUNMAP || record->core.flag & BAM_FSECONDARY || record->core.flag & BAM_FQCFAIL) {
            continue;
        }

        // Get the alignment position
        int pos = record->core.pos;
        if (pos < start) {
            // std::cout << "Read position " << pos << " is less than transcript start " << start << std::endl;
            continue;
        }

        if (pos >= end) {
            // std::cout << "Read position " << pos << " is greater than transcript end " << end << std::endl;
            continue;
        }

        // Add the position to the set
        read_positions.insert(pos);
        read_count++;

        // Break if the minimum read count is reached
        if ((int) read_positions.size() > min_reads) {
            min_reads_met = true;
            break;
        }
    }

    // Destroy the iterator
    hts_itr_destroy(iter);

    // Destroy the record
    bam_destroy1(record);

    return min_reads_met;
}

void calculateTIN(TINStats* tin_stats, const std::string& gene_bed, const std::string& bam_filepath, int min_cov, int sample_size, const std::string& output_folder, int thread_count)
{
    std::cout << "Using TIN minimum coverage " << min_cov << " and sample size " << sample_size << std::endl;

    // Open the BAM file
    htsFile* bam_file = sam_open(bam_filepath.c_str(), "r");
    if (bam_file == NULL) {
        std::cerr << "Error opening BAM file" << std::endl;
        exit(1);
    }

    // Enable multi-threading
    hts_set_threads(bam_file, thread_count);

    // Read the BAM header
    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        std::cerr << "Error reading BAM header" << std::endl;
        exit(1);
    }

    // Get the index for the BAM file
    hts_idx_t* index = sam_index_load(bam_file, bam_filepath.c_str());
    if (index == NULL) {
        std::cerr << "Error loading BAM index" << std::endl;
        exit(1);
    }

    // Read the gene BED file
    std::ifstream gene_bed_file(gene_bed);
    if (!gene_bed_file.is_open()) {
        std::cerr << "Error opening gene BED file" << std::endl;
        exit(1);
    }

    // Vector to store the TIN scores for each transcript (gene ID -> (chrom,
    // tx_start, tx_end, TIN)
    TINMap tin_map;

    // Loop through the gene BED file and calculate the TIN score for each
    // transcript
    std::cout << "Calculating TIN scores for each transcript..." << std::endl;
    std::vector<double> TIN_scores;
    std::vector<std::string> gene_ids;
    std::string line;
    while (std::getline(gene_bed_file, line)) {
        // Parse the gene BED line and get the exon positions for the transcript
        std::istringstream iss(line);
        std::string chrom;
        int start, end;
        std::string name;
        std::string score;
        std::string strand;
        int thick_start, thick_end;
        std::string item_rgb;
        int exon_count;
        std::string exon_sizes_str;
        std::string exon_starts_str;

        iss >> chrom >> start >> end >> name >> score >> strand >> thick_start
            >> thick_end >> item_rgb >> exon_count >> exon_sizes_str
            >> exon_starts_str;

        // Check if the transcript passes the minimum read depth threshold
        if (!checkMinReads(bam_file, index, header, chrom, start, end, min_cov)) {
            // std::cout << "Skipping transcript " << name << " because it does not meet the minimum coverage requirement" << std::endl;
            
            // Set the TIN score to 0
            tin_map[name] = std::make_tuple(chrom, start, end, 0.0);
            gene_ids.push_back(name);

            continue;
        }

        std::cout << "Calculating TIN for transcript " << name << " with strand " << strand << std::endl;
        bool strand_bool = true;
        if (strand == "-") {
            strand_bool = false;
        }

        // Remove the last comma from the exon sizes and starts strings
        if (exon_sizes_str.back() == ',') {
            exon_sizes_str.pop_back();
        }
        if (exon_starts_str.back() == ',') {
            exon_starts_str.pop_back();
        }

        // Get the comma-separated exon sizes for the transcript
        std::vector<int> exon_sizes;
        while (exon_sizes_str.find(",") != std::string::npos) {
            int pos = exon_sizes_str.find(",");
            std::string exon_size_str = exon_sizes_str.substr(0, pos);            
            exon_sizes.push_back(std::stoi(exon_size_str));
            exon_sizes_str.erase(0, pos + 1);
        }
        exon_sizes.push_back(std::stoi(exon_sizes_str));

        // Get the comma-separated exon starts for the transcript
        std::vector<int> exon_starts;
        while (exon_starts_str.find(",") != std::string::npos) {
            int pos = exon_starts_str.find(",");
            std::string exon_start_str = exon_starts_str.substr(0, pos);
            exon_starts.push_back(std::stoi(exon_start_str));
            exon_starts_str.erase(0, pos + 1);
            // std::cout << "Exon start: " << exon_start_str << std::endl;
        }
        exon_starts.push_back(std::stoi(exon_starts_str));

        // Get the read depths and cumulative read depth for each exon
        std::unordered_map<int, int> C;
        for (int i = 0; i < exon_count; i++) {
            int exon_start = start + 1 + exon_starts[i];
            int exon_end = exon_start + exon_sizes[i] - 1;

            // Get the depths and cumulative depths for the region
            std::unordered_map<int, int> exon_depths = getReadDepths(bam_file, index, header, chrom, exon_start, exon_end, strand_bool);
            for (const auto& depth : exon_depths) {
                C[depth.first] = depth.second;
            }
        }

        // Get the read depths for the transcript start+1 and end
        std::vector<int> transcript_positions = {start + 1, end};
        for (const auto& position : transcript_positions) {
            std::unordered_map<int, int> transcript_depths = getReadDepths(bam_file, index, header, chrom, position, position, strand_bool);
            for (const auto& depth : transcript_depths) {
                C[depth.first] = depth.second;
            }
        }

        // Determine the sample size for the transcript (transcript start,
        // end, + exon lengths)
        int transcript_size = 2;
        for (const auto& exon_size : exon_sizes) {
            transcript_size += exon_size;
        }
        // std::cout << "mRNA size: " << transcript_size - 2 << " for transcript " << name << std::endl;

        // Sort C by position
        std::vector<int> positions;
        for (const auto& entry : C) {
            positions.push_back(entry.first);
        }
        std::sort(positions.begin(), positions.end());

        // Sample the values if the mRNA size is greater than the user-specified
        // sample size
        int mRNA_size = transcript_size - 2;
        if (mRNA_size > sample_size) {
            // std::cout << "Sampling " << sample_size << " positions" << std::endl;
            std::vector<int> sampled_positions;
            int step = mRNA_size / sample_size;
            for (size_t i = 0; i < positions.size(); i += step) {
                sampled_positions.push_back(positions[i]);
            }
            positions = sampled_positions;

            // Also add all the exon start and end positions
            for (size_t i = 0; i < exon_starts.size(); i++) {
                positions.push_back(start + 1 + exon_starts[i]);
                positions.push_back(start + 1 + exon_starts[i] + exon_sizes[i] - 1);
                // std::cout << "Added exon start: " << start + 1 + exon_starts[i] << ", end: " << start + 1 + exon_starts[i] + exon_sizes[i] - 1 << std::endl;
            }

            // Remove duplicates
            std::sort(positions.begin(), positions.end());
            positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

            // Create a new map with the sampled positions
            std::unordered_map<int, int> sampled_C;
            for (const auto& position : positions) {
                // std::cout << "Position: " << position << " Coverage: " << C[position] << std::endl;
                sampled_C[position] = C[position];
            }
            C = sampled_C;

            // Update the sample size
            transcript_size = positions.size();
        }

        // // Print the positions and read depth values
        // std::cout << "Read depth values: " << std::endl;
        // for (const auto& position : positions) {
        //     std::string Ci_str = std::to_string(C[position]) + ".0";
        //     std::cout << "Position: " << position << " Coverage: " << Ci_str << std::endl;
        //     // std::cout << "C[" << position << "]: " << C[position] << std::endl;
        // }
        // std::cout << "Length of C: " << C.size() << std::endl;

        // Calculate total coverage, ΣCi
        int sigma_Ci = 0;
        for (const auto& entry : C) {
            sigma_Ci += entry.second;
        }
        // std::cout << "Sigma Ci: " << sigma_Ci << std::endl;

        // Calculate the relative coverage, Pi = Ci / ΣCi
        std::vector<double> Pi;
        if (sigma_Ci > 0) {
            for (const auto& entry : C) {
                double Pi_value = static_cast<double>(entry.second) / sigma_Ci;
                Pi.push_back(Pi_value);
            }
        } else {
            Pi = std::vector<double>(C.size(), 0);
        }

        // Calculate H
        double H = 0;
        for (const auto& Pi_value : Pi) {
            // Use the natural logarithm
            if (Pi_value > 0) {
                H += Pi_value * std::log(Pi_value);
            }
        }
        // std::cout << "H: " << -H << std::endl;

        double TIN = 0;
        if (H != 0) {

            // Calculate U
            double U = std::exp(-H);

            // Calculate the TIN score for the exon
            int k = transcript_size;
            // std::cout << "sample size: " << k << std::endl;
            // std::cout << "TIN calculation: " << U << " / " << k << " = " << std::to_string(U / k) << std::endl;
            TIN = 100.0 * (U / k);
        }

        // Print the TIN score
        // std::cout << "TIN for transcript " << name << ": " << TIN << std::endl;
        TIN_scores.push_back(TIN);
        gene_ids.push_back(name);

        // Store the TIN score for the transcript
        tin_map[name] = std::make_tuple(chrom, start, end, TIN);

        // Log every 1000 transcripts
        if (gene_ids.size() % 1000 == 0) {
            std::cout << "Processed " << gene_ids.size() << " transcripts" << std::endl;
        }
    }

    // Close the BAM file
    sam_close(bam_file);

    // Close the gene BED file
    gene_bed_file.close();

    // Destroy the header
    bam_hdr_destroy(header);

    // Destroy the index
    hts_idx_destroy(index);

    if (TIN_scores.size() == 0) {
        std::cerr << "No TIN scores calculated" << std::endl;
    } else {
        std::cout << "Calculating TIN summary for " << TIN_scores.size() << " transcripts..." << std::endl;

        // Print the TIN mean, median, and standard deviation
        double TIN_sum = 0;
        for (const auto& TIN_score : TIN_scores) {
            TIN_sum += TIN_score;
        }
        double TIN_mean = TIN_sum / TIN_scores.size();

        // Calculate the standard deviation
        double TIN_sum_sq = 0;
        for (const auto& TIN_score : TIN_scores) {
            TIN_sum_sq += (TIN_score - TIN_mean) * (TIN_score - TIN_mean);
        }

        double TIN_variance = TIN_sum_sq / TIN_scores.size();
        double TIN_stddev = std::sqrt(TIN_variance);

        // std::cout << "Number of TIN scores: " << TIN_scores.size() << std::endl;

        // Set the precision for the output
        std::cout << std::fixed << std::setprecision(14);

        // std::cout << "TIN mean: " << TIN_mean << std::endl;

        // Sort the TIN scores
        std::sort(TIN_scores.begin(), TIN_scores.end());

        // Calculate the median
        double TIN_median = TIN_scores[TIN_scores.size() / 2];
        // std::cout << "TIN median: " << TIN_median << std::endl;
        // std::cout << "TIN standard deviation: " << TIN_stddev << std::endl;

        std::cout << "Writing TIN scores to file..." << std::endl;

        // Write the TIN scores to a file
        std::string output_tin_tsv = output_folder + "/tin_scores.tsv";
        std::ofstream output_tin_file(output_tin_tsv);
        output_tin_file << std::fixed << std::setprecision(14);

        if (!output_tin_file.is_open()) {
            std::cerr << "Error opening output TIN file" << std::endl;
            exit(1);
        }

        // Write the header
        output_tin_file << "geneID\tchrom\ttx_start\ttx_end\tTIN" << std::endl;

        // Write the TIN scores to the file in the order of the gene IDs
        for (const auto& gene_id : gene_ids) {
            std::string chrom = std::get<0>(tin_map[gene_id]);
            int tx_start = std::get<1>(tin_map[gene_id]);
            int tx_end = std::get<2>(tin_map[gene_id]);
            double TIN = std::get<3>(tin_map[gene_id]);

            output_tin_file << gene_id << "\t" << chrom << "\t" << tx_start << "\t" << tx_end << "\t" << TIN << std::endl;
        }

        // Close the output TIN file
        output_tin_file.close();

        std::cout << "TIN scores written to " << output_tin_tsv << std::endl;

        // Write the TIN summary to a file
        std::string output_tin_summary_tsv = output_folder + "/tin_summary.tsv";

        std::ofstream output_tin_summary_file(output_tin_summary_tsv);
        output_tin_summary_file << std::fixed << std::setprecision(14);

        if (!output_tin_summary_file.is_open()) {
            std::cerr << "Error opening output TIN summary file" << std::endl;
            exit(1);
        }

        // Write the header
        output_tin_summary_file << "Bam_file\tTIN(mean)\tTIN(median)\tTIN(stddev)" << std::endl;

        // Write the TIN summary to the file
        output_tin_summary_file << bam_filepath << "\t" << TIN_mean << "\t" << TIN_median << "\t" << TIN_stddev << std::endl;

        // Close the output TIN summary file
        output_tin_summary_file.close();

        std::cout << "TIN summary written to " << output_tin_summary_tsv << std::endl;

        // Update the TIN stats struct
        tin_stats->mean = TIN_mean;
        tin_stats->median = TIN_median;
        tin_stats->stddev = TIN_stddev;
        tin_stats->num_transcripts = TIN_scores.size();
    }
}
