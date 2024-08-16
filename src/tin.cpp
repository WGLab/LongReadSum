#include "tin.h"

/// @cond
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
/// @endcond

std::unordered_map<int, int> getReadDepths(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end)
{
    // Set up the region to fetch reads (1-based)
    std::string region = chr + ":" + std::to_string(start) + "-" + std::to_string(end);
    hts_itr_t* iter = sam_itr_querys(idx, header, region.c_str());
    // hts_itr_t* iter = sam_itr_querys(index, header, region.c_str());
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
        if (record->core.flag & BAM_FUNMAP || record->core.flag & BAM_FSECONDARY || record->core.flag & BAM_FQCFAIL) {
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
            if (op == BAM_CMATCH) {
                for (int j = 0; j < len; j++) {

                    // Check if the position has already been counted
                    if (std::find(read_positions.begin(), read_positions.end(), pos + j) != read_positions.end()) {
                        std::cout << "Read position " << pos + j << " already counted" << std::endl;
                        continue;
                    }

                    // Check if the base is N
                    // if (bam_seqi(seq, query_pos + j) == 15) {
                    if (bam_seqi(bam_get_seq(record), query_pos + j) == 15) {
                        base_skip++;
                        continue;
                    }

                    // Skip if base quality is less than 13
                    if (bam_get_qual(record)[query_pos + j] < 13) {
                        base_skip++;
                        continue;
                    }

                    // Check if the position is within the exon, or if
                    // it is equal to the transcript start+1 or end
                    // if ((pos + j >= exon_start && pos + j <= exon_end) || pos + j == start + 1 || pos + j == end) {
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
        // std::cout << "Base skip count: " << base_skip << std::endl;
    }
    // std::cout << "Read count: " << read_count << std::endl;
    // std::cout << "Read skip count: " << skip_count << std::endl;

    // Destroy the iterator
    hts_itr_destroy(iter);

    // Destroy the record
    bam_destroy1(record);

    return C;
}


std::unordered_map<std::string, int> getExonExonJunctions(const std::string& gene_bed)
{
    // Open the gene BED file
    std::ifstream gene_bed_file(gene_bed);

    if (!gene_bed_file.is_open()) {
        std::cerr << "Error opening gene BED file" << std::endl;
        exit(1);
    }

    // Loop through the gene BED file and get the exon-exon junctions for each
    // transcript (gene ID -> exon-exon junctions)
    std::unordered_map<std::string, int> exon_exon_junctions;
    std::string line;
    while (std::getline(gene_bed_file, line)) {
        // Parse the gene BED line (chrom, start, end, transcript_id) (1-based)
        std::istringstream iss(line);
        std::string chrom;
        int start, end;
        std::string gene_id;
        iss >> chrom >> start >> end >> gene_id;

        // Get the exon-exon junctions for the transcript. If the gene ID is not
        // in the map, initialize it to 0
        if (exon_exon_junctions.find(gene_id) == exon_exon_junctions.end()) {
            exon_exon_junctions[gene_id] = 0;
        } else {
            exon_exon_junctions[gene_id]++;
            // std::cout << "Found exon-exon junction for gene " << gene_id << std::endl;
        }
    }

    // Close the gene BED file
    gene_bed_file.close();

    return exon_exon_junctions;
}

std::unordered_map<std::string, std::tuple<std::string, int, int>> getGenePositions(const std::string& gene_bed)
{
    // Open the gene BED file
    std::ifstream gene_bed_file(gene_bed);

    if (!gene_bed_file.is_open()) {
        std::cerr << "Error opening gene BED file" << std::endl;
        exit(1);
    }

    // Loop through the BED file with exon positions and get the gene start and
    // end positions for each transcript (gene ID -> (chrom, tx_start, tx_end))
    std::unordered_map<std::string, std::tuple<std::string, int, int>> gene_positions;
    std::string line;
    std::string test_gene = "ENSG00000227232.5";
    while (std::getline(gene_bed_file, line)) {
        // Parse the gene BED line (chrom, start, end, transcript_id) (1-based)
        std::istringstream iss(line);
        std::string chrom;
        int start, end;
        std::string gene_id;
        iss >> chrom >> start >> end >> gene_id;

        // Check if the gene ID is already in the map
        if (gene_positions.find(gene_id) == gene_positions.end()) {
            gene_positions[gene_id] = std::make_tuple(chrom, start, end);

        // Else, update the gene start and end positions if the current start
        // position is less than the stored start position, or if the current
        // end position is greater than the stored end position
        } else {
            std::tuple<std::string, int, int> current_positions = gene_positions[gene_id];
            std::string current_chrom = std::get<0>(current_positions);
            int current_start = std::get<1>(current_positions);
            int current_end = std::get<2>(current_positions);

            if (gene_id == test_gene) {
                std::cout << "Found gene " << test_gene << std::endl;
                std::cout << "Positions: " << chrom << " " << start << " " << end << std::endl;
            }

            if (start < current_start) {
                gene_positions[gene_id] = std::make_tuple(chrom, start, current_end);
            }

            if (end > current_end) {
                gene_positions[gene_id] = std::make_tuple(chrom, current_start, end);
            }
        }
    }

    // Print the test gene's range
    std::tuple<std::string, int, int> test_gene_positions = gene_positions[test_gene];
    std::string test_chrom = std::get<0>(test_gene_positions);
    int test_start = std::get<1>(test_gene_positions);
    int test_end = std::get<2>(test_gene_positions);
    std::cout << "Test gene positions: " << test_chrom << " " << test_start << " " << test_end << std::endl;

    // Close the gene BED file
    gene_bed_file.close();

    return gene_positions;
}

std::unordered_map<std::string, std::vector<std::tuple<std::string, int, int>>> getExonPositions(const std::string &gene_bed)
{
    // Open the gene BED file
    std::ifstream gene_bed_file(gene_bed);

    if (!gene_bed_file.is_open()) {
        std::cerr << "Error opening gene BED file" << std::endl;
        exit(1);
    }

    // Loop through the BED file with exon positions and get the exon positions
    // for each transcript (gene ID -> (chrom, start, end))
    std::unordered_map<std::string, std::vector<std::tuple<std::string, int, int>>> exon_positions;
    std::string line;
    while (std::getline(gene_bed_file, line)) {
        // Parse the gene BED line (chrom, start, end, transcript_id) (1-based)
        std::istringstream iss(line);
        std::string chrom;
        int start, end;
        std::string gene_id;
        iss >> chrom >> start >> end >> gene_id;

        // Get the exon positions for the transcript
        if (exon_positions.find(gene_id) == exon_positions.end()) {
            std::vector<std::tuple<std::string, int, int>> positions;
            positions.push_back(std::make_tuple(chrom, start, end));
            exon_positions[gene_id] = positions;
        } else {
            std::vector<std::tuple<std::string, int, int>> positions = exon_positions[gene_id];
            positions.push_back(std::make_tuple(chrom, start, end));
            exon_positions[gene_id] = positions;
        }
    }

    // Close the gene BED file
    gene_bed_file.close();

    return exon_positions;
}

std::vector<double> calculateTIN(const std::string& gene_bed, const std::string& bam_filepath)
{
    // Open the BAM file
    htsFile* bam_file = sam_open(bam_filepath.c_str(), "r");
    if (bam_file == NULL) {
        std::cerr << "Error opening BAM file" << std::endl;
        exit(1);
    }

    // Read the BAM header
    bam_hdr_t* header = sam_hdr_read(bam_file);
    if (header == NULL) {
        std::cerr << "Error reading BAM header" << std::endl;
        exit(1);
    }

    // Get the index for the BAM file
    hts_idx_t* index = sam_index_load(bam_file, bam_filepath.c_str());

    // Read the gene BED file
    std::ifstream gene_bed_file(gene_bed);

    if (!gene_bed_file.is_open()) {
        std::cerr << "Error opening gene BED file" << std::endl;
        exit(1);
    }

    // Vector to store the TIN scores for each transcript (gene ID -> (chrom,
    // tx_start, tx_end, TIN)
    // std::unordered_map<std::string, std::tuple<std::string, int, int, double>>
    //     tin_scores;

    // Loop through the gene BED file and calculate the TIN score for each
    // transcript
    std::vector<double> TIN_scores;
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
        

        // Remove the last comma from the exon sizes and starts strings
        if (exon_sizes_str.back() == ',') {
            exon_sizes_str.pop_back();
        }
        if (exon_starts_str.back() == ',') {
            exon_starts_str.pop_back();
        }

        // Get the comma-separated exon sizes for the transcript
        std::vector<int> exon_sizes;
        // std::istringstream exon_sizes_iss(exon_sizes_str);
        // std::istringstream exon_starts_iss(exon_starts_str);
        // std::cout << "Exon sizes: " << exon_sizes_str << std::endl;
        while (exon_sizes_str.find(",") != std::string::npos) {
            int pos = exon_sizes_str.find(",");
            std::string exon_size_str = exon_sizes_str.substr(0, pos);            
            exon_sizes.push_back(std::stoi(exon_size_str));
            exon_sizes_str.erase(0, pos + 1);
            // std::cout << "Exon size: " << exon_size_str << std::endl;
        }
        exon_sizes.push_back(std::stoi(exon_sizes_str));
        // std::cout << "Exon size: " << exon_sizes_str << std::endl;

        // Get the comma-separated exon starts for the transcript
        std::vector<int> exon_starts;
        // std::cout << "Exon starts: " << exon_starts_str << std::endl;
        while (exon_starts_str.find(",") != std::string::npos) {
            int pos = exon_starts_str.find(",");
            std::string exon_start_str = exon_starts_str.substr(0, pos);
            exon_starts.push_back(std::stoi(exon_start_str));
            exon_starts_str.erase(0, pos + 1);
            // std::cout << "Exon start: " << exon_start_str << std::endl;
        }
        exon_starts.push_back(std::stoi(exon_starts_str));
        // std::cout << "Exon start: " << exon_starts_str << std::endl;

        // int exon_size;
        // while (exon_sizes_iss >> exon_size) {
        //     std::cout << "Exon size: " << exon_size << std::endl;
        //     exon_sizes.push_back(exon_size);
        // }

        // int exon_start;
        // while (exon_starts_iss >> exon_start) {
        //     // Add 1 to the exon start to make it 1-based
        //     exon_start++;
        //     std::cout << "Exon start: " << exon_start << std::endl;
        //     exon_starts.push_back(exon_start);
        // }

        std::cout << "Exon count: " << exon_count << std::endl;
        std::cout << "Transcript: " << name << std::endl;
        std::cout << "Chrom: " << chrom << std::endl;
        std::cout << "Start: " << start << std::endl;
        std::cout << "End: " << end << std::endl;

        // Get the read depths and cumulative read depth for each exon
        std::unordered_map<int, int> C;
        for (int i = 0; i < exon_count; i++) {
            // int exon_start = start + exon_starts[i];
            int exon_start = start + 1 + exon_starts[i];
            int exon_end = exon_start + exon_sizes[i] - 1;
            // int exon_end = exon_start + exon_sizes[i];

            // std::cout << "Exon start: " << exon_start << ", Exon end: " << exon_end << std::endl;
            // std::string region = chrom + ":" + std::to_string(exon_start) + "-"
            //     + std::to_string(exon_end);
            // std::string region = chrom + ":" + std::to_string(start+1) + "-"
            //     + std::to_string(end);
            // std::cout << "Region (positions range): " << region << std::endl;

            // Get the depths and cumulative depths for the region
            std::unordered_map<int, int> exon_depths = getReadDepths(bam_file, index, header, chrom, exon_start, exon_end);
            for (const auto& depth : exon_depths) {
                C[depth.first] = depth.second;
            }
        }

        // Get the read depths for the transcript start+1 and end
        std::vector<int> transcript_positions = {start + 1, end};
        for (const auto& position : transcript_positions) {
            std::unordered_map<int, int> transcript_depths = getReadDepths(bam_file, index, header, chrom, position, position);
            for (const auto& depth : transcript_depths) {
                C[depth.first] = depth.second;
                // C[depth.first] = depth.second;
            }
        }

        // Determine the sample size for the transcript (transcript start,
        // end, + exon lengths)
        int sample_size = 2;
        for (const auto& exon_size : exon_sizes) {
            sample_size += exon_size;
        }
        // std::cout << "Sample size: " << sample_size << std::endl;

        // Sort C by position
        std::vector<int> positions;
        for (const auto& entry : C) {
            positions.push_back(entry.first);
        }
        std::sort(positions.begin(), positions.end());

        // Print the positions and read depth values
        // std::cout << "Read depth values: " << std::endl;
        // for (const auto& position : positions) {
        //     std::string Ci_str = std::to_string(C[position]) + ".0";
        //     std::cout << "Position: " << position << " Coverage: " << Ci_str << std::endl;
        //     // std::cout << "C[" << position << "]: " << C[position] << std::endl;
        // }
        std::cout << std::endl;
        std::cout << "Length of C: " << C.size() << std::endl;

        // Calculate total coverage, ΣCi
        int sigma_Ci = 0;
        for (const auto& entry : C) {
            sigma_Ci += entry.second;
        }
        std::cout << "Sigma Ci: " << sigma_Ci << std::endl;

        // Calculate the relative coverage, Pi = Ci / ΣCi
        std::vector<double> Pi;
        if (sigma_Ci > 0) {
            for (const auto& entry : C) {
                double Pi_value = static_cast<double>(entry.second) / sigma_Ci;
                Pi.push_back(Pi_value);
                // std::cout << "Pi[" << entry.first << "]: " << Pi_value << std::endl;
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
            }  // else { H += 0; }
        }
        std::cout << "H: " << -H << std::endl;

        double TIN = 0;
        if (H != 0) {

            // Calculate U
            double U = std::exp(-H);

            // Calculate the TIN score for the exon
            // int k = exon_end - exon_start + 1;
            int k = sample_size;
            std::cout << "sample size: " << k << std::endl;
            std::cout << "TIN calculation: " << U << " / " << k << " = " << std::to_string(U / k) << std::endl;
            TIN = 100.0 * (U / k);
        }

        // Print the TIN score
        std::cout << "TIN score: " << TIN << std::endl;
        TIN_scores.push_back(TIN);
    }

    // Close the BAM file
    sam_close(bam_file);

    // Print the TIN mean, median, and standard deviation
    double TIN_sum = 0;
    double TIN_sum_sq = 0;
    for (const auto& TIN_score : TIN_scores) {
        TIN_sum += TIN_score;
        TIN_sum_sq += TIN_score * TIN_score;
    }
    double TIN_mean = TIN_sum / TIN_scores.size();
    double TIN_stddev = std::sqrt((TIN_sum_sq - TIN_sum * TIN_sum / TIN_scores.size()) / (TIN_scores.size() - 1));
    std::cout << "TIN mean: " << TIN_mean << std::endl;
    std::cout << "TIN median: " << TIN_scores[TIN_scores.size() / 2] << std::endl;
    std::cout << "TIN standard deviation: " << TIN_stddev << std::endl;

    // Close the gene BED file
    gene_bed_file.close();

    // Destroy the header
    bam_hdr_destroy(header);

    // Destroy the index
    hts_idx_destroy(index);

    // Return the TIN scores
    return std::vector<double>();
}
