#include "tin.h"

/// @cond
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <htslib/sam.h>
/// @endcond

std::vector<double> calculateTIN(const std::string& gene_bed, const std::string& bam_filepath) {
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
    std::unordered_map<std::string, std::tuple<std::string, int, int, double>>
        tin_scores;
    
    // std::unordered_map<std::string, std::u

    // Loop through the gene BED file and calculate the TIN score for each
    // transcript
    std::string line;
    while (std::getline(gene_bed_file, line)) {
        // Parse the gene BED line (chrom, start, end, transcript_id)
        std::istringstream iss(line);
        std::string chrom;
        int start, end;
        std::string gene_id;
        iss >> chrom >> start >> end >> gene_id;

        // Format the region string for the BAM index
        std::string region = chrom + ":" + std::to_string(start) + "-" + std::to_string(end);

        // Set up the region to fetch reads (1-based)
        hts_itr_t* iter = sam_itr_querys(index, header, region.c_str());
        if (iter == NULL) {
            std::cerr << "Error creating iterator for region " << chrom << ":"
                      << start << "-" << end << std::endl;
            exit(1);
        }


        // Loop through the reads in the region
        std::unordered_map<int, int> read_depth_values;  // Position (1-based) -> read depth
        bam1_t* record = bam_init1();
        while (sam_itr_next(bam_file, iter, record) >= 0) {
            // Get the alignment position (0-based)
            int pos = record->core.pos;
            for (int i = 0; i < record->core.l_qseq; i++) {
                // Increment the read depth value at the position
                if (pos + i >= start && pos + i < end) {

                    // Increment the read depth value at the position (1-based)
                    read_depth_values[pos + i + 1]++;
                }
            }
        }
        std::cout << "[TEST] TIN calculation for transcript " << gene_id << std::endl;

        // Print the transcript location
        std::cout << "chrom\tstart\tend" << std::endl;
        std::cout << chrom << "\t" << start << "\t" << end << std::endl;

        // Print the read depth values
        std::cout << "Position\tRead Depth" << std::endl;
        for (const auto& entry : read_depth_values) {
            std::cout << entry.first << "\t" << entry.second << std::endl;
        }

        // [TEST] Exit after the first iteration
        break;

        // Calculate the TIN score for the transcript
        // double tin_score = calculateTINScore(read_depth_values);
    }

    // Close the BAM file
    sam_close(bam_file);

    // Close the gene BED file
    gene_bed_file.close();

    // Destroy the header
    bam_hdr_destroy(header);

    // Return the TIN scores
    return std::vector<double>();
}

// Function to calculate the transcript integrity number (TIN) score for the
// transcript
// (Reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z#Sec11)
double calculateTINScore(const std::vector<int>& read_depth_values) {
    int num_values = read_depth_values.size();

    // TODO: Calculate U = e(-SUM(Pi * log(Pi) for i = 1 to k+j)), where k is the
    // number of equally spaced positions, and j is the number of exon-exon
    // junctions
    int sum = 0;

    for (int i = 0; i < num_values; i++) {
        sum += read_depth_values[i];
    }

    double average = static_cast<double>(sum) / num_values;
    double tinScore = average / sum;

    return tinScore;
}

// TIN::~TIN()
// {
// }
