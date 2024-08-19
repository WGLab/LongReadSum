#ifndef TIN_H
#define TIN_H

/// @cond
#include <string>
#include <vector>
#include <unordered_map>

#include <htslib/sam.h>
/// @endcond

// Calculate the TIN score for each transcript in the gene BED file
// (Reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z#Sec11)
std::vector<double> calculateTIN(const std::string& gene_bed, const std::string& bam_filepath, int min_cov, int sample_size);

std::unordered_map<std::string, int> getExonExonJunctions(const std::string& gene_bed);

std::unordered_map<std::string, std::tuple<std::string, int, int>> getGenePositions(const std::string& gene_bed);

std::unordered_map<std::string, std::vector<std::tuple<std::string, int, int>>> getExonPositions(const std::string &gene_bed);

std::unordered_map<int, int> getReadDepths(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end);

bool checkMinReads(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end, int min_reads);

#endif
