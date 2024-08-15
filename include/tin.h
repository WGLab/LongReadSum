#ifndef TIN_H
#define TIN_H

/// @cond
#include <string>
#include <vector>
#include <unordered_map>

#include <htslib/sam.h>
/// @endcond

std::vector<double> calculateTIN(const std::string& gene_bed, const std::string& bam_filepath);

double calculateTINScore(const std::vector<int>& read_depth_values, int junctions);

std::unordered_map<std::string, int> getExonExonJunctions(const std::string& gene_bed);

std::unordered_map<std::string, std::tuple<std::string, int, int>> getGenePositions(const std::string& gene_bed);

std::unordered_map<std::string, std::vector<std::tuple<std::string, int, int>>> getExonPositions(const std::string &gene_bed);

std::pair<std::unordered_map<int, int>, int> getReadDepths(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end);

#endif
