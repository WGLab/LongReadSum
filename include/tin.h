#ifndef TIN_H
#define TIN_H

/// @cond
#include <string>
#include <vector>
/// @endcond

std::vector<double> calculateTIN(const std::string& gene_bed, const std::string& bam_filepath);

double calculateTINScore(const std::vector<int>& read_depth_values);

#endif
