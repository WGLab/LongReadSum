#ifndef TIN_H
#define TIN_H

/// @cond
#include <string>
#include <vector>
#include <unordered_map>

#include <htslib/sam.h>
/// @endcond

#include "tin_stats.h"

typedef std::unordered_map<std::string, std::tuple<std::string, int, int, double>> TINMap;

// Calculate the TIN score for each transcript in the gene BED file
// (Reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0922-z#Sec11)
void calculateTIN(TINStats* tin_stats, const std::string& gene_bed, const std::string& bam_filepath, int min_cov, int sample_size, const std::string& output_folder, int thread_count);

std::unordered_map<int, int> getReadDepths(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end);

bool checkMinReads(htsFile* bam_file, hts_idx_t* idx, bam_hdr_t* header, std::string chr, int start, int end, int min_reads);

#endif
