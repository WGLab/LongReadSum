/*
input_parameters.h:
Define the Python bindings from our C++ modules
*/

#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <vector>
#include <string>
#include <unordered_set>  // For RRMS read ID filtering
#define MAX_INPUT_FILES 2048


/*
Input_Para:
Define the input parameters for the program
*/
class Input_Para{
public: 
    // Parameters
    int threads;
    size_t num_input_files;
    std::string out_prefix;
    int64_t other_flags;
    int32_t user_defined_fastq_base_qual_offset;
    std::string output_folder;  // Output folder
    std::string input_files[MAX_INPUT_FILES];  // Input files
    std::string read_ids;  // Read IDs comma-separated (FAST5 signal module)
    std::string rrms_csv;  // CSV file with accepted/rejected read IDs (RRMS module)
    bool rrms_filter;  // Generate RRMS stats for accepted (true) or rejected (false) reads
    std::unordered_set<std::string> rrms_read_ids;  // List of read IDs from RRMS CSV file (accepted or rejected)
    std::string ref_genome;  // Reference genome file for BAM base modification analysis
    double base_mod_threshold;  // Base modification threshold for BAM base modification analysis
    std::string gene_bed;  // Gene BED file for RNA-Seq transcript quantification (TIN)
    bool mod_analysis;  // Perform base modification analysis on BAM file
    int tin_sample_size;  // Number of equally spaced samples for TIN calculation
    int tin_min_coverage;  // Minimum coverage for TIN calculation

    // Functions
    std::string add_input_file(const std::string& input_filepath);

    Input_Para();

    ~Input_Para();
};

#endif
