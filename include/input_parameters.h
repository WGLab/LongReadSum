/*
input_parameters.h:
Define the Python bindings from our C++ modules
*/

#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <string>
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
    int rdm_seed;
    float downsample_percentage;
    int32_t user_defined_fastq_base_qual_offset;
    std::string output_folder;  // Output folder
    std::string input_files[MAX_INPUT_FILES];  // Input files

    // Functions
    std::string add_input_file(const std::string& _ip_file);

    Input_Para();

    ~Input_Para();
};

#endif
