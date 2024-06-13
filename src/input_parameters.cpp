#include <string>

#include "input_parameters.h"

Input_Para::Input_Para(){
    // Set default parameters
    num_input_files = 0;
    threads = 1;
    rdm_seed = 1;
    downsample_percentage = 100;
    other_flags = 0;
    user_defined_fastq_base_qual_offset = -1;
    rrms_csv = "";
    ref_genome = "";
}

Input_Para::~Input_Para(){
}


std::string Input_Para::add_input_file(const std::string& _ip_file){
   if ( num_input_files < MAX_INPUT_FILES){
      input_files[ num_input_files ] = _ip_file;
      num_input_files++;
      return "";
   }else{
       return "Only "+std::to_string(MAX_INPUT_FILES)+" input files are supported!!";
   }
}
