#include <string>

#include "input_parameters.h"

Input_Para::Input_Para(){
    // Set default parameters
    this->num_input_files = 0;
    this->threads = 1;
    this->rdm_seed = 1;
    this->downsample_percentage = 100;
    this->other_flags = 0;
    this->user_defined_fastq_base_qual_offset = -1;
    this->rrms_csv = "";
    this->ref_genome = "";
    this->base_mod_threshold = 0.5;
    this->gene_bed = "";
    this->mod_analysis = false;
}

Input_Para::~Input_Para(){
}


std::string Input_Para::add_input_file(const std::string& input_filepath){
   if (this->num_input_files < MAX_INPUT_FILES){
      this->input_files[ this->num_input_files ] = input_filepath;
      this->num_input_files++;
      return "";
   }else{
       return "Only "+std::to_string(MAX_INPUT_FILES)+" input files are supported!!";
   }
}
