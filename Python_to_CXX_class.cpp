#include "Python_to_CXX_class.h"


Input_Para::Input_Para(){
   input_files = new std::string[MAX_INPUT_FILES];
   num_input_files = 0;
   rdm_seed = 1;
   downsample_percentage = 100;
}

Input_Para::~Input_Para(){
   delete [] input_files;
}

std::string Input_Para::add_input_file(std::string _ip_file){
   if ( num_input_files < MAX_INPUT_FILES){
      input_files[ num_input_files ] = _ip_file;
      num_input_files++;
      return "";
   }else{
       return "Only "+std::to_string(MAX_INPUT_FILES)+" input files are supported!!";
   }
}

Input_Para::Input_Para(const Input_Para& ip1){
   threads = ip1.threads;
   output_folder = ip1.output_folder;
   
   out_prefix = ip1.out_prefix;
   other_flags = ip1.other_flags;

   rdm_seed = ip1.rdm_seed;
   downsample_percentage = ip1.downsample_percentage;

   num_input_files = ip1.num_input_files;
   input_files = new std::string[MAX_INPUT_FILES];
   for (size_t nip=0; nip<num_input_files; nip++){
      input_files[ nip ] = ip1.input_files[ nip ];
   }
}
