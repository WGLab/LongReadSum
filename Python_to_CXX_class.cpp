#include "Python_to_CXX_class.h"


Input_Para::Input_Para(){
   input_files = new std::string[MAX_INPUT_FILES];
   num_input_files = 0;
}

Input_Para::Input_Para(){
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


