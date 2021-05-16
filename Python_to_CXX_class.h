#ifndef PYTHON_TO_CXX_CLASS_H
#define PYTHON_TO_CXX_CLASS_H
// Python_to_CXX_class

#include <string>
//#include <vector>

#define MAX_INPUT_FILES 2048

class Input_Para{
public:
   int threads;
   
   std::string output_folder;
   //std::vector<std::string> input_files;
   //std::string *input_files;
   std::string input_files[MAX_INPUT_FILES];
   size_t num_input_files;

   std::string out_prefix;

   int64_t other_flags;  

   /// functions list;
   std::string add_input_file(std::string _ip_file);

   int rdm_seed;
   float downsample_percentage;

   Input_Para();
   Input_Para(const Input_Para& ip1);
   ~Input_Para();
};

#endif
