#ifndef PYTHON_TO_CXX_CLASS_H
#define PYTHON_TO_CXX_CLASS_H
// Python_to_CXX_class

#include <string>
#include <vector>


class Input_Para{
public:
   int threads;
   
   std::string output_folder;
   std::vector<std::string> input_files;

   std::string out_prefix;

   int64_t other_flags;  
};

#endif
