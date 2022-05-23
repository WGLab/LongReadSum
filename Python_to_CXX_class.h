/*
Python_to_CXX_class.h:
Define the Python bindings from our C++ modules
TODO: Not C++-specific functionality so this wrapper may not be needed by the python module
*/

#ifndef PYTHON_TO_CXX_CLASS_H
#define PYTHON_TO_CXX_CLASS_H
// Python_to_CXX_class

#include <string>
//#include <vector>

//#include <stdint.h>
//#include <cstddef>

#define MAX_INPUT_FILES 2048
//#define PATH_LENGTH 1024


/*
Input_Para:
Python-invoked C++ class
*/
class Input_Para{
private:
   /*char output_folder[PATH_LENGTH];
   char out_prefix[PATH_LENGTH];
   */
public: 
   /*char* get_output_folder(){ return output_folder; }
   char* get_out_prefix(){ return out_prefix; }
   void set_output_folder(const char* p_str);
   void set_out_prefix(const char* p_str);
   */

   int threads;
   
   std::string output_folder;
   //std::vector<std::string> input_files;
   //std::string *input_files;
   std::string input_files[MAX_INPUT_FILES];
   int32_t user_defined_fastq_base_qual_offset = -1; //
   /*//char output_folder[PATH_LENGTH];
   char input_files[MAX_INPUT_FILES][PATH_LENGTH];
   //char out_prefix[PATH_LENGTH];
   int add_input_file(const char* _ip_file);   
   */

   size_t num_input_files;

   std::string out_prefix;

   int64_t other_flags;

   /// functions list;
   std::string add_input_file(const std::string& _ip_file);

   int rdm_seed;
   float downsample_percentage;

   Input_Para();
   // Input_Para(const Input_Para& ip1);
   ~Input_Para();
};

#endif
