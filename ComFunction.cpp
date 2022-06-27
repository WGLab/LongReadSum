#include "ComFunction.h"

#include "glob.h"

#include <fstream>
#include <math.h>
#include <algorithm>

#include <sys/stat.h>
#include <unistd.h>

#include <iostream>


#include <iostream>

// Returns false if the input is larger than the specified buffer size
bool isExpectedLength(char * destination_buffer, const char * input_string, int max_buffer_size){
   int buffer_count = snprintf(destination_buffer, max_buffer_size, "%s", input_string);
   if (buffer_count<0 || buffer_count>=max_buffer_size){
       fprintf(stderr, "Input (%s) is larger than the expected buffer size (%d)", input_string, max_buffer_size);
       return false;
   }
   return true;
}

// Convert the repeat pattern string into a vector containing its data (Should be 4 elements).
std::vector<std::string> readRepeatDataFromString(const std::string & m_str, std::string multi_delimiters, bool contain_delimiter){
    std::vector<std::string> m_substr_list;
    m_substr_list.reserve(10);

    size_t cur_pos = 0;
    size_t last_pos = 0;
    while ((cur_pos=m_str.find_first_of(multi_delimiters, last_pos))!=std::string::npos){
        m_substr_list.push_back(m_str.substr(last_pos, cur_pos-last_pos+(contain_delimiter?1:0)));
        last_pos = cur_pos+1;
        if (last_pos==std::string::npos || last_pos>=m_str.size()){
            break;
        }
    }
    if (!(last_pos==std::string::npos || last_pos>m_str.size())){
        m_substr_list.push_back(m_str.substr(last_pos));
    }

    return m_substr_list;
}
