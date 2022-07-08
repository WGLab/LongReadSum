#ifndef COMFUNCTION_H_
#define COMFUNCTION_H_

#include <stdio.h>

#include <map>
#include <vector>
#include <string>

#include <limits>
#include <iomanip>

#include "ComStruct.h"

#define get_array_size(m_a) sizeof(m_a)/sizeof(m_a[0])

// Rounding function used for timing the modules
#define round3(d) (((double)((int)(d*1000+0.5)))/1000)

#define UNUSED(expr) do { (void)(expr); } while (0)

#ifdef WINDOWS
    #include <direct.h>
#else
    #include <unistd.h>
#endif

bool isExpectedLength(char * destination_buffer, const char * input_string, int max_buffer_size);

std::vector<std::string> readRepeatDataFromString(const std::string & m_str, std::string multi_delimiters=WhiteSpace, bool contain_delimiter=false);

#endif
