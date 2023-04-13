/*
lrst.i: SWIG module defining the Python wrapper for our C++ modules
[Ref.: https://www.swig.org/Doc1.3/Python.html]
 */

%module lrst  // Module name

// Headers
%{
#include <string>
#include <vector>
#include "InputStructure.h"
#include "OutputStructures.h"
#include "ModuleCaller.h"
%}

// typedef long int int64_t;
// %apply long int { int64_t };
%apply long long { int64_t };  // Maps int64_t to long long in Python
%apply unsigned long long { uint64_t };  // Maps uint64_t to unsigned long long in Python
%apply long { long int };  // Maps long int to long in Python

%include <std_string.i>
%include <stdint.i>
%include <std_vector.i>

// These are the header functions wrapped by our lrst module (Like an 'import')
%include "InputStructure.h"  // Contains InputPara for passing parameters to C++
%include "OutputStructures.h"  // Contains data structures for storing statistics for each file type
%include "ModuleCaller.h"  // Functions for calling the C++ statistics computation modules
