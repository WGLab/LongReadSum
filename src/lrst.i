/*
lrst.i: SWIG module defining the Python wrapper for our C++ modules
[Ref.: https://www.swig.org/Doc1.3/Python.html]
 */

%module lrst  // Module name

// Headers
%{
#include <string>
#include <vector>
#include "input_parameters.h"
#include "output_data.h"
#include "module_caller.h"
%}

typedef long int int64_t;

%include <std_string.i>
%include <stdint.i>
%include <std_vector.i>

namespace std{
    %template(IntVector) vector<int>;
    %template(DoubleVector) vector<double>;
    %template(Int64Vector) vector<int64_t>;
    %template(Int2DVector) vector<vector<int>>;
};

// These are the header functions wrapped by our lrst module (Like an 'import')
%include "input_parameters.h"  // Contains InputPara for passing parameters to C++
%include "output_data.h"  // Contains data structures for storing statistics for each file type
%include "module_caller.h"  // Functions for calling the C++ statistics computation modules
