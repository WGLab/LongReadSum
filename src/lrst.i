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
%include "InputStructure.h"  // Contains InputPara for passing parameters to C++
%include "OutputStructures.h"  // Contains data structures for storing statistics for each file type
%include "ModuleCaller.h"  // Functions for calling the C++ statistics computation modules
