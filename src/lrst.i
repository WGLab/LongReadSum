/*
lrst.i: SWIG module defining the Python wrapper for our C++ modules
[Ref.: https://www.swig.org/Doc1.3/Python.html]
 */

%module lrst  // Module name

// Headers
%{
#include <string>
#include <vector>
#include "Python_to_CXX_class.h"
#include "CXX_to_python_class.h"
#include "LRST_function.h"
#include "FAST5_module.h"
%}

int runExample(void);

typedef long int int64_t;
%apply long int { int64_t };
%include <std_string.i>
%include <stdint.i>
%include <std_vector.i>

namespace std{
  %template(IntVector) vector<int>;
  %template(DoubleVector) vector<double>;
  %template(Int64Vector) vector<int64_t>;
};

// These are the header functions wrapped by our lrst module (Like an 'import')
%include "Python_to_CXX_class.h"  // Contains InputPara for passing parameters to C++
%include "CXX_to_python_class.h"  // Contains data structures for storing statistics for each file type
%include "LRST_function.h"  // Functions for calling the C++ statistics computation modules
%include "LRST_function.h"  // Functions for calling the C++ statistics computation modules
