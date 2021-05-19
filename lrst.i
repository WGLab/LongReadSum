%module lrst
%{
#include "Python_to_CXX_class.h"
#include "CXX_to_python_class.h"
#include "LRST_function.h"
#include <string>
%}

%include <std_string.i>
%include <stdint.i>

%include "Python_to_CXX_class.h"
%include "CXX_to_python_class.h"
%include "LRST_function.h"

