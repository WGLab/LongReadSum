%module lrst
%{
#include "CXX_to_python_class.h"
#include "Python_to_CXX_class.h"
#include "LRST_function.h"
%}

%include "CXX_to_python_class.h"
%include "Python_to_CXX_class.h"
%include std_string.i
%include stdint.i

Output_BAM generate_statistic_from_bam( Input_Para& _input_data );


