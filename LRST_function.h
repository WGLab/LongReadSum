#include "Python_to_CXX_class.h"
#include "CXX_to_python_class.h"


int callBAMModule( Input_Para& _input_data, Output_BAM& py_output_bam );

int callFASTQModule( Input_Para& _input_data, Output_FQ& py_output_fq );

int callFASTAModule(Input_Para& _input_data, Output_FA& py_output_fa);

int callFAST5Module( Input_Para& _input_data, Output_F5& py_output_f5 );
