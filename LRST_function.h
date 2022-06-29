#include "Python_to_CXX_class.h"
#include "CXX_to_python_class.h"


int generate_statistic_from_bam( Input_Para& _input_data, Output_BAM& py_output_bam );

int generate_statistic_from_fq( Input_Para& _input_data, Output_FQ& py_output_fq );

int generate_statistic_from_fa( Input_Para& _input_data, Output_FA& py_output_fa );

//int callFASTAModule(Input_Para& _input_data);
int callFASTAModule(Input_Para& _input_data, Output_FA& py_output_fa);

int generate_statistic_from_f5( Input_Para& _input_data, Output_F5& py_output_f5 );
