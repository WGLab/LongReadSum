#include "Python_to_CXX_class.h"
#include "CXX_to_python_class.h"


int callBAMModule( Input_Para& _input_data, Output_BAM& py_output_bam );

int callFASTQModule( Input_Para& _input_data, Output_FQ& py_output_fq );

int callFASTAModule( Input_Para& _input_data, Output_FA& py_output_fa );

int callSeqTxtModule( Input_Para& _input_data, Output_SeqTxt& py_output_seqtxt );

int callFAST5Module( Input_Para &_input_data, Output_FAST5 &py_output_FAST5 );
