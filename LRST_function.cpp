#include <iostream>
#include "LRST_function.h"
#include "BAM_module.h"
#include "F5_module.h"
#include "FASTQ_module.h"
#include "FASTA_module.h"

int callBAMModule(Input_Para &_input_data, Output_BAM &py_output_bam)
{
   BAM_Module _bam_module(_input_data);
   int exit_code = _bam_module.calculateStatistics(py_output_bam);
   return exit_code;
}

int callFASTQModule(Input_Para &_input_data, Output_FQ &py_output_fq)
{
   int exit_code = qc_fastq_files(_input_data, py_output_fq);
   return exit_code;
}

int callFASTAModule(Input_Para& _input_data, Output_FA &py_output_fa)
{
    int exit_code = qc_fasta_files(_input_data, py_output_fa);
    return exit_code;
}

int callFAST5Module(Input_Para &_input_data, Output_F5 &py_output_f5)
{
    // Initialize the FAST5 module with parameters
    F5_Module _F5_mod(_input_data);
    int exit_code = _F5_mod.F5_st(py_output_f5);
    return exit_code;
}
