#include <iostream>
#include "LRST_function.h"
#include "BAM_module.h"
#include "F5_module.h"
#include "FASTQ_module.h"
#include "FASTA_module.h"

int generate_statistic_from_bam(Input_Para &_input_data, Output_BAM &py_output_bam)
{
   BAM_Module _bam_mod(_input_data);
   _bam_mod.bam_st(py_output_bam);
   return 0;
}

int generate_statistic_from_fq(Input_Para &_input_data, Output_FQ &py_output_fq)
{
   qc_fastq_files(_input_data, py_output_fq);
   return 0;
}

int generate_statistic_from_fa(Input_Para &_input_data, Output_FA &py_output_fa)
{
//   std::cout << "FA Module entered." << std:: endl;
//   qc_fasta_files(_input_data, py_output_fa);
   return 0;
}

int callFASTAModule(Input_Para& _input_data, Output_FA &py_output_fa)
{
    int exit_code;
    exit_code = qc_fasta_files(_input_data, py_output_fa);
    return exit_code;
}

//int callFASTAModule(Input_Para& _input_data)
//{
//    return 0;
//}

int generate_statistic_from_f5(Input_Para &_input_data, Output_F5 &py_output_f5)
{
   F5_Module _F5_mod(_input_data);
   _F5_mod.F5_st(py_output_f5);
   return 0;
}
