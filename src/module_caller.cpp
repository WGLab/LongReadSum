#include <iostream>
#include "module_caller.h"
#include "bam_module.h"
#include "seqtxt_module.h"
#include "fastq_module.h"
#include "fasta_module.h"
#include "fast5_module.h"
#include "pod5_module.h"

int callBAMModule(Input_Para &_input_data, Output_BAM &py_output_bam)
{
//    BAM_Module _bam_module(_input_data);
    BAM_Module module;
    int exit_code = module.run(_input_data, py_output_bam);
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

int callSeqTxtModule(Input_Para &_input_data, Output_SeqTxt &py_output_SeqTxt)
{
    // Initialize the sequencing_summary.txt module with parameters
    SeqTxt_Module stm(_input_data);
    int exit_code = stm.generateStatistics(py_output_SeqTxt);
    return exit_code;
}

int callFAST5Module(Input_Para &_input_data, Output_FAST5 &py_output_FAST5)
{
    int exit_code = generateQCForFAST5(_input_data, py_output_FAST5);
    return exit_code;
}

int callPOD5Module(Input_Para &_input_data, Output_FAST5 &py_output_FAST5)
{
    int exit_code = generateQCForPOD5(_input_data, py_output_FAST5);
    return exit_code;
}
