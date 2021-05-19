#include "LRST_function.h"

/*
Output_BAM generate_statistic_from_bam( Input_Para& _input_data ){
   BAM_Module _bam_mod( _input_data );
   //return _bam_mod.bam_st();
   Output_BAM ob = _bam_mod.bam_st();
   return ob;
}*/

int generate_statistic_from_bam( Input_Para& _input_data, Output_BAM& py_output_bam ){
   BAM_Module _bam_mod( _input_data);
   _bam_mod.bam_st(py_output_bam);
   return 0;
}
