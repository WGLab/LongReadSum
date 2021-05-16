#include "LRST_function.h"

Output_BAM generate_statistic_from_bam( Input_Para& _input_data ){
   BAM_Module _bam_mod( _input_data );
   return _bam_mod.bam_st();
}

