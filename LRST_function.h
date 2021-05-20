#include "BAM_module.h"
#include "Python_to_CXX_class.h"

int generate_statistic_from_bam( Input_Para& _input_data, Output_BAM& py_output_bam );


int generate_statistic_from_fq( Input_Para& _input_data, Output_FQ& py_output_fq );


int generate_statistic_from_fa( Input_Para& _input_data, Output_FA& py_output_fa );


int generate_statistic_from_f5( Input_Para& _input_data, Output_F5& py_output_f5 );


/*
 * class Output_Info{
class Basic_Seq_Statistics{
class Output_FA : public Output_Info {
class Basic_Seq_Quality_Statistics{
class Output_FQ: public Output_FA {
class Output_BAM: public Output_FQ {
class Basic_F5_Statistics{
class Output_F5 : public Output_Info {
*/

