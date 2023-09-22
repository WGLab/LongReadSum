#ifndef BAM_MODULE_H_
#define BAM_MODULE_H_

#include <thread>
#include <mutex>
#include <map>

#include "input_parameters.h"
#include "output_data.h"
#include "hts_reader.h"

class BAM_Module{
public:
    static const size_t records_per_thread = 3000;
    size_t file_index = 0;
    std::map<std::string, bool> secondary_alignment;
    std::map<std::string, bool> supplementary_alignment;

    int run(Input_Para& input_params, Output_BAM& final_output);
    int calculateStatistics(Input_Para& input_params, Output_BAM& final_output);
    static void batchStatistics(HTSReader& reader, int batch_size, Input_Para& input_params, Output_BAM& ref_output, std::mutex& bam_mutex, std::mutex& output_mutex, std::mutex& cout_mutex);

    // RRMS
    // Read the RRMS CSV file and store the read IDs (accepted or rejected)
    std::unordered_set<std::string> readRRMSFile(std::string rrms_csv_file, bool accepted_reads);
};

#endif
