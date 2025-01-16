#include "fastq_module.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>  // std::sort
#include <cmath>  // std::round

#include <fstream>
#include <iostream>
#include <sstream>

#include <sys/stat.h>
#include <sys/types.h>

#include "utils.h"

int qc1fastq(const char *input_file, char fastq_base_qual_offset, Output_FQ &output_data, FILE *read_details_fp)
{
    int exit_code = 0;
    int read_len;
    double read_gc_cnt;
    // double read_mean_base_qual;
    Basic_Seq_Statistics &long_read_info = output_data.long_read_info;
    Basic_Seq_Quality_Statistics &seq_quality_info = output_data.seq_quality_info;
    long_read_info.total_num_reads = ZeroDefault; // total number of long reads
    long_read_info.longest_read_length = ZeroDefault; // the length of longest reads

    std::ifstream input_file_stream(input_file);
    if (!input_file_stream.is_open())
    {
        fprintf(stderr, "Failed to open file for reading: %s\n", input_file);
        exit_code = 3;
    } else {
        std::string line, read_seq, read_name, raw_read_qual;
        while (std::getline(input_file_stream, line) && !line.empty())
        {
            if (line[0] == '@')
            {
                read_name = line.substr(1);
                read_name = read_name.substr(0, read_name.find_first_of(" \t"));
                std::getline(input_file_stream, read_seq);
                std::getline(input_file_stream, line);
                std::getline(input_file_stream, raw_read_qual);
                read_len = read_seq.size();
                if (read_len == 0) {continue;}
                if (read_len > long_read_info.longest_read_length)
                {
                    long_read_info.longest_read_length = read_len;
                }
                long_read_info.total_num_reads += 1;
                long_read_info.total_num_bases += (uint64_t)read_len;

                // Update the read length counts
                if ((uint64_t)read_len < long_read_info.read_length_count.size()) {
                    long_read_info.read_length_count[read_len] += 1;
                } else {
                    long_read_info.read_length_count.resize(read_len + 1000, 0);
                    long_read_info.read_length_count[read_len] += 1;
                }

                // Store the read length
                long_read_info.read_lengths.push_back(read_len);

                // Access base quality data
                // printMessage("[TEST1] Base quality string: " + raw_read_qual);
                char value;
                std::vector<int> base_quality_values;
                // std::string base_quality_str = raw_read_qual;
                std::istringstream iss(raw_read_qual);
                while (iss >> value)
                {
                    int base_quality_value = value - '!';
                    base_quality_values.push_back(base_quality_value);
                    // printMessage("[TEST1] Base quality value: " + std::to_string(base_quality_value));
                }

                // Ensure that the base quality string has the same length as
                // the read sequence
                if (base_quality_values.size() != read_len)
                {
                    printError("Error: Base quality string length does not match read sequence length");
                    exit_code = 1;
                    break;
                }

                // Process base and quality information
                read_gc_cnt = 0;
                // read_mean_base_qual = 0;
                int base_quality_value;
                double cumulative_base_prob = 0;  // Read cumulative base quality probability
                for (int i = 0; i < read_len; i++)
                {
                    if (read_seq[i] == 'A' || read_seq[i] == 'a')
                    {
                        long_read_info.total_a_cnt += 1;
                    }
                    else if (read_seq[i] == 'G' || read_seq[i] == 'g')
                    {
                        long_read_info.total_g_cnt += 1;
                        read_gc_cnt += 1;
                    }
                    else if (read_seq[i] == 'C' || read_seq[i] == 'c')
                    {
                        long_read_info.total_c_cnt += 1;
                        read_gc_cnt += 1;
                    }
                    else if (read_seq[i] == 'T' || read_seq[i] == 't' || read_seq[i] == 'U' || read_seq[i] == 'u')
                    {
                        long_read_info.total_tu_cnt += 1;
                    }

                    // Get the base quality (Phred) value
                    base_quality_value = base_quality_values[i];
                    // base_quality_value = (uint64_t)raw_read_qual[i] - (uint64_t)fastq_base_qual_offset;
                    try {
                        seq_quality_info.base_quality_distribution[base_quality_value] += 1;
                    } catch (const std::out_of_range& oor) {
                        printError("Warning: Base quality value " + std::to_string(base_quality_value) + " exceeds maximum value");
                    }
                    // read_mean_base_qual += (double) base_quality_value;

                    // Convert the Phred quality value to a probability
                    double base_quality_prob = pow(10, -base_quality_value / 10.0);
                    cumulative_base_prob += base_quality_prob;
                }

                // Calculate the mean base quality probability
                cumulative_base_prob /= (double)read_len;

                // Convert the mean base quality probability to a Phred quality
                // value
                double read_mean_base_qual = -10.0 * log10(cumulative_base_prob);
                // printMessage("Mean Q Score for read ID " + read_name + " is " + std::to_string(read_mean_base_qual));

                // Update the per-read GC content distribution
                double gc_content_pct = (100.0 * read_gc_cnt) / static_cast<double>(read_len);
                int gc_content_int = static_cast<int>(std::round(gc_content_pct));
                try {
                    long_read_info.read_gc_content_count[gc_content_int] += 1;
                } catch (const std::out_of_range& oor) {
                    printError("Warning: Invalid GC content value " + std::to_string(gc_content_int));
                }
                
                // Update the per-read base quality distribution
                // double read_mean_base_qual_pct = read_mean_base_qual / static_cast<double>(read_len);
                // unsigned int read_mean_base_qual_int = static_cast<unsigned
                // int>(std::round(read_mean_base_qual_pct));
                int read_mean_base_qual_int = static_cast<int>(std::round(read_mean_base_qual));

                // printMessage("Rounded Mean Q Score for read ID " + read_name + " is " + std::to_string(read_mean_base_qual_int));

                try {
                    seq_quality_info.read_quality_distribution[read_mean_base_qual_int] += 1;
                } catch (const std::out_of_range& oor) {
                    printError("Warning: Base quality value " + std::to_string(read_mean_base_qual_int) + " exceeds maximum value");
                }
                
                // try {
                //     seq_quality_info.read_average_base_quality_distribution[read_mean_base_qual_int] += 1;
                // } catch (const std::out_of_range& oor) {
                //     printError("Warning: Base quality value " + std::to_string(read_mean_base_qual_int) + " exceeds maximum value");
                // }

                fprintf(read_details_fp, "%s\t%d\t%.2f\t%.2f\n", read_name.c_str(), read_len, gc_content_pct, read_mean_base_qual);  // Write to file
            }
        }
        input_file_stream.close();
    }
    return exit_code;
}

int qc_fastq_files(Input_Para &_input_data, Output_FQ &output_data)
{
    int exit_code = 0;
    const char *input_file = NULL;
    char fastq_base_qual_offset;
    std::string read_details_file, read_summary_file;
    FILE *read_details_fp, *read_summary_fp;

    read_details_file = _input_data.output_folder + "/FASTQ_details.txt";
    read_summary_file = _input_data.output_folder + "/FASTQ_summary.txt";
    
    output_data.long_read_info.total_num_reads = ZeroDefault; // total number of long reads
    output_data.long_read_info.total_num_bases = ZeroDefault; // total number of bases

    output_data.long_read_info.longest_read_length = ZeroDefault; // the length of longest reads
    output_data.long_read_info.n50_read_length = MoneDefault;     // N50
    output_data.long_read_info.n95_read_length = MoneDefault;     // N95
    output_data.long_read_info.n05_read_length = MoneDefault;     // N05;
    output_data.long_read_info.mean_read_length = MoneDefault;    // mean of read length

    output_data.long_read_info.NXX_read_length.clear();
    output_data.long_read_info.median_read_length = MoneDefault; // median of read length

    output_data.long_read_info.total_a_cnt = ZeroDefault;  // A content
    output_data.long_read_info.total_c_cnt = ZeroDefault;  // C content
    output_data.long_read_info.total_g_cnt = ZeroDefault;  // G content
    output_data.long_read_info.total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
    output_data.long_read_info.total_n_cnt = ZeroDefault;  // N content
    output_data.long_read_info.gc_cnt = ZeroDefault;       // GC ratio

    output_data.long_read_info.read_gc_content_count.clear();
    output_data.long_read_info.read_length_count.clear();

    output_data.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
    // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

    output_data.long_read_info.read_gc_content_count.resize(101, 0);
    // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

    output_data.long_read_info.NXX_read_length.resize(101, 0);
    // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;

    output_data.seq_quality_info.read_average_base_quality_distribution.resize(MAX_BASE_QUALITY, 0);

    if (_input_data.user_defined_fastq_base_qual_offset > 0) {
        fastq_base_qual_offset = _input_data.user_defined_fastq_base_qual_offset;
    } else {
        fastq_base_qual_offset = 33;
    }

    read_details_fp = fopen(read_details_file.c_str(), "w");

    if (NULL == read_details_fp)
    {
        std::cerr << "Failed to write output file: " << read_details_file << std::endl;
        exit_code = 3;
    } else {
        fprintf(read_details_fp, "#read_name\tlength\tGC\taverage_baseq_quality\n");

        for (size_t i = 0; i < _input_data.num_input_files; i++)
        {
            input_file = _input_data.input_files[i].c_str();
            qc1fastq(input_file, fastq_base_qual_offset, output_data, read_details_fp);
        }
        fclose(read_details_fp);

        // Add the G + C bases
        double g_c = output_data.long_read_info.total_g_cnt + output_data.long_read_info.total_c_cnt;

        // Add all bases
        double a_tu_g_c = g_c + output_data.long_read_info.total_a_cnt + output_data.long_read_info.total_tu_cnt;

        // Calculate read length statistics if base counts are not zero
        uint64_t total_num_bases = output_data.long_read_info.total_num_bases;
        if (total_num_bases == 0) {
            std::cerr << "No bases found in input files." << std::endl;
            exit_code = 3;
        } else {
            // Calculate GC-content
            output_data.long_read_info.gc_cnt = g_c / a_tu_g_c;

            // Get the max read length
            std::vector<int> read_lengths = output_data.long_read_info.read_lengths;
            if (read_lengths.size() == 0) {
                std::cerr << "No reads found in input files." << std::endl;
                exit_code = 3;
                return exit_code;
            } else {
                // Sort the read lengths in descending order
                std::sort(read_lengths.begin(), read_lengths.end(), std::greater<int64_t>());

                // Get the max read length
                int64_t max_read_length = read_lengths.at(0);
                output_data.long_read_info.longest_read_length = max_read_length;
            }

            // Get the median read length
            int64_t median_read_length = read_lengths[read_lengths.size() / 2];
            output_data.long_read_info.median_read_length = median_read_length;

            // Get the mean read length
            float mean_read_length = (double)total_num_bases / (double)read_lengths.size();
            output_data.long_read_info.mean_read_length = mean_read_length;

            // Calculate N50 and other N-scores
            for (int percent_value = 1; percent_value <= 100; percent_value++)
            {
                // Get the base percentage threshold for this N-score
                double base_threshold = (double)total_num_bases * (percent_value / 100.0);

                // Calculate the NXX score
                double current_base_count = 0;
                int current_read_index = -1;
                while (current_base_count < base_threshold) {
                    current_read_index ++;
                    current_base_count += read_lengths.at(current_read_index);
                }
                int nxx_read_length = read_lengths.at(current_read_index);
                output_data.long_read_info.NXX_read_length[percent_value] = nxx_read_length;
            }

            // Set common score variables
            output_data.long_read_info.n50_read_length = output_data.long_read_info.NXX_read_length[50];
            output_data.long_read_info.n95_read_length = output_data.long_read_info.NXX_read_length[95];
            output_data.long_read_info.n05_read_length = output_data.long_read_info.NXX_read_length[5];

            read_summary_fp = fopen(read_summary_file.c_str(), "w");
            fprintf(read_summary_fp, "total number of reads\t%d\n", output_data.long_read_info.total_num_reads);
            fprintf(read_summary_fp, "total number of bases\t%ld\n", output_data.long_read_info.total_num_bases);
            fprintf(read_summary_fp, "longest read length\t%d\n", output_data.long_read_info.longest_read_length);
            fprintf(read_summary_fp, "N50 read length\t%d\n", output_data.long_read_info.n50_read_length);
            fprintf(read_summary_fp, "mean read length\t%.2f\n", output_data.long_read_info.mean_read_length);
            fprintf(read_summary_fp, "median read length\t%d\n", output_data.long_read_info.median_read_length);
            fprintf(read_summary_fp, "GC%%\t%.2f\n", output_data.long_read_info.gc_cnt * 100);
            fprintf(read_summary_fp, "\n\n");
            for (int percent = 5; percent < 100; percent += 5)
            {
                fprintf(read_summary_fp, "N%02d read length\t%.d\n", percent, output_data.long_read_info.NXX_read_length[percent]);
            }

            fprintf(read_summary_fp, "\n\n");

            fprintf(read_summary_fp, "GC content\tnumber of reads\n");
            for (int gc_ratio = 0; gc_ratio < 100; gc_ratio++)
            {
                fprintf(read_summary_fp, "GC=%d%%\t%d\n", gc_ratio, output_data.long_read_info.read_gc_content_count[gc_ratio]);
            }

            fprintf(read_summary_fp, "\n\n");
            fprintf(read_summary_fp, "base quality\tnumber of bases\n");
            for (int baseq = 0; baseq <= 60; baseq++)
            {
                fprintf(read_summary_fp, "%d\t%ld\n", baseq, output_data.seq_quality_info.base_quality_distribution[baseq]);
            }

            fprintf(read_summary_fp, "\n\n");
            fprintf(read_summary_fp, "read average base quality\tnumber of reads\n");
            for (int baseq = 0; baseq <= 60; baseq++)
            {
                fprintf(read_summary_fp, "%d\t%d\n", baseq, output_data.seq_quality_info.read_average_base_quality_distribution[baseq]);
            }
            fclose(read_summary_fp);
        }
    }

    return exit_code;
}
