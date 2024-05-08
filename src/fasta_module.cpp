/*
FASTA_module.cpp:

*/
#include <stdio.h>
#include <stdlib.h>
// #include <zlib.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include <fstream>

#include "fasta_module.h"


// Helper function for saving formatted read statistics to an output text file
static int qc1fasta(const char *input_file, Output_FA &py_output_fa, FILE *read_details_fp)
{
    int exit_code = 0;
    // gzFile input_fp;
    FILE *input_fp;
    // kseq_t *seq;
    // kseq_t *seq;
    char *read_seq;
    char *read_name;
    double read_gc_cnt;

    Basic_Seq_Statistics &long_read_info = py_output_fa.long_read_info;
    std::ifstream input_file_stream(input_file);
    // input_fp = fopen(input_file, "r");
    // if (!input_fp)
    // {
    //     fprintf(stderr, "\nERROR! Failed to open file for reading: %s\n", input_file);
    //     exit_code = 3;
    // } else {
    if (!input_file_stream.is_open()) {
        fprintf(stderr, "Failed to open file for reading: %s\n", input_file);
        exit_code = 3;
    } else {
        std::string line, sequence, sequence_data_str;
        int gc_count = 0;
        while (std::getline(input_file_stream, line))
        {
            if (line.empty())
            {
                continue;
            }

            if (line[0] == '>') {  // Header line
                // Save the last sequence's statistics
                if (!sequence.empty())
                {
                    // Get the base counts
                    uint64_t base_count = 0;
                    uint64_t n_count = 0;
                    for (int i = 0; i < sequence.size(); i++)
                    {
                        if (sequence[i] == 'A' || sequence[i] == 'a')
                        {
                            long_read_info.total_a_cnt += 1;
                            base_count += 1;
                        }
                        else if (sequence[i] == 'G' || sequence[i] == 'g')
                        {
                            long_read_info.total_g_cnt += 1;
                            gc_count += 1;
                            base_count += 1;
                        }
                        else if (sequence[i] == 'C' || sequence[i] == 'c')
                        {
                            long_read_info.total_c_cnt += 1;
                            gc_count += 1;
                            base_count += 1;
                        }
                        else if (sequence[i] == 'T' || sequence[i] == 't' || sequence[i] == 'U' || sequence[i] == 'u')
                        {
                            long_read_info.total_tu_cnt += 1;
                            base_count += 1;
                        } else {
                            n_count += 1;
                        }
                    }

                    // Save sequence length statistics
                    if (base_count > long_read_info.longest_read_length)
                    {
                        long_read_info.longest_read_length = base_count;
                    }
                    long_read_info.total_num_reads += 1;

                    // Get the sequence length distribution
                    // long_read_info.total_num_bases += (uint64_t)read_len;
                    if (base_count < long_read_info.read_length_count.size()) {
                        long_read_info.read_length_count[(int)base_count] += 1;
                    } else {
                        long_read_info.read_length_count.resize(base_count + 1000, 0);
                        long_read_info.read_length_count[(int)base_count] += 1;
                    }

                    // long_read_info.total_n_cnt += read_len - base_count;
                    long_read_info.total_num_bases += base_count;
                    long_read_info.total_n_cnt += n_count;
                    read_gc_cnt = 100.0 * gc_count / (double)base_count;
                    long_read_info.read_gc_content_count[(int)(read_gc_cnt + 0.5)] += 1;
                    fprintf(read_details_fp, "%s\t%ld\t%.2f\n", sequence_data_str, base_count, read_gc_cnt);
                    sequence.clear();
                    gc_count = 0;
                }

                // Update the new sequence's name
                sequence_data_str = line.substr(1);

                // Clear the sequence
                sequence.clear();

            } else {
                sequence += line;  // T
            }

        }

        // Save the last sequence's statistics
        if (!sequence.empty())
        {
            // Get the base counts
            uint64_t base_count = 0;
            uint64_t n_count = 0;
            for (int i = 0; i < sequence.size(); i++)
            {
                if (sequence[i] == 'A' || sequence[i] == 'a')
                {
                    long_read_info.total_a_cnt += 1;
                    base_count += 1;
                }
                else if (sequence[i] == 'G' || sequence[i] == 'g')
                {
                    long_read_info.total_g_cnt += 1;
                    gc_count += 1;
                    base_count += 1;
                }
                else if (sequence[i] == 'C' || sequence[i] == 'c')
                {
                    long_read_info.total_c_cnt += 1;
                    gc_count += 1;
                    base_count += 1;
                }
                else if (sequence[i] == 'T' || sequence[i] == 't' || sequence[i] == 'U' || sequence[i] == 'u')
                {
                    long_read_info.total_tu_cnt += 1;
                    base_count += 1;
                } else {
                    n_count += 1;
                }
            }

            // Save sequence length statistics
            if (base_count > long_read_info.longest_read_length)
            {
                long_read_info.longest_read_length = base_count;
            }
            long_read_info.total_num_reads += 1;

            // Get the sequence length distribution
            // long_read_info.total_num_bases += (uint64_t)read_len;
            if (base_count < long_read_info.read_length_count.size()) {
                long_read_info.read_length_count[(int)base_count] += 1;
            } else {
                long_read_info.read_length_count.resize(base_count + 1000, 0);
                long_read_info.read_length_count[(int)base_count] += 1;
            }

            // long_read_info.total_n_cnt += read_len - base_count;
            long_read_info.total_num_bases += base_count;
            long_read_info.total_n_cnt += n_count;
            read_gc_cnt = 100.0 * gc_count / (double)base_count;
            long_read_info.read_gc_content_count[(int)(read_gc_cnt + 0.5)] += 1;
            fprintf(read_details_fp, "%s\t%ld\t%.2f\n", sequence_data_str, base_count, read_gc_cnt);
            sequence.clear();
            gc_count = 0;
        }

        // Close the input file
        input_file_stream.close();
    }

    return exit_code;
}

// Save summary statistics to the output file
int qc_fasta_files(Input_Para &_input_data, Output_FA &py_output_fa)
{
    int exit_code = 0;
    const char *input_file = NULL;
    std::string read_details_file, read_summary_file;
    FILE *read_details_fp, *read_summary_fp;

    read_details_file = _input_data.output_folder + "/FASTA_details.txt";
    read_summary_file = _input_data.output_folder + "/FASTA_summary.txt";

    // =============

    py_output_fa.long_read_info.total_num_reads = ZeroDefault; // total number of long reads
    py_output_fa.long_read_info.total_num_bases = ZeroDefault; // total number of bases

    py_output_fa.long_read_info.longest_read_length = ZeroDefault; // the length of longest reads
    py_output_fa.long_read_info.n50_read_length = MoneDefault;     // N50
    py_output_fa.long_read_info.n95_read_length = MoneDefault;     // N95
    py_output_fa.long_read_info.n05_read_length = MoneDefault;     // N05;
    py_output_fa.long_read_info.mean_read_length = MoneDefault;    // mean of read length

    py_output_fa.long_read_info.NXX_read_length.clear();
    py_output_fa.long_read_info.median_read_length = MoneDefault; // median of read length

    py_output_fa.long_read_info.total_a_cnt = ZeroDefault;  // A content
    py_output_fa.long_read_info.total_c_cnt = ZeroDefault;  // C content
    py_output_fa.long_read_info.total_g_cnt = ZeroDefault;  // G content
    py_output_fa.long_read_info.total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
    py_output_fa.long_read_info.total_n_cnt = ZeroDefault;  // N content
    py_output_fa.long_read_info.gc_cnt = ZeroDefault;       // GC ratio

    // =============

    //int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

    py_output_fa.long_read_info.read_gc_content_count.clear();
    py_output_fa.long_read_info.read_length_count.clear();

    py_output_fa.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
    // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

    py_output_fa.long_read_info.read_gc_content_count.resize(101, 0);
    // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

    py_output_fa.long_read_info.NXX_read_length.resize(101, 0);
    // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;

    // =============

    read_details_fp = fopen(read_details_file.c_str(), "w");
    if (NULL == read_details_fp)
    {
        fclose(read_details_fp);
        std::cerr << "Failed to write output file: " << read_details_file << std::endl;
        exit_code = 3;
    } else {
        fprintf(read_details_fp, "#read_name\tlength\tGC\n");

        // Check that the input files are valid
        for (size_t i = 0; i < _input_data.num_input_files; i++)
        {
            input_file = _input_data.input_files[i].c_str();
            exit_code = qc1fasta(input_file, py_output_fa, read_details_fp);
        }
        fclose(read_details_fp);

        // Calculate statistics if the input files are valid
        if (exit_code == 0)
        {
            double g_c = py_output_fa.long_read_info.total_g_cnt + py_output_fa.long_read_info.total_c_cnt;
            double a_tu_g_c = g_c + py_output_fa.long_read_info.total_a_cnt + py_output_fa.long_read_info.total_tu_cnt;

            // Calculate read length statistics if base counts are not zero
            uint64_t total_num_bases = py_output_fa.long_read_info.total_num_bases;
            if (total_num_bases == 0) {
                std::cerr << "No bases found in input files." << std::endl;
                exit_code = 3;
            } else {
                py_output_fa.long_read_info.gc_cnt = g_c / a_tu_g_c;

                int percent = 1;
                uint64_t num_bases_sum = 0;
                int64_t num_reads_sum = 0;
                py_output_fa.long_read_info.median_read_length = -1;
                for (int read_len = py_output_fa.long_read_info.read_length_count.size() - 1; read_len > 0; read_len--)
                {
                    num_reads_sum += py_output_fa.long_read_info.read_length_count[read_len];
                    num_bases_sum += (uint64_t) (py_output_fa.long_read_info.read_length_count[read_len] * read_len);
                    if (num_reads_sum * 2 > py_output_fa.long_read_info.total_num_reads && py_output_fa.long_read_info.median_read_length < 0)
                    {
                        py_output_fa.long_read_info.median_read_length = read_len;
                    }
                    if (num_bases_sum * 100 > py_output_fa.long_read_info.total_num_bases * percent)
                    {
                        py_output_fa.long_read_info.NXX_read_length[percent] = read_len;
                        percent += 1;
                        if (percent > 100)
                        {
                            break;
                        }
                    }
                }

                py_output_fa.long_read_info.n50_read_length = py_output_fa.long_read_info.NXX_read_length[50];
                py_output_fa.long_read_info.n95_read_length = py_output_fa.long_read_info.NXX_read_length[95];
                py_output_fa.long_read_info.n05_read_length = py_output_fa.long_read_info.NXX_read_length[5];
                py_output_fa.long_read_info.mean_read_length = (double)py_output_fa.long_read_info.total_num_bases / (double)py_output_fa.long_read_info.total_num_reads;

                read_summary_fp = fopen(read_summary_file.c_str(), "w");
                fprintf(read_summary_fp, "total number of reads\t%d\n", py_output_fa.long_read_info.total_num_reads);
                fprintf(read_summary_fp, "total number of bases\t%ld\n", py_output_fa.long_read_info.total_num_bases);
                fprintf(read_summary_fp, "longest read length\t%d\n", py_output_fa.long_read_info.longest_read_length);
                fprintf(read_summary_fp, "N50 read length\t%d\n", py_output_fa.long_read_info.n50_read_length);
                fprintf(read_summary_fp, "mean read length\t%.2f\n", py_output_fa.long_read_info.mean_read_length);
                fprintf(read_summary_fp, "median read length\t%d\n", py_output_fa.long_read_info.median_read_length);
                fprintf(read_summary_fp, "GC%%\t%.2f\n", py_output_fa.long_read_info.gc_cnt * 100);
                fprintf(read_summary_fp, "\n\n");
                for (int percent = 5; percent < 100; percent += 5)
                {
                    fprintf(read_summary_fp, "N%02d read length\t%.d\n", percent, py_output_fa.long_read_info.NXX_read_length[percent]);
                }

                fprintf(read_summary_fp, "\n\n");

                fprintf(read_summary_fp, "GC content\tnumber of reads\n");
                for (int gc_ratio = 0; gc_ratio <= 100; gc_ratio++)
                {
                    fprintf(read_summary_fp, "GC=%d%%\t%d\n", gc_ratio, py_output_fa.long_read_info.read_gc_content_count[gc_ratio]);
                }
                fclose(read_summary_fp);
            }
        }
    }
    return exit_code;
}
