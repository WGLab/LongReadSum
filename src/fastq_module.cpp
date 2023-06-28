#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <algorithm>

// "Seqtk is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format":
// https://github.com/lh3/seqtk
#include "kseq.h"
#include "kseq.h"
#include "fastq_module.h"

KSEQ_INIT(gzFile, gzread) // this is a macro defined in kseq.h

static int qc1fastq(const char *input_file, char fastq_base_qual_offset, Output_FQ &output_data, FILE *read_details_fp)
{
    int exit_code = 0;
    gzFile input_fp;
    kseq_t *seq;
    char *read_seq;
    char *raw_read_qual;
    char *read_name;
    int read_len;
    double read_gc_cnt;
    double read_mean_base_qual;

    Basic_Seq_Statistics &long_read_info = output_data.long_read_info;
    Basic_Seq_Quality_Statistics &seq_quality_info = output_data.seq_quality_info;
    input_fp = gzopen(input_file, "r");
    if (!input_fp)
    {
        std::cerr << "Failed to open file for reading: " << input_file << std::endl;
        exit_code = 3;
    } else {
        seq = kseq_init(input_fp);
        while ((read_len = kseq_read(seq)) >= 0)
        {
            if (read_len == 0) {continue;}
            read_name = seq->name.s;
            read_seq = seq->seq.s;
            raw_read_qual = seq->qual.s;
            if (read_len > long_read_info.longest_read_length)
            {
                long_read_info.longest_read_length = read_len;
            }
            long_read_info.total_num_reads += 1;
            long_read_info.total_num_bases += read_len;
            if ((uint64_t)read_len < long_read_info.read_length_count.size()) {
                long_read_info.read_length_count[read_len] += 1;
            } else {
                long_read_info.read_length_count.resize(read_len + 1000, 0);
                long_read_info.read_length_count[read_len] += 1;
            }

            // Store the read length
            long_read_info.read_lengths.push_back(read_len);

            read_gc_cnt = 0;
            read_mean_base_qual = 0;
            uint64_t base_quality_value;
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
                base_quality_value = (uint64_t)raw_read_qual[i] - (uint64_t)fastq_base_qual_offset;
                seq_quality_info.base_quality_distribution[base_quality_value] += 1;
                read_mean_base_qual += (double) base_quality_value;
            }
            read_gc_cnt = 100.0 * read_gc_cnt / (double)read_len;
            long_read_info.read_gc_content_count[(int)(read_gc_cnt + 0.5)] += 1;
            read_mean_base_qual /= (double) read_len;
            seq_quality_info.read_average_base_quality_distribution[(uint)(read_mean_base_qual + 0.5)] += 1;
            fprintf(read_details_fp, "%s\t%d\t%.2f\t%.2f\n", read_name, read_len, read_gc_cnt, read_mean_base_qual);
        }
        kseq_destroy(seq);
    }
    gzclose(input_fp);
    return exit_code;
}

static uint8_t predict_base_quality_offset(Input_Para &_input_data)
{
    const char * input_file;
    gzFile input_fp;
    kseq_t *seq;
    char *raw_read_qual;
    int read_len;
    char fastq_base_qual_offset = 0;
    int max_num_reads;
    int num_processed_reads;

    if (_input_data.num_input_files == 0)
    {
        std::cerr << "No input files provided.\n" << std::endl;
    } else {
        max_num_reads = 100000;
        num_processed_reads = 0;
        for (size_t i = 0; i < _input_data.num_input_files; i++)
        {
            if (fastq_base_qual_offset){break;}
            input_file = _input_data.input_files[i].c_str();
            input_fp = gzopen(input_file, "r");
            if (!input_fp)
            {
                std::cerr <<  "Failed to open file for reading: %s" << input_file << std::endl;
            } else {
                seq = kseq_init(input_fp);
                while ((read_len = kseq_read(seq)) >= 0)
                {
                    if (fastq_base_qual_offset) {break;}
                    if (read_len == 0) {continue;}
                    raw_read_qual = seq->qual.s;
                    for (int j = 0; j < read_len; j++)
                    {
                        if (raw_read_qual[j] < 64) {
                            fastq_base_qual_offset = 33;
                            break;
                        }
                    }
                    num_processed_reads++;
                    if (num_processed_reads > max_num_reads){
                        fastq_base_qual_offset = 64;
                        break;
                    }
                }
                kseq_destroy(seq);
            }
            gzclose(input_fp);
        }
    }

    if (fastq_base_qual_offset == 0){
        fastq_base_qual_offset = 64;
    }
//    std::cout << "NOTICE: predicted FASTQ base quality offset is " << fastq_base_qual_offset << std::endl;

    return fastq_base_qual_offset;
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

    //int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

    output_data.long_read_info.read_gc_content_count.clear();
    output_data.long_read_info.read_length_count.clear();
    //output_data.seq_quality_info.base_quality_distribution.clear();
    //output_data.seq_quality_info.read_average_base_quality_distribution.clear();

    output_data.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
    // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

    output_data.long_read_info.read_gc_content_count.resize(101, 0);
    // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

    output_data.long_read_info.NXX_read_length.resize(101, 0);
    // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;

    //output_data.seq_quality_info.base_quality_distribution.resize(256, 0);
    // base_quality_distribution[x] means number of bases that quality = x.

    output_data.seq_quality_info.read_average_base_quality_distribution.resize(256, 0);
    // base_quality_distribution[x] means number of reads that average base quality = x.

    if (_input_data.user_defined_fastq_base_qual_offset > 0) {
        fastq_base_qual_offset = _input_data.user_defined_fastq_base_qual_offset;
    } else {
        fastq_base_qual_offset = predict_base_quality_offset(_input_data);
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

            // Sort the read lengths in descending order
            std::vector<int> read_lengths = output_data.long_read_info.read_lengths;
            std::sort(read_lengths.begin(), read_lengths.end(), std::greater<int64_t>());

            // Get the max read length
            int64_t max_read_length = read_lengths.at(0);
            output_data.long_read_info.longest_read_length = max_read_length;

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