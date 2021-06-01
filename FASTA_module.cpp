#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include "kseq.h"
#include "FASTA_module.h"

KSEQ_INIT(gzFile, gzread) // this is a macro defined in kseq.h

static int qc1fasta(const char *input_file, Output_FA &py_output_fa, FILE *read_details_fp)
{
    gzFile input_fp;
    kseq_t *seq;
    char *read_seq;
    char *read_name;
    int read_len;
    double read_gc_cnt;

    Basic_Seq_Statistics &long_read_info = py_output_fa.long_read_info;
    input_fp = gzopen(input_file, "r");
    if (!input_fp)
    {
        fprintf(stderr, "ERROR! Failed to open file for reading: %s", input_file);
        exit(1);
    }
    seq = kseq_init(input_fp);
    while ((read_len = kseq_read(seq)) >= 0)
    {
        if (read_len == 0) {continue;}
        read_name = seq->name.s;
        read_seq = seq->seq.s;
        if ((uint64_t)read_len > long_read_info.longest_read_length)
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
        read_gc_cnt = 0;

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
        }
        read_gc_cnt = 100.0 * read_gc_cnt / (double)read_len;
        long_read_info.read_gc_content_count[(int)(read_gc_cnt + 0.5)] += 1;
        fprintf(read_details_fp, "%s\t%d\t%.2f\n", read_name, read_len, read_gc_cnt);
    }

    kseq_destroy(seq);
    gzclose(input_fp);

    return 0;
}


int qc_fasta_files(Input_Para &_input_data, Output_FA &py_output_fa)
{
    const char *input_file = NULL;
    std::string read_details_file, read_summary_file;
    FILE *read_details_fp, *read_summary_fp;

    read_details_file = _input_data.output_folder + "/FASTA_details.txt";
    read_summary_file = _input_data.output_folder + "/FASTA_summary.txt";
    
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

    //int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

    py_output_fa.long_read_info.read_gc_content_count.clear();
    py_output_fa.long_read_info.read_length_count.clear();

    py_output_fa.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
    // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

    py_output_fa.long_read_info.read_gc_content_count.resize(101, 0);
    // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

    py_output_fa.long_read_info.NXX_read_length.resize(101, 0);
    // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;


    read_details_fp = fopen(read_details_file.c_str(), "w");
    if (NULL == read_details_fp)
    {
        std::cerr << "ERROR! failed to write output file: " << read_details_file << std::endl;
        exit(1);
    }
    fprintf(read_details_fp, "#read_name\tlength\tGC\n");

    for (size_t i = 0; i < _input_data.num_input_files; i++)
    {
        input_file = _input_data.input_files[i].c_str();
        qc1fasta(input_file, py_output_fa, read_details_fp);
    }
    fclose(read_details_fp);

    double g_c = py_output_fa.long_read_info.total_g_cnt + py_output_fa.long_read_info.total_c_cnt;
    double a_tu_g_c = g_c + py_output_fa.long_read_info.total_a_cnt + py_output_fa.long_read_info.total_tu_cnt;
    if (a_tu_g_c != (double)py_output_fa.long_read_info.total_num_bases)
    {
        fprintf(stderr, "ERROR! total_num_bases is not consistent! this is a bug!");
        exit(1);
    }
    py_output_fa.long_read_info.gc_cnt = g_c / a_tu_g_c;

    int percent = 1;
    int64_t num_bases_sum = 0;
    int64_t num_reads_sum = 0;
    py_output_fa.long_read_info.median_read_length = -1;
    for (int read_len = py_output_fa.long_read_info.read_length_count.size() - 1; read_len > 0; read_len--)
    {
        num_reads_sum += py_output_fa.long_read_info.read_length_count[read_len];
        num_bases_sum += py_output_fa.long_read_info.read_length_count[read_len] * read_len;
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
    fprintf(read_summary_fp, "total number of reads\t%ld\n", py_output_fa.long_read_info.total_num_reads);
    fprintf(read_summary_fp, "total number of bases\t%ld\n", py_output_fa.long_read_info.total_num_bases);
    fprintf(read_summary_fp, "longest read length\t%lu\n", py_output_fa.long_read_info.longest_read_length);
    fprintf(read_summary_fp, "N50 read length\t%ld\n", py_output_fa.long_read_info.n50_read_length);
    fprintf(read_summary_fp, "mean read length\t%.2f\n", py_output_fa.long_read_info.mean_read_length);
    fprintf(read_summary_fp, "median read length\t%ld\n", py_output_fa.long_read_info.median_read_length);
    fprintf(read_summary_fp, "GC%%\t%.2f\n", py_output_fa.long_read_info.gc_cnt * 100);
    fprintf(read_summary_fp, "\n\n");
    for (int percent = 5; percent < 100; percent += 5)
    {
        fprintf(read_summary_fp, "N%02d read length\t%.ld\n", percent, py_output_fa.long_read_info.NXX_read_length[percent]);
    }

    fprintf(read_summary_fp, "\n\n");

    fprintf(read_summary_fp, "GC content\tnumber of reads\n");
    for (int gc_ratio = 0; gc_ratio <= 100; gc_ratio++)
    {
        fprintf(read_summary_fp, "GC=%d%%\t%ld\n", gc_ratio, py_output_fa.long_read_info.read_gc_content_count[gc_ratio]);
    }

    fclose(read_summary_fp);

    return 0;
}