/*
F5_module.cpp:
Class for calling FAST5 statistics modules.

*/

#include <iostream>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>

#include "FAST5_module.h"
#include "ComFunction.h"
#include "H5Cpp.h"
using namespace H5;


std::vector<std::string> splitString(const std::string& str)
{
    std::vector<std::string> tokens;

    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, '\n')) {
        tokens.push_back(token);
    }

    return tokens;
}


std::string getFastq(H5::H5File f5, std::string basecall_group)
{
        H5::DataSet   dataset    = f5.openDataSet("/Analyses/"+basecall_group+"/BaseCalled_template/Fastq");
        H5::DataType  datatype   = dataset.getDataType();

        // allocate output
        std::string data;

        // read output
        dataset.read(data, datatype);
        std::vector<std::string> tokens = splitString(data);

    return tokens[1];
}


// Calculate GC content for the file
static int calculateFileGCContent(const char *input_file, char quality_value_offset, Output_FAST5 &output_data, FILE *read_details_fp)
{
    int exit_code = 0;
    gzFile input_fp;
    char *read_seq;
    char *raw_read_qual;
    char *read_name;
    char baseq;
    int read_len;
    double read_gc_cnt;
    double read_mean_base_qual;
    int guanine_cytosine_count = 0;  // Count GC bases
    int total_base_count = 0;  // Count total number of bases
    int percent_gc_content;  // Final percent GC content

    Basic_Seq_Statistics &long_read_info = output_data.long_read_info;
    Basic_Seq_Quality_Statistics &seq_quality_info = output_data.seq_quality_info;
    input_fp = gzopen(input_file, "r");
    if (!input_fp)
    {
        std::cerr << "Failed to open file for reading: " << input_file << std::endl;
        exit_code = 3;
    } else {
        // TODO: Read the FAST5 file here
        // Loop through each read

        // Get percent guanine (G) or cytosine (C) in the read


//        seq = kseq_init(input_fp);
//        while ((read_len = kseq_read(seq)) >= 0)
//        {
//            if (read_len == 0) {continue;}
//            read_name = seq->name.s;
//            read_seq = seq->seq.s;
//            raw_read_qual = seq->qual.s;
//            if ((uint64_t)read_len > long_read_info.longest_read_length)
//            {
//                long_read_info.longest_read_length = read_len;
//            }
//            long_read_info.total_num_reads += 1;
//            long_read_info.total_num_bases += read_len;
//            if ((uint64_t)read_len < long_read_info.read_length_count.size()) {
//                long_read_info.read_length_count[read_len] += 1;
//            } else {
//                long_read_info.read_length_count.resize(read_len + 1000, 0);
//                long_read_info.read_length_count[read_len] += 1;
//            }
//            read_gc_cnt = 0;
//            read_mean_base_qual = 0;
//            for (int i = 0; i < read_len; i++)
//            {
//                if (read_seq[i] == 'A' || read_seq[i] == 'a')
//                {
//                    long_read_info.total_a_cnt += 1;
//                }
//                else if (read_seq[i] == 'G' || read_seq[i] == 'g')
//                {
//                    long_read_info.total_g_cnt += 1;
//                    read_gc_cnt += 1;
//                }
//                else if (read_seq[i] == 'C' || read_seq[i] == 'c')
//                {
//                    long_read_info.total_c_cnt += 1;
//                    read_gc_cnt += 1;
//                }
//                else if (read_seq[i] == 'T' || read_seq[i] == 't' || read_seq[i] == 'U' || read_seq[i] == 'u')
//                {
//                    long_read_info.total_tu_cnt += 1;
//                }
//                baseq = raw_read_qual[i] - quality_value_offset;
//                seq_quality_info.base_quality_distribution[baseq] += 1;
//                read_mean_base_qual += baseq;
//            }
//            read_gc_cnt = 100.0 * read_gc_cnt / (double)read_len;
//            long_read_info.read_gc_content_count[(int)(read_gc_cnt + 0.5)] += 1;
//            read_mean_base_qual /= (double) read_len;
//            seq_quality_info.read_average_base_quality_distribution[(uint)(read_mean_base_qual + 0.5)] += 1;
//            fprintf(read_details_fp, "%s\t%d\t%.2f\t%.2f\n", read_name, read_len, read_gc_cnt, read_mean_base_qual);
//        }
//        kseq_destroy(seq);
    }
    gzclose(input_fp);
    return exit_code;
}

static uint8_t predict_base_quality_offset(Input_Para &_input_data)
{
    const char * input_file;
    gzFile input_fp;
    char *raw_read_qual;
    int read_len;
    char quality_value_offset = 0;
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
            if (quality_value_offset){break;}
            input_file = _input_data.input_files[i].c_str();
            input_fp = gzopen(input_file, "r");
            if (!input_fp)
            {
                std::cerr <<  "Failed to open file for reading: %s" << input_file << std::endl;
            } else {
                // TODO: Read the FAST5 file here
//                seq = kseq_init(input_fp);
//                while ((read_len = kseq_read(seq)) >= 0)
//                {
//                    if (quality_value_offset) {break;}
//                    if (read_len == 0) {continue;}
//                    raw_read_qual = seq->qual.s;
//                    for (int j = 0; j < read_len; j++)
//                    {
//                        if (raw_read_qual[j] < 64) {
//                            quality_value_offset = 33;
//                            break;
//                        }
//                    }
//                    num_processed_reads++;
//                    if (num_processed_reads > max_num_reads){
//                        quality_value_offset = 64;
//                        break;
//                    }
//                }
//                kseq_destroy(seq);
            }
            gzclose(input_fp);
        }
    }

    if (quality_value_offset == 0){
        quality_value_offset = 64;
    }
    std::cout << "NOTICE: predicted FAST5 base quality offset is " << quality_value_offset << std::endl;

    return quality_value_offset;
}

int generateQCForFAST5(Input_Para &_input_data, Output_FAST5 &output_data)
{
    // HDF5 section
    //  ==================
    const char* f5_path= "/home/ahsanm1/repos/new_deepmod/guppy_data/kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch101_read110620_strand.fast5";
    std::string basecall_group = "Basecall_1D_000";
    H5::H5File f5 = H5::H5File(f5_path, H5F_ACC_RDONLY);


    std::string fq = getFastq(f5, basecall_group);
    std::cout << fq << std::endl;

    H5::DataSet signal_ds = f5.openDataSet("/Raw/Reads/Read_110620/Signal");
    H5::DataType mdatatype= signal_ds.getDataType();
    H5::DataSpace dataspace  = signal_ds.getSpace();
    hsize_t rank;
    hsize_t dims[2];
    rank = dataspace.getSimpleExtentDims(dims, NULL); // rank = 1
    std::cout<<"Datasize: "<<dims[0]<<std::endl; // this is the correct number of values
    int16_t f5signals [MAX_F5SIGNAL_SIZE];
    signal_ds.read(f5signals, mdatatype);

    H5::DataSet mv_ds = f5.openDataSet("/Analyses/Basecall_1D_000/BaseCalled_template/Move");
    H5::DataType mvdatatype= mv_ds.getDataType();
    H5::DataSpace mvdataspace  = mv_ds.getSpace();
    hsize_t mvrank;
    hsize_t mvdims[2];

    mvrank = mvdataspace.getSimpleExtentDims(mvdims, NULL); // rank = 1
    std::cout<<"Datasize: "<<mvdims[0]<<std::endl; // this is the correct number of values

    int8_t * f5move;
    f5move = new int8_t[MAX_F5SIGNAL_SIZE];

    mv_ds.read(f5move, mvdatatype, mvdataspace);

    for(int i = 0; i < mvdims[0]; i++) {
       std::cout<< *f5move <<" | ";
        f5move++;
    }
    //  ==================

    int exit_code = 0;
    const char *input_file = NULL;
    char quality_value_offset;
    std::string read_details_file, read_summary_file;
    FILE *read_details_fp, *read_summary_fp;

    read_details_file = _input_data.output_folder + "/FAST5_details.txt";
    read_summary_file = _input_data.output_folder + "/FAST5_summary.txt";

    output_data.long_read_info.total_num_reads = ZeroDefault; // total number of long reads
    output_data.long_read_info.total_num_bases = ZeroDefault; // total number of bases

    output_data.long_read_info.longest_read_length = ZeroDefault; // the length of longest reads
    output_data.long_read_info.n50_read_length = MoneDefault;     // N50
    output_data.long_read_info.n95_read_length = MoneDefault;     // N95
    output_data.long_read_info.n05_read_length = MoneDefault;     // N05;
    output_data.long_read_info.mean_read_length = MoneDefault;    // mean of read length

    output_data.long_read_info.NXX_read_length.clear();
    output_data.long_read_info.median_read_length = MoneDefault; // median of read length

    // TODO: Don't need the commented block below
//    output_data.long_read_info.total_a_cnt = ZeroDefault;  // A content
//    output_data.long_read_info.total_c_cnt = ZeroDefault;  // C content
//    output_data.long_read_info.total_g_cnt = ZeroDefault;  // G content
//    output_data.long_read_info.total_tu_cnt = ZeroDefault; // T content for DNA, or U content for RNA
//    output_data.long_read_info.total_n_cnt = ZeroDefault;  // N content
    output_data.long_read_info.gc_cnt = ZeroDefault;       // GC ratio

    //int64_t *read_length_count; // statistics of read length: each position is the number of reads with the length of the index;

    output_data.long_read_info.read_gc_content_count.clear();
    output_data.long_read_info.read_length_count.clear();
    output_data.seq_quality_info.base_quality_distribution.clear();
    output_data.seq_quality_info.read_average_base_quality_distribution.clear();

    output_data.long_read_info.read_length_count.resize(MAX_READ_LENGTH + 1, 0);
    // read_length_count[x] is the number of reads that length is equal to x. MAX_READ_LENGTH is a initial max value, the vector can expand if thre are reads longer than MAX_READ_LENGTH.

    output_data.long_read_info.read_gc_content_count.resize(101, 0);
    // read_gc_content_count[x], x is a integer in the range of [0, 101). read_gc_content_count[x] means number of reads that average GC is x%.

    output_data.long_read_info.NXX_read_length.resize(101, 0);
    // NXX_read_length[50] means N50 read length; NXX_read_length[95] means N95 read length;

    output_data.seq_quality_info.base_quality_distribution.resize(256, 0);
    // base_quality_distribution[x] means number of bases that quality = x.

    output_data.seq_quality_info.read_average_base_quality_distribution.resize(256, 0);
    // base_quality_distribution[x] means number of reads that average base quality = x.

    if (_input_data.user_defined_fastq_base_qual_offset > 0) {
        quality_value_offset = _input_data.user_defined_fastq_base_qual_offset;
    } else {
        // TODO: Implement below function for FAST5
        ;
        //quality_value_offset = predict_base_quality_offset(_input_data);
    }

    read_details_fp = fopen(read_details_file.c_str(), "w");
    if (NULL == read_details_fp)
    {
        std::cerr << "Failed to write output file: " << read_details_file << std::endl;
        exit_code = 3;
    } else {
        fprintf(read_details_fp, "#read_name\tlength\tGC\taverage_baseq_quality\n");

        // Loop through each input file and get the QC data
        for (size_t i = 0; i < _input_data.num_input_files; i++)
        {
            input_file = _input_data.input_files[i].c_str();
            exit_code = calculateFileGCContent(input_file, quality_value_offset, output_data, read_details_fp);
        }
        fclose(read_details_fp);

        // Check if the GC content was calculated successfully
        if (exit_code != 0) {

            // Add the G + C bases
            double g_c = output_data.long_read_info.total_g_cnt + output_data.long_read_info.total_c_cnt;

            // Add all bases
            double a_tu_g_c = g_c + output_data.long_read_info.total_a_cnt + output_data.long_read_info.total_tu_cnt;

            // Check that our total base counts match what was stored (That our code works)
            if (a_tu_g_c != (double)output_data.long_read_info.total_num_bases)
            {
                std::cerr << "Total number of bases is not consistent." << std::endl;
                exit_code = 4;
            } else {
                // Calculate GC-content
                output_data.long_read_info.gc_cnt = g_c / a_tu_g_c;

                int percent = 1;
                int64_t num_bases_sum = 0;
                int64_t num_reads_sum = 0;
                output_data.long_read_info.median_read_length = -1;
                for (int read_len = output_data.long_read_info.read_length_count.size() - 1; read_len > 0; read_len--)
                {
                    num_reads_sum += output_data.long_read_info.read_length_count[read_len];
                    num_bases_sum += output_data.long_read_info.read_length_count[read_len] * read_len;
                    if (num_reads_sum * 2 > output_data.long_read_info.total_num_reads && output_data.long_read_info.median_read_length < 0)
                    {
                        output_data.long_read_info.median_read_length = read_len;
                    }
                    if (num_bases_sum * 100 > output_data.long_read_info.total_num_bases * percent)
                    {
                        output_data.long_read_info.NXX_read_length[percent] = read_len;
                        percent += 1;
                        if (percent > 100)
                        {
                            break;
                        }
                    }
                }

                output_data.long_read_info.n50_read_length = output_data.long_read_info.NXX_read_length[50];
                output_data.long_read_info.n95_read_length = output_data.long_read_info.NXX_read_length[95];
                output_data.long_read_info.n05_read_length = output_data.long_read_info.NXX_read_length[5];
                output_data.long_read_info.mean_read_length = (double)output_data.long_read_info.total_num_bases / (double)output_data.long_read_info.total_num_reads;

                read_summary_fp = fopen(read_summary_file.c_str(), "w");
                fprintf(read_summary_fp, "total number of reads\t%ld\n", output_data.long_read_info.total_num_reads);
                fprintf(read_summary_fp, "total number of bases\t%ld\n", output_data.long_read_info.total_num_bases);
                fprintf(read_summary_fp, "longest read length\t%lu\n", output_data.long_read_info.longest_read_length);
                fprintf(read_summary_fp, "N50 read length\t%ld\n", output_data.long_read_info.n50_read_length);
                fprintf(read_summary_fp, "mean read length\t%.2f\n", output_data.long_read_info.mean_read_length);
                fprintf(read_summary_fp, "median read length\t%ld\n", output_data.long_read_info.median_read_length);
                fprintf(read_summary_fp, "GC%%\t%.2f\n", output_data.long_read_info.gc_cnt * 100);
                fprintf(read_summary_fp, "\n\n");
                for (int percent = 5; percent < 100; percent += 5)
                {
                    fprintf(read_summary_fp, "N%02d read length\t%.ld\n", percent, output_data.long_read_info.NXX_read_length[percent]);
                }

                fprintf(read_summary_fp, "\n\n");

                fprintf(read_summary_fp, "GC content\tnumber of reads\n");
                for (int gc_ratio = 0; gc_ratio < 100; gc_ratio++)
                {
                    fprintf(read_summary_fp, "GC=%d%%\t%ld\n", gc_ratio, output_data.long_read_info.read_gc_content_count[gc_ratio]);
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
                    fprintf(read_summary_fp, "%d\t%ld\n", baseq, output_data.seq_quality_info.read_average_base_quality_distribution[baseq]);
                }
                fclose(read_summary_fp);
            }
        }
    }

    return exit_code;
}
