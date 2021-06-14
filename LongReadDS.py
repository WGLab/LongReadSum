#!/usr/bin/env python

import lrst_global
import os
import sys
import string
import math
from glob import glob
import argparse
import logging
from argparse import RawTextHelpFormatter

import faulthandler
faulthandler.enable()


def get_log_level(log_l):
    if log_l == lrst_global.LOG_ALL:
        return logging.ALL
    elif log_l == lrst_global.LOG_DEBUG:
        return logging.DEBUG
    elif log_l == lrst_global.LOG_INFO:
        return logging.INFO
    elif log_l == lrst_global.LOG_WARN:
        return logging.WARN
    elif log_l == lrst_global.LOG_ERROR:
        return logging.ERROR
    elif log_l == lrst_global.LOG_FATAL:
        return logging.FATAL
    elif log_l == lrst_global.LOG_OFF:
        return logging.OFF
    else:
        return logging.ERROR


def get_common_param(margs, para_dict):
    this_error_str = ""

    if (margs.input == None or margs.input == "") and (margs.inputs == None or margs.inputs == "") and (margs.inputPattern == None or margs.inputPattern == ""):
        this_error_str += "No input file(s) are provided. \n"
    else:
        para_dict["input_files"] = []
        if not (margs.input == None or margs.input == ""):
            para_dict["input_files"].append(margs.input)
        if not (margs.inputs == None or margs.inputs == ""):
            para_dict["input_files"].extend(margs.inputs.split(','))
        if not (margs.inputPattern == None or margs.inputPattern == ""):
            pat_split = margs.inputPattern.split("*")
            para_dict["input_files"].extend(
                glob(os.path.join("*".join(pat_split[:-1]), "*"+pat_split[-1])))
        if len(para_dict["input_files"]) == 0:
            this_error_str += "No input file(s) can be found. \n"
        else:
            for _af in para_dict["input_files"]:
                if not os.path.isfile(_af):
                    this_error_str += "Cannot find "+_af+". \n"

    if (margs.outputfolder == None or margs.outputfolder == ""):
        this_error_str += "No output file is provided. \n"
    else:
        para_dict["output_folder"] = margs.outputfolder
        try:
            if not os.path.isdir(para_dict["output_folder"]):
                os.makedirs(para_dict["output_folder"])
            if not os.path.isdir(para_dict["output_folder"]+'/'+lrst_global.default_image_path):
                os.makedirs(para_dict["output_folder"]+'/' +
                        lrst_global.default_image_path)

        except OSError as e:
            this_error_str += "Cannot create folder for " + \
                para_dict["output_folder"]+" \n"
    para_dict["out_prefix"] = margs.outprefix

    if (margs.log == None or margs.log == ""):
        this_error_str += "No log file is provided. \n"
    else:
        para_dict["log_file"] = margs.log
        para_dict["log_level"] = margs.Log_level
        logging.basicConfig(filename=margs.log, level=get_log_level(margs.Log_level),
                            filemode='w', format="%(levelname)s: %(message)s")

    para_dict["downsample_percentage"] = margs.downsample_percentage

    para_dict["threads"] = margs.thread

    para_dict["random_seed"] = margs.seed

    para_dict["detail"] = margs.detail;

    return this_error_str


def fq_module(margs):
    para_dict = {}
    errorStr = lrst_global.originalError
    errorStr += get_common_param(margs, para_dict)

    if not errorStr == lrst_global.originalError:
        print(errorStr)
        print("#############################################\n")
        # parser.print_help()
        parser.parse_args(['fq', '--help'])
        sys.exit(1001)
    else:
        logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]))
        para_dict["out_prefix"] += "fq_";
        import lrst
        input_para = lrst.Input_Para()
        input_para.threads = para_dict["threads"]
        input_para.rdm_seed = para_dict["random_seed"]
        input_para.downsample_percentage = para_dict["downsample_percentage"]

        input_para.other_flags = 0
        input_para.user_defined_fastq_base_qual_offset = margs.udqual; 

        input_para.output_folder = str(para_dict["output_folder"])
        input_para.out_prefix = str(para_dict["out_prefix"])

        for _ipf in para_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fq_output = lrst.Output_FQ()
        lrst.generate_statistic_from_fq(input_para, fq_output)
        import plot_for_FQ
        plot_for_FQ.fq_plot(fq_output, para_dict)
        import generate_html
        fq_html_gen = generate_html.ST_HTML_Generator(
            [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for FQ", para_dict], static=True)
        fq_html_gen.generate_st_html()
        
        fq_html_gen = generate_html.ST_HTML_Generator(
            [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for FQ", para_dict], static=False)
        fq_html_gen.generate_st_html()
        
        print("Call FQ-module done!")


def fa_module(margs):
    para_dict = {}
    errorStr = lrst_global.originalError
    errorStr += get_common_param(margs, para_dict)

    if not errorStr == lrst_global.originalError:
        print(errorStr)
        print("#############################################\n")
        # parser.print_help()
        parser.parse_args(['fa', '--help'])
        sys.exit(1002)
    else:
        logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]))
        para_dict["out_prefix"] += "fa_";
        import lrst
        input_para = lrst.Input_Para()
        input_para.threads = para_dict["threads"]
        input_para.rdm_seed = para_dict["random_seed"]
        input_para.downsample_percentage = para_dict["downsample_percentage"]

        input_para.other_flags = 0

        input_para.output_folder = str(para_dict["output_folder"])
        input_para.out_prefix = str(para_dict["out_prefix"])

        for _ipf in para_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fa_output = lrst.Output_FA()
        lrst.generate_statistic_from_fa(input_para, fa_output)
        import plot_for_FA
        plot_for_FA.fa_plot(fa_output, para_dict)
        import generate_html
        fa_html_gen = generate_html.ST_HTML_Generator(
            [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for FA", para_dict], static=True)
        fa_html_gen.generate_st_html()
        
        fa_html_gen = generate_html.ST_HTML_Generator(
            [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for FA", para_dict], static=False)
        fa_html_gen.generate_st_html()
        print("Call FA-module done!")


def bam_module(margs):
    para_dict = {}
    errorStr = lrst_global.originalError
    errorStr += get_common_param(margs, para_dict)

    if not errorStr == lrst_global.originalError:
        print(errorStr)
        print("#############################################\n")
        # parser.print_help()
        parser.parse_args(['bam', '--help'])
        sys.exit(1003)
    else:
        logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]))
        para_dict["out_prefix"] += "bam_";
        import lrst
        input_para = lrst.Input_Para()
        input_para.threads = para_dict["threads"]
        input_para.rdm_seed = para_dict["random_seed"]
        input_para.downsample_percentage = para_dict["downsample_percentage"]
        input_para.other_flags =  (1 if para_dict["detail"]>0 else 0) ;
        input_para.output_folder = str(para_dict["output_folder"])
        input_para.out_prefix = str(para_dict["out_prefix"])

        for _ipf in para_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        bam_output = lrst.Output_BAM()
        lrst.generate_statistic_from_bam(input_para, bam_output)
        import plot_for_BAM
        plot_for_BAM.bam_plot(bam_output, para_dict)
        import generate_html
        bam_html_gen = generate_html.ST_HTML_Generator(
            [["basic_st","map_st", "err_st", "read_length_st", "read_length_hist", "base_st", "basic_info"], "The statistics for BAM", para_dict], static=True)
        bam_html_gen.generate_st_html()
        
        bam_html_gen = generate_html.ST_HTML_Generator(
            [["basic_st","map_st", "err_st", "read_length_st", "read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for BAM", para_dict], static=False)
        bam_html_gen.generate_st_html()
        
        
        print("Call BAM-module done!")


def f5_module(margs):
    para_dict = {}
    errorStr = lrst_global.originalError
    errorStr += get_common_param(margs, para_dict)

    if not errorStr == lrst_global.originalError:
        print(errorStr)
        print("#############################################\n")
        # parser.print_help()
        parser.parse_args(['f5', '--help'])
        sys.exit(1004)
    else:
        logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]))
        para_dict["out_prefix"] += "f5_";
        import lrst
        input_para = lrst.Input_Para()
        input_para.threads = para_dict["threads"]
        input_para.rdm_seed = para_dict["random_seed"]
        input_para.downsample_percentage = para_dict["downsample_percentage"]
        input_para.other_flags = margs.seq
        input_para.other_flags = (input_para.other_flags << 4);
        input_para.other_flags += (1 if para_dict["detail"]>0 else 0) ;
        input_para.other_flags = (input_para.other_flags << 4);
        input_para.other_flags += int(margs.sum_type)

        input_para.output_folder = str(para_dict["output_folder"])
        input_para.out_prefix = str(para_dict["out_prefix"])

        for _ipf in para_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        f5_output = lrst.Output_F5()
        lrst.generate_statistic_from_f5(input_para, f5_output)
        import plot_for_F5
        plot_for_F5.f5_plot(f5_output, para_dict)
        import generate_html
        if margs.seq==0:
           f5_html_gen = generate_html.ST_HTML_Generator([["basic_st", "read_length_st","read_length_hist","base_st","basic_info"], "The statistics for F5", para_dict ]);
        else:
           f5_html_gen = generate_html.ST_HTML_Generator(
               [["basic_st", "read_length_st","read_length_hist", "basic_info"], "The statistics for F5", para_dict])
        f5_html_gen.generate_st_html()
        print("Call F5-module done!")


parser = argparse.ArgumentParser(description="Data analysis tools for long-read sequencing data",
                                 epilog="For example, \n \
                                         \tpython %(prog)s fq: with fq or fastq input\n\
                                         \tpython %(prog)s fa: with fa or fasta input\n\
                                         \tpython %(prog)s bam: with bam input\n\
                                         \tpython %(prog)s f5: with fast5 input\n\
                                         ",
                                 formatter_class=RawTextHelpFormatter)


##############################################################################
#
# Different functional modules starting from here.
#
##############################################################################


subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

common_grp_param = parent_parser.add_argument_group(
    "Common parameters for %(prog)s")
#common_grp_param.add_argument("-h", "--help", default="", help="Show this help documents");
input_files_group = common_grp_param.add_mutually_exclusive_group()
input_files_group.add_argument(
    "-i", "--input", type=str, default=None, help="The input file for the analysis")
input_files_group.add_argument(
    "-I", "--inputs", type=str, default=None, help="The input files for the analysis. Files are separated by ','.")
input_files_group.add_argument(
    "-P", "--inputPattern", type=str, default=None, help="The pattern of input files with *. The format is \"patter*n\" where \" is required. ")
input_files_group.add_argument("-p", "--downsample_percentage", type=float, default=1.0,
                               help="The percentage of downsampling for quick run. Default: 1.0 without downsampling")

common_grp_param.add_argument(
    "-g", "--log", type=str, default="log_output.log", help="Log file")
common_grp_param.add_argument("-G", "--Log_level", type=int, default=lrst_global.LOG_ERROR,
                              help="Level for logging: ALL(0) < DEBUG(1) < INFO(2) < WARN(3) < ERROR(4) < FATAL(5) < OFF(6). Default: 4 (ERROR)")

common_grp_param.add_argument("-o", "--outputfolder", type=str,
                              default="output_"+lrst_global.prg_name, help="The output folder.")
common_grp_param.add_argument("-t", "--thread", type=int,
                              default=1, help="The number of threads used. Default: 1.")
common_grp_param.add_argument("-Q", "--outprefix", type=str,
                              default="st_", help="The prefix of output. Default: `st_`.")
common_grp_param.add_argument(
    "-s", "--seed", type=int, default=1, help="The number for random seed. Default: 1.")
common_grp_param.add_argument("-d", "--detail", type=int, default=0, help="Will output detail in files? Default: 0(no).")

fq_parsers = subparsers.add_parser('fq',
                                   parents=[parent_parser],
                                   help="Show data analysis for fq files",
                                   description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                   formatter_class=RawTextHelpFormatter)
fq_parsers.add_argument("-u", "--udqual", type=int, default=-1,
                        help="User defined quality offset for bases in fq. Default: -1.")
fq_parsers.set_defaults(func=fq_module)


fa_parsers = subparsers.add_parser('fa',
                                   parents=[parent_parser],
                                   help="Show data analysis for fa files",
                                   description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                   formatter_class=RawTextHelpFormatter)
fa_parsers.set_defaults(func=fa_module)

bam_parsers = subparsers.add_parser('bam',
                                    parents=[parent_parser],
                                    help="Show data analysis for bam files",
                                    description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                    formatter_class=RawTextHelpFormatter)
bam_parsers.set_defaults(func=bam_module)

f5_parsers = subparsers.add_parser('f5',
                                   parents=[parent_parser],
                                   help="Show data analysis for f5 files",
                                   description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                   formatter_class=RawTextHelpFormatter)
f5_parsers.add_argument("-S", "--seq", type=int, default=1,
                        help="Sequence_summary.txt only? Default: 1(yes).")
f5_parsers.add_argument("-m", "--sum_type", type=int, default=1, choices=[1,2,3],
                        help="Different fields in Sequence_summary.txt. Default: 1.")

f5_parsers.set_defaults(func=f5_module)

if sys.version_info[0] < 2:
    print(lrst_global.prg_name +
          " could not be run with lower version than python 2.7.")
else:
    if sys.version_info[1] < 6:
        print(lrst_global.prg_name+" could be run with python 2.7.")
    else:
        if len(sys.argv) < 2:
            parser.print_help()
        else:
            args = parser.parse_args()
            args.func(args)
