#!/usr/bin/env python
"""
cli.py:
Parse arguments and run the filetype-specific module.
"""

import os
import sys
from glob import glob
import argparse
import logging
from argparse import RawTextHelpFormatter
from src import generate_html
from src import lrst_global
from lib import lrst

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
    """
    Format a dict object para_dict to contain our input files/output folder parameters
    Output a string containing all parse argument errors.
    Inputs:
    margs: Command line input arguments (dict).
    """
    # TODO: The para_dict input is always empty in this function. Thus, it can be defined in this function.
    this_error_str = ""

    if (margs.input == None or margs.input == "") and (margs.inputs == None or margs.inputs == "") and (margs.inputPattern == None or margs.inputPattern == ""):
        this_error_str += "No input file(s) are provided. \n"
    else:
        # Populate our
        para_dict["input_files"] = []
        if not (margs.input == None or margs.input == ""):
            para_dict["input_files"].append(margs.input)
        if not (margs.inputs == None or margs.inputs == ""):
            file_list = [file_str.strip() for file_str in margs.inputs.split(',')]
            para_dict["input_files"].extend(file_list)
        if not (margs.inputPattern == None or margs.inputPattern == ""):
            pat_split = margs.inputPattern.split("*")
            para_dict["input_files"].extend(
                glob(os.path.join("*".join(pat_split[:-1]), "*"+pat_split[-1])))
        if len(para_dict["input_files"]) == 0:
            this_error_str += "No input file(s) can be found. \n"
        else:
            for input_filepath in para_dict["input_files"]:
                if not os.path.isfile(input_filepath):
                    this_error_str += "Cannot find the input file: " + input_filepath + "\n"

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

        # Import the SWIG Python wrapper for our C++ module
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
        exit_code = lrst.callFASTQModule(input_para, fq_output)
        if exit_code == 0:
            from src import plot_for_FQ
            plot_for_FQ.fq_plot(fq_output, para_dict)
            for static in [True, False]:
                fq_html_gen = generate_html.ST_HTML_Generator(
                    [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for FQ", para_dict], static=static)
                fq_html_gen.generate_st_html()

            print("Call FQ-module done!")


def fa_module(margs):
    """
    Run the FASTA filetype module.
    """
    para_dict = {}
    errorStr = lrst_global.originalError

    # Format a dict object to contain our input files/output folder parameters
    errorStr += get_common_param(margs, para_dict)  # TODO: para_dict is always an empty dict in these calls

    # If original error string has not been appended to by get_common_params(), this means that there were no parse
    # errors. TODO: originalError can be an empty string in lrst_global.py
    #  TODO: Define these global variables at the top of this file (or module function) instead
    if not errorStr == lrst_global.originalError:
        # If there are parse errors, display the filetype-specific help instructions
        print(errorStr)
        print("#############################################\n")
        # parser.print_help()
        parser.parse_args(['fa', '--help'])
        sys.exit(1002)
    else:
        # If there are no parse errors, run the filetype-specific module
        logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]))
        para_dict["out_prefix"] += "fa_";

        # TODO: Multiple calls to the Python-wrapped C++ module are made (the lrst calls), but only a single
        #  module call to the filetype-specific function should be needed. Could improve performance
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
        exit_code = lrst.callFASTAModule(input_para, fa_output)
        if exit_code == 0:
            from src import plot_for_FA
            plot_for_FA.fa_plot(fa_output, para_dict)

            # TODO: Unused 'static' variable results in redundant function call
            for static in [True, False]:
                fa_html_gen = generate_html.ST_HTML_Generator(
                    [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info"], "The statistics for FA", para_dict], static=True)
                fa_html_gen.generate_st_html()

            print("Call FA-module done!")


def bam_module(margs):
    """
    Run the BAM filetype module.
    """
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
        exit_code = lrst.callBAMModule(input_para, bam_output)
        if exit_code == 0:
            from src import plot_for_BAM
            plot_for_BAM.bam_plot(bam_output, para_dict)

            for static in [True, False]:
                bam_html_gen = generate_html.ST_HTML_Generator(
                    [["basic_st","map_st", "err_st", "read_length_st", "read_length_hist", "base_st", "basic_info", "base_quality"], "The statistics for BAM", para_dict], static=static)
                bam_html_gen.generate_st_html()


            print("Call BAM-module done!")

def seqtxt_module(margs):
    """
    Run the sequencing_summary.txt filetype module.
    """
    para_dict = {}
    errorStr = lrst_global.originalError
    errorStr += get_common_param(margs, para_dict)

    if not errorStr == lrst_global.originalError:
        print(errorStr)
        print("#############################################\n")
        # parser.print_help()
        parser.parse_args(['seqtxt', '--help'])
        sys.exit(1004)
    else:
        logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]))
        para_dict["out_prefix"] += "seqtxt_";
        input_para = lrst.Input_Para()
        input_para.threads = para_dict["threads"]
        input_para.rdm_seed = para_dict["random_seed"]
        input_para.downsample_percentage = para_dict["downsample_percentage"]

        # TODO: 'seq' flag is 1 by default, which in the argument description states that only the summary text will be
        #  output (no HTML).
        input_para.other_flags = margs.seq  # Default = 1
        input_para.other_flags = (input_para.other_flags << 4)
        input_para.other_flags += (1 if para_dict["detail"]>0 else 0)
        input_para.other_flags = (input_para.other_flags << 4)
        input_para.other_flags += int(margs.sum_type)

        input_para.output_folder = str(para_dict["output_folder"])
        input_para.out_prefix = str(para_dict["out_prefix"])

        for _ipf in para_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        seqtxt_output = lrst.Output_SeqTxt()
        exit_code = lrst.callSeqTxtModule(input_para, seqtxt_output)
        if exit_code == 0:
            print("Generating output files...")
            from src import plot_for_SeqTxt
            plot_for_SeqTxt.plot(seqtxt_output, para_dict)
            for static in [True, False]:
                if margs.seq==0:
                    f5_html_gen = generate_html.ST_HTML_Generator([["basic_st", "read_length_st","read_length_hist","base_st","basic_info"], "QC for sequencing_summary.txt", para_dict], static=static);
                else:
                    f5_html_gen = generate_html.ST_HTML_Generator(
                       [["basic_st", "read_length_st","read_length_hist", "basic_info"], "QC for sequencing_summary.txt", para_dict], static=static)
                f5_html_gen.generate_st_html()
            print("Done.")

def fast5_module(margs):
    """
    Run the FAST5 filetype module.
    """
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
        para_dict["out_prefix"] += "f5_"
        input_para = lrst.Input_Para()
        input_para.threads = para_dict["threads"]
        input_para.rdm_seed = para_dict["random_seed"]
        input_para.downsample_percentage = para_dict["downsample_percentage"]
        input_para.output_folder = str(para_dict["output_folder"])
        input_para.out_prefix = str(para_dict["out_prefix"])

        for _ipf in para_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fast5_output = lrst.Output_FAST5()
        exit_code = lrst.callFAST5Module(input_para, fast5_output)
        if exit_code == 0:
            print("Generating output files...")
            from src import plot_for_FAST5
            plot_for_FAST5.plot(fast5_output, para_dict)
            for static in [True, False]:
                fast5_html_obj = generate_html.ST_HTML_Generator(
                    [["basic_st", "read_length_st","read_length_hist", "base_st", "basic_info", "base_quality", "read_avg_base_quality"], "The statistics for FQ", para_dict], static=static)
                fast5_html_obj.generate_st_html()
            print("Done.")


# =====
# Set up the argument parser
parser = argparse.ArgumentParser(description="QC tools for long-read sequencing data",
                                 epilog="Example with single inputs:\n"
                                        "\tpython %(prog)s bam -i path/to/input.bam -o /output_directory/\n"
                                        "\nExample with multiple inputs:\n"
                                        "\tpython %(prog)s bam -I \"path/to/input1.bam, path/to/input2.bam\" -o /output_directory/\n"
                                        "\tpython %(prog)s bam -P \"path/to/*.bam\" -o /output_directory/\n",
                                 formatter_class=RawTextHelpFormatter)

# The subparser will determine our filetype-specific modules
subparsers = parser.add_subparsers()

# The parent parser contains parameters common to all modules (input/output file, etc.)
parent_parser = argparse.ArgumentParser(add_help=False)

common_grp_param = parent_parser.add_argument_group(
    "Common parameters for %(prog)s")
input_files_group = common_grp_param.add_mutually_exclusive_group()
input_files_group.add_argument(
    "-i", "--input", type=str, default=None, help="Single input filepath")
input_files_group.add_argument(
    "-I", "--inputs", type=str, default=None,
    help="Multiple comma-separated input filepaths",)
input_files_group.add_argument(
    "-P", "--inputPattern", type=str, default=None,
    help="Use pattern matching (*) to specify multiple input files. Enclose the pattern in double quotes.")
input_files_group.add_argument("-p", "--downsample_percentage", type=float, default=1.0,
                               help="The percentage of downsampling for quick run. Default: 1.0 without downsampling")

common_grp_param.add_argument(
    "-g", "--log", type=str, default="log_output.log", help="Log file")
common_grp_param.add_argument("-G", "--Log_level", type=int, default=lrst_global.LOG_ERROR,
                              help="Level for logging: ALL(0) < DEBUG(1) < INFO(2) < WARN(3) < ERROR(4) < FATAL(5) < OFF(6). Default: 4 (ERROR)")

common_grp_param.add_argument("-o", "--outputfolder", type=str,
                              default="output_" + lrst_global.prg_name, help="The output folder.")
common_grp_param.add_argument("-t", "--thread", type=int,
                              default=1, help="The number of threads used. Default: 1.")
common_grp_param.add_argument("-Q", "--outprefix", type=str,
                              default="st_", help="The prefix of output. Default: `st_`.")
common_grp_param.add_argument(
    "-s", "--seed", type=int, default=1, help="The number for random seed. Default: 1.")
common_grp_param.add_argument("-d", "--detail", type=int, default=0,
                              help="Will output detail in files? Default: 0(no).")

# FASTA inputs
fa_parsers = subparsers.add_parser('fa',
                                   parents=[parent_parser],
                                   help="FASTA file input",
                                   description="For example:\n"
                                               "python %(prog)s -i input.fasta -o /output_directory/",
                                   formatter_class=RawTextHelpFormatter)
fa_parsers.set_defaults(func=fa_module)

# FASTQ inputs
fq_parsers = subparsers.add_parser('fq',
                                   parents=[parent_parser],
                                   help="FASTQ file input",
                                   description="For example:\n"
                                               "python %(prog)s -i input.fastq -o /output_directory/",
                                   formatter_class=RawTextHelpFormatter)
fq_parsers.add_argument("-u", "--udqual", type=int, default=-1,
                        help="User defined quality offset for bases in fq. Default: -1.")
fq_parsers.set_defaults(func=fq_module)

# FAST5 inputs
fast5_parser = subparsers.add_parser('f5',
                                     parents=[parent_parser],
                                     help="FAST5 file input",
                                     description="For example:\n"
                                                 "python %(prog)s -i input.fast5 -o /output_directory/",
                                     formatter_class=RawTextHelpFormatter)
fast5_parser.set_defaults(func=fast5_module)

# sequencing_summary.txt inputs
seqtxt_parsers = subparsers.add_parser('seqtxt',
                                       parents=[parent_parser],
                                       help="sequencing_summary.txt input",
                                       description="For example:\n"
                                                   "python %(prog)s -i sequencing_summary.txt -o /output_directory/",
                                       formatter_class=RawTextHelpFormatter)
seqtxt_parsers.add_argument("-S", "--seq", type=int, default=1,
                            help="sequencing_summary.txt only? Default: 1(yes).")
seqtxt_parsers.add_argument("-m", "--sum_type", type=int, default=1, choices=[1, 2, 3],
                            help="Different fields in sequencing_summary.txt. Default: 1.")

seqtxt_parsers.set_defaults(func=seqtxt_module)

# BAM inputs
bam_parsers = subparsers.add_parser('bam',
                                    parents=[parent_parser],
                                    help="BAM file input",
                                    description="For example:\n"
                                                "python %(prog)s -i input.bam -o /output_directory/",
                                    formatter_class=RawTextHelpFormatter)
bam_parsers.set_defaults(func=bam_module)


# =====
def main():
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

                # Run the filetype-specific module
                args.func(args)
