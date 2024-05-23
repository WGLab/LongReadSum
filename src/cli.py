#!/usr/bin/env python
# CLI.py: Parse arguments and run the filetype-specific module.

import os
import sys
import logging
from glob import glob
import argparse
from argparse import RawTextHelpFormatter


# Print the package name
if __name__ == 'src.cli':
    # logging.debug("Running locally.")
    from lib import lrst  # For debugging
    from src import generate_html
    from src.plot_utils import *
    from src.pod5_module import generate_pod5_qc
else:
    # logging.debug("Running from installed package.")
    import lrst
    import generate_html
    from plot_utils import *
    from pod5_module import generate_pod5_qc

prg_name = "LongReadSum"


# Return the log level type from the code (Default: ERROR)
def get_log_level(level_code):
    switch = {
        1: logging.DEBUG,
        2: logging.INFO,
        3: logging.WARN,
        4: logging.ERROR,
        5: logging.FATAL,
    }
    return switch.get(level_code)


def get_common_param(margs):
    """
    Format a dict object param_dict to contain our input files/output folder parameters
    Output a string containing all parse argument errors.
    Inputs:
    margs: Command line input arguments (dict).
    """
    param_dict = {}
    param_dict["prg_name"] = prg_name
    parsing_error_msg = ""

    if (margs.input == None or margs.input == "") and (margs.inputs == None or margs.inputs == "") and (
            margs.pattern == None or margs.pattern == ""):
        parsing_error_msg += "No input file(s) are provided. \n"
    else:
        # Group parameters into an array
        param_dict["input_files"] = []
        if not (margs.input == None or margs.input == ""):
            input_filepath = margs.input.name
            param_dict["input_files"].append(input_filepath)
        if not (margs.inputs == None or margs.inputs == ""):
            file_list = [file_str.strip() for file_str in margs.inputs.split(',')]
            param_dict["input_files"].extend(file_list)
        if not (margs.pattern == None or margs.pattern == ""):
            pat_split = margs.pattern.split("*")
            param_dict["input_files"].extend(
                glob(os.path.join("*".join(pat_split[:-1]), "*" + pat_split[-1])))

        # Number of reads to sample
        read_count = int(margs.read_count[0])
        param_dict["read_count"] = read_count
        logging.info("Number of reads to sample: %d", read_count)

        if len(param_dict["input_files"]) == 0:
            parsing_error_msg += "No input file(s) can be found. \n"
        else:
            for input_filepath in param_dict["input_files"]:
                if not os.path.isfile(input_filepath):
                    parsing_error_msg += "Cannot find the input file: " + input_filepath + "\n"

    if (margs.outputfolder is None or margs.outputfolder == ""):
        parsing_error_msg += "No output file is provided. \n"
    else:
        output_dir = margs.outputfolder
        param_dict["output_folder"] = output_dir
        try:
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)

        except OSError:
            parsing_error_msg += "Cannot create folder for " + param_dict["output_folder"] + " \n"

    param_dict["out_prefix"] = margs.outprefix

    # Set up logging to file and stdout
    if margs.log is None or margs.log == "":
        parsing_error_msg += "No log file is provided. \n"

        # Set up logging to stdout
        logging.basicConfig(stream=sys.stdout,
                            level=get_log_level(margs.log_level),
                            format="%(asctime)s [%(levelname)s] %(message)s")
    else:
        logging.basicConfig(level=get_log_level(margs.log_level),
                            format="%(asctime)s [%(levelname)s] %(message)s",
                            handlers=[
                                logging.FileHandler(margs.log),
                                logging.StreamHandler(sys.stdout)
                                ]
                            )

        logging.info("Log file is %s", margs.log)
        param_dict["log_file"] = margs.log
        param_dict["log_level"] = margs.log_level

    param_dict["downsample_percentage"] = margs.downsample_percentage

    param_dict["threads"] = margs.threads

    param_dict["random_seed"] = margs.seed

    param_dict["detail"] = margs.detail

    # Plot style parameters
    param_dict["fontsize"] = margs.fontsize
    param_dict["markersize"] = margs.markersize

    # Reset the param_dict if there are parsing errors
    if parsing_error_msg != "":
        param_dict = {}
        logging.error(parsing_error_msg)

    return param_dict


def fq_module(margs):
    # Run the FASTQ filetype module.

    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['fq', '--help'])
        sys.exit(0)

    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "fastq"

        # Import the SWIG Python wrapper for our C++ module
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]

        input_para.other_flags = 0
        input_para.user_defined_fastq_base_qual_offset = margs.udqual;

        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])

        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fq_output = lrst.Output_FQ()
        exit_code = lrst.callFASTQModule(input_para, fq_output)
        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(fq_output, param_dict, 'FASTQ')
            fq_html_gen = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist", "base_counts", "base_quality",
                  "read_avg_base_quality"], "FASTQ QC", param_dict], plot_filepaths, static=False)
            fq_html_gen.generate_st_html()

            logging.info("Done. Output files are in %s", param_dict["output_folder"])
        else:
            logging.error("QC did not generate.")


def fa_module(margs):
    # Run the FASTA filetype module.

    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['fa', '--help'])
        sys.exit(0)
        
    else:
        # If there are no parse errors, run the filetype-specific module
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "fasta"
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]
        input_para.other_flags = 0
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])

        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fa_output = lrst.Output_FA()
        exit_code = lrst.callFASTAModule(input_para, fa_output)
        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(fa_output, param_dict, 'FASTA')
            fa_html_gen = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist", "base_counts"], "FASTA QC",
                 param_dict], plot_filepaths, static=True)
            fa_html_gen.generate_st_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")

def bam_module(margs):
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['bam', '--help'])
        sys.exit(0)

    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "bam";
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]
        input_para.other_flags = (1 if param_dict["detail"] > 0 else 0);
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])
        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        bam_output = lrst.Output_BAM()
        exit_code = lrst.callBAMModule(input_para, bam_output)
        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(bam_output, param_dict, 'BAM')

            # TODO: Add read average base quality plot (not currently generated by bam_plot.plot)
            bam_html_gen = generate_html.ST_HTML_Generator(
                [["basic_st", "read_alignments_bar", "base_alignments_bar", "read_length_bar", "read_length_hist", "base_counts", "basic_info",
                  "base_quality"], "BAM QC", param_dict], plot_filepaths, static=False)
            bam_html_gen.generate_st_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")


def rrms_module(margs):
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['rrms', '--help'])
        sys.exit(0)
    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]
        input_para.other_flags = (1 if param_dict["detail"] > 0 else 0);
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])
        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        # Set the RRMS input CSV file
        input_para.rrms_csv = margs.csv
        logging.info("RRMS CSV file is " + input_para.rrms_csv)

        # Get the output prefix
        output_prefix = param_dict["out_prefix"]

        # Run QC for both accepted and rejected reads
        rrms_filter = [True, False]
        for filter_type in rrms_filter:

            # Set the RRMS filter type
            input_para.rrms_filter = filter_type

            # Set the output prefix
            param_dict["out_prefix"] = output_prefix + "rrms_" + ("accepted" if filter_type else "rejected")

            # Run the QC module
            logging.info("Running QC for " + ("accepted" if filter_type else "rejected") + " reads...")
            bam_output = lrst.Output_BAM()
            exit_code = lrst.callBAMModule(input_para, bam_output)
            if exit_code == 0:
                logging.info("QC generated.")
                logging.info("Generating HTML report...")
                plot_filepaths = plot(bam_output, param_dict, 'BAM')

                # Generate the HTML report
                bam_html_gen = generate_html.ST_HTML_Generator(
                    [["basic_st", "read_alignments_bar", "base_alignments_bar", "read_length_bar", "read_length_hist", "base_counts", "basic_info",
                    "base_quality"], "BAM QC", param_dict], plot_filepaths, static=False)
                bam_html_gen.generate_st_html()
                logging.info("Done. Output files are in %s", param_dict["output_folder"])

            else:
                logging.error("QC did not generate.")
        

def seqtxt_module(margs):
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['seqtxt', '--help'])
        sys.exit(0)
        
    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "seqtxt"
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]
        input_para.other_flags = margs.seq  # Default = 1
        input_para.other_flags = (input_para.other_flags << 4)
        input_para.other_flags += (1 if param_dict["detail"] > 0 else 0)
        input_para.other_flags = (input_para.other_flags << 4)
        input_para.other_flags += int(margs.sum_type)

        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])

        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        seqtxt_output = lrst.Output_SeqTxt()
        exit_code = lrst.callSeqTxtModule(input_para, seqtxt_output)
        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(seqtxt_output, param_dict, 'SeqTxt')

            if margs.seq == 0:
                seqtxt_html_gen = generate_html.ST_HTML_Generator(
                    [["basic_st", "read_length_bar", "read_length_hist", "base_counts", "base_quality", "basic_info"],
                     "sequencing_summary.txt QC", param_dict], plot_filepaths, static=False)
            else:
                seqtxt_html_gen = generate_html.ST_HTML_Generator(
                    [["basic_st", "read_length_bar", "read_length_hist", "basic_info"], "sequencing_summary.txt QC",
                     param_dict], plot_filepaths, static=False)
            seqtxt_html_gen.generate_st_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])
        else:
            logging.error("QC did not generate.")


def fast5_module(margs):
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['f5', '--help'])
        sys.exit(0)

    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "FAST5"
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])
        input_para.other_flags = 0  # 0 for normal QC, 1 for signal statistics output

        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fast5_output = lrst.Output_FAST5()
        exit_code = lrst.callFAST5Module(input_para, fast5_output)
        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(fast5_output, param_dict, 'FAST5')
            fast5_html_obj = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist", "base_counts", "basic_info", "base_quality",
                  "read_avg_base_quality"], "FAST5 QC", param_dict], plot_filepaths, static=False)
            fast5_html_obj.generate_st_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")


def fast5_signal_module(margs):
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['f5s', '--help'])
        sys.exit(0)

    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "fast5_signal"
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.rdm_seed = param_dict["random_seed"]
        input_para.downsample_percentage = param_dict["downsample_percentage"]
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])
        input_para.other_flags = 1  # 0 for normal QC, 1 for signal statistics output

        # Get the read ID list if specified
        read_ids = margs.read_ids
        if read_ids != "" and read_ids is not None:
            input_para.read_ids = read_ids
            #print("Read ID list is " + str(read_ids))

        for _ipf in param_dict["input_files"]:
            input_para.add_input_file(str(_ipf))

        fast5_output = lrst.Output_FAST5()
        exit_code = lrst.callFAST5Module(input_para, fast5_output)

        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(fast5_output, param_dict, 'FAST5s')
            fast5_html_obj = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist", "base_counts", "basic_info", "base_quality",
                  "read_avg_base_quality", "ont_signal"], "FAST5 QC", param_dict], plot_filepaths, static=False)
            fast5_html_obj.generate_st_html(signal_plots=True)
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")


def pod5_module(margs):
    """POD5 file input module."""
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['pod5', '--help'])
        sys.exit(0)

    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "POD5"
        input_para = {}
        input_para['threads'] = param_dict["threads"]
        input_para['rdm_seed'] = param_dict["random_seed"]
        input_para['downsample_percentage'] = param_dict["downsample_percentage"]
        input_para['output_folder'] = str(param_dict["output_folder"])
        input_para['out_prefix'] = str(param_dict["out_prefix"])
        input_para['other_flags'] = 0  # 0 for normal QC, 1 for signal statistics output
        input_para['input_files'] = []
        for input_file in param_dict["input_files"]:
            input_para['input_files'].append(str(input_file))

        # Get the read ID list if specified
        read_ids = margs.read_ids
        if read_ids != "" and read_ids is not None:
            input_para['read_ids'] = read_ids
        else:
            input_para['read_ids'] = ""

        read_signal_dict = generate_pod5_qc(input_para)
        if read_signal_dict is not None:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot_pod5(read_signal_dict, param_dict)
            # plot_filepaths = plot(read_signal_dict, param_dict, 'POD5')
            webpage_title = "POD5 QC"
            fast5_html_obj = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist", "base_counts", "basic_info", "base_quality",
                  "read_avg_base_quality", "ont_signal"], webpage_title, param_dict], plot_filepaths, static=False)
            fast5_html_obj.generate_st_html(signal_plots=True)
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")


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

# File input parameter
input_files_group = common_grp_param.add_mutually_exclusive_group()
input_files_group.add_argument(
    "-i", "--input", type=argparse.FileType('r'), default=None, help="Single input filepath")

input_files_group.add_argument(
    "-I", "--inputs", type=str, default=None,
    help="Multiple comma-separated input filepaths", )

input_files_group.add_argument(
    "-P", "--pattern", type=str, default=None,
    help="Use pattern matching (*) to specify multiple input files. Enclose the pattern in double quotes.")

# Plot style parameters
common_grp_param.add_argument("--fontsize", type=int, default=14,
                              help="Font size for plots. Default: 14")

common_grp_param.add_argument("--markersize", type=int, default=10,
                              help="Marker size for plots. Default: 10")

# Number of reads to sample
common_grp_param.add_argument(
    "-R", "--read-count", nargs="+", type=int, default=8,
    help="Set the number of reads to randomly sample from the file. Default: 3.")

# Misc. parameters
input_files_group.add_argument("-p", "--downsample_percentage", type=float, default=1.0,
                               help="The percentage of downsampling for quick run. Default: 1.0 without downsampling")

common_grp_param.add_argument(
    "-g", "--log", type=str, default="log_output.log", help="Log file")
common_grp_param.add_argument("-G", "--log_level", type=int, default=2,
                              help="Logging level. 1: DEBUG, 2: INFO, 3: WARNING, 4: ERROR, 5: CRITICAL. Default: 2.")

common_grp_param.add_argument("-o", "--outputfolder", type=str,
                              default="output_" + prg_name, help="The output folder.")
common_grp_param.add_argument("-t", "--threads", type=int,
                              default=1, help="The number of threads used. Default: 1.")
common_grp_param.add_argument("-Q", "--outprefix", type=str,
                              default="QC_", help="The prefix of output. Default: `QC_`.")
common_grp_param.add_argument(
    "-s", "--seed", type=int, default=1, help="The number for random seed. Default: 1.")
common_grp_param.add_argument("-d", "--detail", type=int, default=0,
                              help="Will output detail in files? Default: 0(no).")

# FASTA input file
fa_parser = subparsers.add_parser('fa',
                                   parents=[parent_parser],
                                   help="FASTA file input",
                                   description="For example:\n"
                                               "python %(prog)s -i input.fasta -o /output_directory/",
                                   formatter_class=RawTextHelpFormatter)
fa_parser.set_defaults(func=fa_module)

# FASTQ input file
fq_parser = subparsers.add_parser('fq',
                                   parents=[parent_parser],
                                   help="FASTQ file input",
                                   description="For example:\n"
                                               "python %(prog)s -i input.fastq -o /output_directory/",
                                   formatter_class=RawTextHelpFormatter)
fq_parser.add_argument("-u", "--udqual", type=int, default=-1,
                        help="User defined quality offset for bases in fq. Default: -1.")
fq_parser.set_defaults(func=fq_module)

# FAST5 input file
fast5_parser = subparsers.add_parser('f5',
                                     parents=[parent_parser],
                                     help="FAST5 file input",
                                     description="For example:\n"
                                                 "python %(prog)s -i input.fast5 -o /output_directory/",
                                     formatter_class=RawTextHelpFormatter)
fast5_parser.set_defaults(func=fast5_module)

# FAST5 input file in signal statistics mode
fast5_signal_parser = subparsers.add_parser('f5s',
                                            parents=[parent_parser],
                                            help="FAST5 file input with signal statistics output",
                                            description="For example:\n"
                                                        "python %(prog)s -R 5 10 -i input.fast5 -o /output_directory/",
                                            formatter_class=RawTextHelpFormatter)
fast5_signal_parser.set_defaults(func=fast5_signal_module)

# Add an argument for specifying the read names to extract
fast5_signal_parser.add_argument("-r", "--read_ids", type=str, default=None,
                                 help="A comma-separated list of read IDs to extract from the file.")

# POD5 input file
pod5_parser = subparsers.add_parser('pod5',
                                    parents=[parent_parser],
                                    help="POD5 file input",
                                    description="For example:\n"
                                                "python %(prog)s -i input.pod5 -o /output_directory/",
                                    formatter_class=RawTextHelpFormatter)
pod5_parser.set_defaults(func=pod5_module)

# Add an argument for specifying the read names to extract
pod5_parser.add_argument("-r", "--read_ids", type=str, default=None,
                            help="A comma-separated list of read IDs to extract from the file.")

# Sequencing summary text file input
seqtxt_parser = subparsers.add_parser('seqtxt',
                                       parents=[parent_parser],
                                       help="sequencing_summary.txt input",
                                       description="For example:\n"
                                                   "python %(prog)s -i sequencing_summary.txt -o /output_directory/",
                                       formatter_class=RawTextHelpFormatter)
seqtxt_parser.add_argument("-S", "--seq", type=int, default=1,
                            help="sequencing_summary.txt only? Default: 1(yes).")
seqtxt_parser.add_argument("-m", "--sum_type", type=int, default=1, choices=[1, 2, 3],
                            help="Different fields in sequencing_summary.txt. Default: 1.")

seqtxt_parser.set_defaults(func=seqtxt_module)

# BAM file input
bam_parser = subparsers.add_parser('bam',
                                    parents=[parent_parser],
                                    help="BAM file input",
                                    description="For example:\n"
                                                "python %(prog)s -i input.bam -o /output_directory/",
                                    formatter_class=RawTextHelpFormatter)
bam_parser.set_defaults(func=bam_module)

# RRMS BAM file input (Splits accepted and rejected reads)
rrms_bam_parser = subparsers.add_parser('rrms',
                                         parents=[parent_parser],
                                         help="RRMS BAM file input",
                                         description="For example:\n"
                                                     "python %(prog)s -i input.bam -c input.csv -o /output_directory/",
                                         formatter_class=RawTextHelpFormatter)

# Add required input CSV file for RRMS
rrms_help_str = "CSV file containing read IDs to extract from the BAM file.\n" \
                "The CSV file should contain a read id column with the header 'read_id' as well as a column " \
                "containing the accepted/rejected status of the read with the header 'decision'.\n" \
                "Accepted reads should have a value of 'stop_receiving' in the 'decision' column, while rejected reads " \
                "should have a value of 'unblock'."
                
rrms_bam_parser.add_argument("-c", "--csv", type=str, required=True,
                                help=rrms_help_str)

rrms_bam_parser.set_defaults(func=rrms_module)


def main():
    if sys.version_info[0] < 2:
        logging.info("%s could not be run with lower version than python 2.7.", prg_name)
    else:
        if sys.version_info[1] < 6:
            logging.info("%s could be run with python 2.7.", prg_name)
        else:
            if len(sys.argv) < 2:
                parser.print_help()
            else:
                args = parser.parse_args()

                # Run the filetype-specific module
                args.func(args)
