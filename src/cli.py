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
        parsing_error_msg += "No input file(s) were provided.\n"
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

        # read_count = int(margs.read_count)
        # param_dict["read_count"] = read_count
        
        if len(param_dict["input_files"]) == 0:
            parsing_error_msg += "No input file(s) were provided.\n"
        else:
            for input_filepath in param_dict["input_files"]:
                if not os.path.isfile(input_filepath):
                    parsing_error_msg += "Cannot find the input file: " + input_filepath + "\n"

    if (margs.outputfolder is None or margs.outputfolder == ""):
        parsing_error_msg += "No output directory was provided.\n"
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
        parsing_error_msg += "No log file was provided.\n"

        # Set up logging to stdout
        logging.basicConfig(stream=sys.stdout,
                            level=get_log_level(margs.log_level),
                            format="%(asctime)s %(message)s")
                            # format="%(asctime)s [%(levelname)s] %(message)s")
    else:
        logging.basicConfig(level=get_log_level(margs.log_level),
                            format="%(asctime)s %(message)s",
                            handlers=[
                                logging.FileHandler(margs.log),
                                logging.StreamHandler(sys.stdout)
                                ]
                            )

        logging.info("Log file is %s", margs.log)
        param_dict["log_file"] = margs.log
        param_dict["log_level"] = margs.log_level

    param_dict["threads"] = margs.threads

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

        input_para.other_flags = 0
        input_para.user_defined_fastq_base_qual_offset = margs.udqual

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
                [["basic_st", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts", "base_quality",
                  "read_avg_base_quality"], "FASTQ QC", param_dict], plot_filepaths, static=False)
            fq_html_gen.generate_html()

            logging.info("Done. Output files are in %s", param_dict["output_folder"])
        else:
            logging.error("QC did not generate.")


def fa_module(margs):
    """FASTA file input module."""

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
                [["basic_st", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts"], "FASTA QC",
                 param_dict], plot_filepaths, static=True)
            fa_html_gen.generate_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")

def bam_module(margs):
    """BAM file input module."""
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
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])

        # Set the reference genome file and base modification threshold
        ref_genome = margs.ref if margs.ref != "" or margs.ref is not None else ""
        param_dict["ref"] = input_para.ref_genome = ref_genome

        # Set the base modification flag, and filtering threshold
        # param_dict["mod"] = input_para.mod_analysis = margs.mod
        if margs.mod:
            param_dict["mod"] = input_para.mod_analysis = True
        else:
            param_dict["mod"] = input_para.mod_analysis = False
            
        mod_prob = margs.modprob
        param_dict["modprob"] = mod_prob
        input_para.base_mod_threshold = mod_prob
        logging.info("Base modification threshold is set to " + str(input_para.base_mod_threshold))

        # Set the gene BED file for RNA-seq transcript analysis
        input_para.gene_bed = margs.genebed if margs.genebed != "" or margs.genebed is not None else ""
        param_dict["genebed"] = input_para.gene_bed

        # Set the minimum coverage and sample size for TIN calculation
        input_para.tin_sample_size = margs.sample_size
        input_para.tin_min_coverage = margs.min_coverage

        # Add the input files to the input parameter
        for input_file in param_dict["input_files"]:
            input_para.add_input_file(str(input_file))

        bam_output = lrst.Output_BAM()
        exit_code = lrst.callBAMModule(input_para, bam_output)
        if exit_code == 0:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")
            plot_filepaths = plot(bam_output, param_dict, 'BAM')

            # Set the list of QC information to display
            qc_info_list = ["basic_st", "read_alignments_bar", "base_alignments_bar", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts", "base_quality", "read_avg_base_quality"]

            # If base modifications were found, add the base modification plots
            # after the first table
            if bam_output.sample_modified_base_count > 0:
                # logging.info("Base modifications found. Adding base modification plots to the HTML report.")
                qc_info_list.insert(1, "read_length_mod_rates")  # Read length modification rates
                qc_info_list.insert(1, "base_mods")

            # If gene BED file was provided, add the TIN plots
            if input_para.gene_bed != "":
                qc_info_list.insert(1, "tin")

            # Generate the HTML report
            bam_html_gen = generate_html.ST_HTML_Generator(
                [qc_info_list, "BAM QC", param_dict], plot_filepaths, static=False)
            bam_html_gen.generate_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")


def rrms_module(margs):
    """RRMS BAM file input module."""
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['rrms', '--help'])
        sys.exit(0)
    else:
        logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
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
            param_dict["mod"] = input_para.mod_analysis = False  # Disable base modification analysis for RRMS (use BAM module for this)

            # Run the QC module
            logging.info("Running QC for " + ("accepted" if filter_type else "rejected") + " reads...")
            bam_output = lrst.Output_BAM()
            exit_code = lrst.callBAMModule(input_para, bam_output)
            if exit_code == 0:
                logging.info("QC generated.")
                logging.info("Generating HTML report...")
                plot_filepaths = plot(bam_output, param_dict, 'BAM')

                # Set the list of QC information to display
                qc_info_list = ["basic_st", "read_alignments_bar", "base_alignments_bar", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts", "base_quality"]

                # If base modifications were found, add the base modification
                # plots
                if bam_output.sample_modified_base_count > 0:
                    logging.info("Base modifications found. Adding base modification plots to the HTML report.")
                    qc_info_list.insert(1, "read_length_mod_rates")
                    qc_info_list.insert(1, "base_mods")

                # Generate the HTML report
                bam_html_gen = generate_html.ST_HTML_Generator(
                    [qc_info_list, "BAM QC", param_dict], plot_filepaths, static=False)
                bam_html_gen.generate_html()
                logging.info("Done. Output files are in %s", param_dict["output_folder"])

            else:
                logging.error("QC did not generate.")
        

def seqtxt_module(margs):
    """Basecall summary text file input module."""
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

            report_title = "Basecall Summary QC"
            seqtxt_html_gen = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist"],
                    report_title, param_dict], plot_filepaths, static=False)
                
            seqtxt_html_gen.generate_html()
            logging.info("Done. Output files are in %s", param_dict["output_folder"])
        else:
            logging.error("QC did not generate.")


def fast5_module(margs):
    """FAST5 file input module."""
    # Get the filetype-specific parameters
    param_dict = get_common_param(margs)
    if param_dict == {}:
        parser.parse_args(['f5', '--help'])
        sys.exit(0)

    else:
        # logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "FAST5"
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
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
                [["basic_st", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts", "base_quality"], 
                 "FAST5 QC", param_dict], plot_filepaths, static=False)
            fast5_html_obj.generate_html()
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
        # logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "fast5_signal"
        input_para = lrst.Input_Para()
        input_para.threads = param_dict["threads"]
        input_para.output_folder = str(param_dict["output_folder"])
        input_para.out_prefix = str(param_dict["out_prefix"])
        input_para.other_flags = 1  # 0 for normal QC, 1 for signal statistics output

        # Get the read count if specified
        read_count = int(margs.read_count)
        param_dict["read_count"] = read_count

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
                [["basic_st", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts", "ont_signal"], "FAST5 QC", param_dict], plot_filepaths, static=False)
            fast5_html_obj.generate_html(signal_plots=True)
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
        # logging.info('Input file(s) are:\n%s', '\n'.join(param_dict["input_files"]))
        param_dict["out_prefix"] += "POD5"
        input_para = {}
        input_para['threads'] = param_dict["threads"]
        input_para['output_folder'] = str(param_dict["output_folder"])
        input_para['out_prefix'] = str(param_dict["out_prefix"])
        input_para['other_flags'] = 0  # 0 for normal QC, 1 for signal statistics output
        input_para['input_files'] = []
        for input_file in param_dict["input_files"]:
            input_para['input_files'].append(str(input_file))

        # Get the read count if specified
        read_count = int(margs.read_count)
        param_dict["read_count"] = read_count

        # Get the read ID list if specified
        read_ids = margs.read_ids
        if read_ids != "" and read_ids is not None:
            input_para['read_ids'] = read_ids
        else:
            input_para['read_ids'] = ""

        # Get the basecalled BAM file if specified, and run the BAM module
        basecall_data = False
        bam_output = None
        basecalls = margs.basecalls
        if basecalls != "" and basecalls is not None:
            basecalls_input = lrst.Input_Para()
            basecalls_input.threads = param_dict["threads"]
            basecalls_input.output_folder = str(param_dict["output_folder"])
            basecalls_input.out_prefix = str(param_dict["out_prefix"])
            basecalls_input.add_input_file(basecalls)
            bam_output = lrst.Output_BAM()
            exit_code = lrst.callBAMModule(basecalls_input, bam_output)
            if exit_code == 0:
                basecall_data = True
                logging.info("Basecalled BAM QC generated.")

        read_signal_dict = generate_pod5_qc(input_para)
        if read_signal_dict is not None:
            logging.info("QC generated.")
            logging.info("Generating HTML report...")

            if basecall_data:
                plot_filepaths = plot_pod5(read_signal_dict, param_dict, bam_output)
            else:
                plot_filepaths = plot(read_signal_dict, param_dict, None)
                
            # plot_filepaths = plot(read_signal_dict, param_dict, 'POD5')
            webpage_title = "POD5 QC"
            fast5_html_obj = generate_html.ST_HTML_Generator(
                [["basic_st", "read_length_bar", "read_length_hist", "gc_content_hist", "base_counts", "ont_signal"], webpage_title, param_dict], plot_filepaths, static=False)
            fast5_html_obj.generate_html(signal_plots=True)
            logging.info("Done. Output files are in %s", param_dict["output_folder"])

        else:
            logging.error("QC did not generate.")
            

def set_file_parser_defaults(file_parser):
    """Create a parser with default arguments for a specific filetype."""
    file_parser.add_argument("-i", "--input", type=argparse.FileType('r'), default=None,
                        help="Single input filepath")
    file_parser.add_argument("-I", "--inputs", type=str, default=None,
                        help="Multiple comma-separated input filepaths")
    file_parser.add_argument("-P", "--pattern", type=str, default=None,
                        help="Use pattern matching (*) to specify multiple input files. Enclose the pattern in double quotes.")
    file_parser.add_argument("-g", "--log", type=str, default="log_output.log",
                        help="Log file")
    file_parser.add_argument("-G", "--log-level", type=int, default=2,
                        help="Logging level. 1: DEBUG, 2: INFO, 3: WARNING, 4: ERROR, 5: CRITICAL. Default: 2.")
    file_parser.add_argument("-o", "--outputfolder", type=str, default="output_" + prg_name,
                        help="The output folder.")
    file_parser.add_argument("-t", "--threads", type=int, default=1,
                        help="The number of threads used. Default: 1.")
    file_parser.add_argument("-Q", "--outprefix", type=str, default="QC_",
                        help="The prefix for output filenames. Default: `QC_`.")


# Set up the argument parser
parser = argparse.ArgumentParser(description="QC summary tool for long-read sequencing data.",
                                 formatter_class=RawTextHelpFormatter)

# The subparser will determine our filetype-specific modules
subparsers = parser.add_subparsers(title="File types", dest="filetype")

# The parent parser contains parameters common to all modules (input/output file, etc.)
parent_parser = argparse.ArgumentParser(add_help=False)

# FASTA input file
fa_parser = subparsers.add_parser('fa',
                                   parents=[parent_parser],
                                   help="FASTA file input",
                                   description="For example:\n"
                                               "python %(prog)s -i input.fasta -o /output_directory/",
                                   formatter_class=RawTextHelpFormatter)

set_file_parser_defaults(fa_parser)
fa_parser.set_defaults(func=fa_module)

# FASTQ input file
fq_parser = subparsers.add_parser('fq',
                                   parents=[parent_parser],
                                   help="FASTQ file input",
                                   description="For example:\n"
                                               "python %(prog)s -i input.fastq -o /output_directory/",
                                   formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(fq_parser)
fq_parser.add_argument("-u", "--udqual", type=int, default=33,
                        help="User defined quality offset for bases in fq. Default: 33.")
fq_parser.set_defaults(func=fq_module)

# FAST5 input file
fast5_parser = subparsers.add_parser('f5',
                                     parents=[parent_parser],
                                     help="FAST5 file input",
                                     description="For example:\n"
                                                 "python %(prog)s -i input.fast5 -o /output_directory/",
                                     formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(fast5_parser)
fast5_parser.set_defaults(func=fast5_module)

# FAST5 input file in signal statistics mode
fast5_signal_parser = subparsers.add_parser('f5s',
                                            parents=[parent_parser],
                                            help="FAST5 file input with signal statistics output",
                                            description="For example:\n"
                                                        "python %(prog)s -R 5 10 -i input.fast5 -o /output_directory/",
                                            formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(fast5_signal_parser)
fast5_signal_parser.set_defaults(func=fast5_signal_module)

# Add an argument for specifying the read names to extract
fast5_signal_parser.add_argument("-r", "--read_ids", type=str, default=None,
                                 help="A comma-separated list of read IDs to extract from the file.")

# Add an argument for specifying the maximum number of reads to extract
fast5_signal_parser.add_argument("-R", "--read-count", type=int, default=3,
                                    help="Set the number of reads to randomly sample from the file. Default: 3.")

# POD5 input file
pod5_parser = subparsers.add_parser('pod5',
                                    parents=[parent_parser],
                                    help="POD5 file input",
                                    description="For example:\n"
                                                "python %(prog)s -i input.pod5 -o /output_directory/",
                                    formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(pod5_parser)
pod5_parser.set_defaults(func=pod5_module)

# Add an argument for specifying the basecalled BAM file
pod5_parser.add_argument("-b", "--basecalls", type=str, default=None,
                            help="The basecalled BAM file to use for signal extraction.")

# Add an argument for specifying the read names to extract
pod5_parser.add_argument("-r", "--read_ids", type=str, default=None,
                            help="A comma-separated list of read IDs to extract from the file.")

# Add an argument for specifying the maximum number of reads to extract
pod5_parser.add_argument("-R", "--read-count", type=int, default=3,
                            help="Set the number of reads to randomly sample from the file. Default: 3.")

# Sequencing summary text file input
seqtxt_parser = subparsers.add_parser('seqtxt',
                                       parents=[parent_parser],
                                       help="Basecall summary (sequencing_summary.txt) input",
                                       description="For example:\n"
                                                   "python %(prog)s -i sequencing_summary.txt -o /output_directory/",
                                       formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(seqtxt_parser)
seqtxt_parser.set_defaults(func=seqtxt_module)

# BAM file input
bam_parser = subparsers.add_parser('bam',
                                    parents=[parent_parser],
                                    help="BAM file input",
                                    description="For example:\n"
                                                "python %(prog)s -i input.bam -o /output_directory/",
                                    formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(bam_parser)

bam_parser.add_argument("--mod", action="store_true",
                        help="Run base modification analysis on the BAM file. Default: False.")

# Add argument for gene BED file required for RNA-seq transcript analysis (TIN, etc.)
bam_parser.add_argument("--genebed", type=str, default="",
                        help="Gene BED12 file required for calculating TIN scores from RNA-seq BAM files. Default: None.")

bam_parser.add_argument("--modprob", type=float, default=0.5,
                        help="Base modification filtering threshold. Above/below this value, the base is considered modified/unmodified. Default: 0.5.")

bam_parser.add_argument("--ref", type=str, default="",
                        help="The reference genome FASTA file to use for identifying CpG sites.")

# Add TIN sample size argument
bam_parser.add_argument("--sample-size", type=int, default=100,
                        help="Sample size for TIN calculation. Default: 100.")

# Add TIN minimum coverage argument
bam_parser.add_argument("--min-coverage", type=int, default=10,
                        help="Minimum coverage for TIN calculation. Default: 10.")

bam_parser.set_defaults(func=bam_module)

# RRMS BAM file input (Splits accepted and rejected reads)
rrms_bam_parser = subparsers.add_parser('rrms',
                                         parents=[parent_parser],
                                         help="RRMS BAM file input",
                                         description="For example:\n"
                                                     "python %(prog)s -i input.bam -c input.csv -o /output_directory/",
                                         formatter_class=RawTextHelpFormatter)
set_file_parser_defaults(rrms_bam_parser)

# Add required input CSV file for RRMS
rrms_help_str = "CSV file containing read IDs to extract from the BAM file. See README for formatting specifics."

# "The CSV file should contain a read id column with the header 'read_id' as well as a column " \
# "containing the accepted/rejected status of the read with the header 'decision'." \
# "Accepted reads should have a value of 'stop_receiving' in the 'decision' column, while rejected reads " \
# "should have a value of 'unblock'."
                
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
