#!/usr/bin/env python

import os,sys,string,math

import argparse, logging
from argparse import RawTextHelpFormatter


prg_name = "LongReadDS"
originalError = '!!!Error: !!!!!! \n'

LOG_ALL = 0;
LOG_DEBUG = 1;
LOG_INFO = 2;
LOG_WARN = 3;
LOG_ERROR = 4;
LOG_FATAL = 5;
LOG_OFF = 6;

def get_log_level(log_l):
   if log_l==LOG_ALL: return logging.ALL;
   elif log_l==LOG_DEBUG: return logging.DEBUG;
   elif log_l==LOG_INFO: return logging.INFO;
   elif log_l==LOG_WARN: return logging.WARN;
   elif log_l==LOG_ERROR: return logging.ERROR;
   elif log_l==LOG_FATAL: return logging.FATAL;
   elif log_l==LOG_OFF: return logging.OFF;
   else: return logging.ERROR;

def get_common_param(margs, para_dict):
   this_error_str = "";

   if (margs.input==None or margs.input=="") and (margs.inputs==None or margs.inputs==""):
      this_error_str += "No input file(s) are provided. \n"
   else:
      para_dict["input_files"] = []
      if not (margs.input==None or margs.input==""):
         para_dict["input_files"].append(margs.input)
      if not (margs.inputs==None or margs.inputs==""):
         para_dict["input_files"].extend(margs.inputs.split(','))

   if (margs.outputfolder==None or margs.outputfolder==""):
      this_error_str += "No output file is provided. \n"
   else: 
      para_dict["output_folder"] = margs.outputfolder;
      try:
         os.makedirs(para_dict["output_folder"]);
      except OSError as e:
         this_error_str += "Cannot create folder for "+para_dict["output_folder"]+" \n"

   if (margs.log==None or margs.log==""):
      this_error_str += "No log file is provided. \n"
   else:
      para_dict["log_file"] = margs.log
      para_dict["log_level"] = margs.Log_level
      logging.basicConfig(filename=margs.log, level=get_log_level(margs.Log_level),
                          filemode='w', format="%(levelname)s: %(message)s")
      logging.info('Input file(s) are ' + para_dict["input_files"])

   return this_error_str;


def fq_module(margs):
   para_dict = {};
   errorStr = originalError; 
   errorStr += get_common_param(margs, para_dict);
   
   if not errorStr == originalError:
      print(errorStr)  
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['fq', '--help'])
      sys.exit(1001)
   else:
      pass

def fa_module(margs):
   para_dict = {};
   errorStr = originalError;
   errorStr += get_common_param(margs, para_dict);

   if not errorStr == originalError:
      print(errorStr)
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['fa', '--help'])
      sys.exit(1002)
   else:
      pass

def bam_module(margs):
   para_dict = {};
   errorStr = originalError;
   errorStr += get_common_param(margs, para_dict);

   if not errorStr == originalError:
      print(errorStr)
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['bam', '--help'])
      sys.exit(1003)
   else:
      pass

def f5_module(margs):
   para_dict = {};
   errorStr = originalError;
   errorStr += get_common_param(margs, para_dict);

   if not errorStr == originalError:
      print(errorStr)
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['f5', '--help'])
      sys.exit(1004)
   else:
      pass


parser = argparse.ArgumentParser(description="Data analysis tools for long-read sequencing data", 
                                 epilog="For example, \n \
                                         \tpython %(prog)s fq: with fq or fastq input\n\
                                         \tpython %(prog)s fa: with fa or fasta input\n\
                                         \tpython %(prog)s bam: with bam input\n\
                                         \tpython %(prog)s f5: with fast5 input\n\
                                         ",
                                 formatter_class=RawTextHelpFormatter);


##############################################################################
#
# Different functional modules starting from here.
#
##############################################################################


subparsers = parser.add_subparsers()
parent_parser = argparse.ArgumentParser(add_help=False)

common_grp_param = parent_parser.add_argument_group("Common parameters for %(prog)s")
#common_grp_param.add_argument("-h", "--help", default="", help="Show this help documents");
input_files_group = common_grp_param.add_mutually_exclusive_group();
input_files_group.add_argument("-i", "--input", type=str, default=None, help="The input file for the analysis");
input_files_group.add_argument("-I", "--inputs", type=str, default=None, help="The input files for the analysis");

common_grp_param.add_argument("-g", "--log", type=str, default="", help="Log file")
common_grp_param.add_argument("-G", "--Log_level", type=int, default=LOG_ERROR, help="Level for logging: ALL(0) < DEBUG(1) < INFO(2) < WARN(3) < ERROR(4) < FATAL(5) < OFF(6). Default: 4 (ERROR)")

common_grp_param.add_argument("-o", "--outputfolder", type=str, default="output_"+prg_name, help="The output folder.")
common_grp_param.add_argument("-t", "--thread", type=int, default=1, help="The number of threads used. Default: 1.")

fq_parsers = subparsers.add_parser('fq', 
                                    parents=[parent_parser], 
                                    help="Show data analysis for fq files", 
                                    description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                    formatter_class=RawTextHelpFormatter)
fq_parsers.set_defaults(func=fq_module);


fa_parsers = subparsers.add_parser('fa',
                                    parents=[parent_parser],
                                    help="Show data analysis for fa files",
                                    description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                    formatter_class=RawTextHelpFormatter)
fa_parsers.set_defaults(func=fa_module);

bam_parsers = subparsers.add_parser('bam',
                                    parents=[parent_parser],
                                    help="Show data analysis for bam files",
                                    description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                    formatter_class=RawTextHelpFormatter)
bam_parsers.set_defaults(func=bam_module);

f5_parsers = subparsers.add_parser('f5',
                                    parents=[parent_parser],
                                    help="Show data analysis for f5 files",
                                    description="For example, \n \
                                                 \tpython %(prog)s  \n \
                                                ",
                                    formatter_class=RawTextHelpFormatter)
f5_parsers.set_defaults(func=f5_module);

if sys.version_info[0] < 2:
    print(prg_name+" could not be run with lower version than python 2.7.")
else:
    if sys.version_info[1] < 6:
        print(prg_name+" could be run with python 2.7.")
    else:
        if len(sys.argv) < 2:
            parser.print_help()
        else:
            args = parser.parse_args()
            args.func(args)



