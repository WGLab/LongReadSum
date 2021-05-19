#!/usr/bin/env python

import os,sys,string,math
from glob import glob
import argparse, logging
from argparse import RawTextHelpFormatter

import faulthandler; faulthandler.enable()

import lrst_global

def get_log_level(log_l):
   if log_l==lrst_global.LOG_ALL: return logging.ALL;
   elif log_l==lrst_global.LOG_DEBUG: return logging.DEBUG;
   elif log_l==lrst_global.LOG_INFO: return logging.INFO;
   elif log_l==lrst_global.LOG_WARN: return logging.WARN;
   elif log_l==lrst_global.LOG_ERROR: return logging.ERROR;
   elif log_l==lrst_global.LOG_FATAL: return logging.FATAL;
   elif log_l==lrst_global.LOG_OFF: return logging.OFF;
   else: return logging.ERROR;

def get_common_param(margs, para_dict):
   this_error_str = "";

   if (margs.input==None or margs.input=="") and (margs.inputs==None or margs.inputs=="") and (margs.inputPattern==None or margs.inputPattern==""):
      this_error_str += "No input file(s) are provided. \n"
   else:
      para_dict["input_files"] = []
      if not (margs.input==None or margs.input==""):
         para_dict["input_files"].append(margs.input)
      if not (margs.inputs==None or margs.inputs==""):
         para_dict["input_files"].extend(margs.inputs.split(','))
      if not (margs.inputPattern==None or margs.inputPattern==""):
         pat_split = margs.inputPattern.split("*")
         para_dict["input_files"].extend( glob(os.path.join("*".join(pat_split[:-1]), "*"+pat_split[-1])) );
      if len( para_dict["input_files"] )==0: 
         this_error_str += "No input file(s) can be found. \n"

   if (margs.outputfolder==None or margs.outputfolder==""):
      this_error_str += "No output file is provided. \n"
   else: 
      para_dict["output_folder"] = margs.outputfolder;
      try:
         if not os.path.isdir( para_dict["output_folder"] ):
            os.makedirs(para_dict["output_folder"]);
      except OSError as e:
         this_error_str += "Cannot create folder for "+para_dict["output_folder"]+" \n"
   para_dict["out_prefix"] = margs.outprefix;

   if (margs.log==None or margs.log==""):
      this_error_str += "No log file is provided. \n"
   else:
      para_dict["log_file"] = margs.log
      para_dict["log_level"] = margs.Log_level
      logging.basicConfig(filename=margs.log, level=get_log_level(margs.Log_level),
                          filemode='w', format="%(levelname)s: %(message)s")

   para_dict["downsample_percentage"] = margs.downsample_percentage;

   para_dict["threads"] = margs.thread;

   para_dict["random_seed"] = margs.seed;

   return this_error_str;


def fq_module(margs):
   para_dict = {};
   errorStr = lrst_global.originalError; 
   errorStr += get_common_param(margs, para_dict);
   
   if not errorStr == lrst_global.originalError:
      print(errorStr)  
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['fq', '--help'])
      sys.exit(1001)
   else:
      pass

def fa_module(margs):
   para_dict = {};
   errorStr = lrst_global.originalError;
   errorStr += get_common_param(margs, para_dict);

   if not errorStr == lrst_global.originalError:
      print(errorStr)
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['fa', '--help'])
      sys.exit(1002)
   else:
      pass

def bam_module(margs):
   para_dict = {};
   errorStr = lrst_global.originalError;
   errorStr += get_common_param(margs, para_dict);

   if not errorStr == lrst_global.originalError:
      print(errorStr)
      print("#############################################\n")
      #parser.print_help()
      parser.parse_args(['bam', '--help'])
      sys.exit(1003)
   else:
      logging.info('Input file(s) are ' + ';'.join(para_dict["input_files"]) )
      import lrst;
      input_para = lrst.Input_Para();
      input_para.threads = para_dict["threads"];
      input_para.rdm_seed = para_dict["random_seed"];
      input_para.downsample_percentage = para_dict["downsample_percentage"];

      input_para.other_flags = 0 ;

      input_para.output_folder = str(para_dict["output_folder"]) ;
      input_para.out_prefix = str(para_dict["out_prefix"]);
      #input_para.set_output_folder( para_dict["output_folder"]) ;
      #input_para.set_out_prefix( para_dict["out_prefix"]);

      for _ipf in para_dict["input_files"]:
         input_para.add_input_file( str(_ipf) ); 

      bam_output = lrst.Output_BAM();
      lrst.generate_statistic_from_bam( input_para, bam_output );
      import plot_for_BAM;
      plot_for_BAM.bam_plot(bam_output)
      import generate_html
      bam_html_gen = generate_html.ST_HTML_Generator([["read_length_distr", "map_st"], "The statistics for BAM", para_dict ]);
      bam_html_gen.generate_st_html();
      print("Call BAM-module done!")

def f5_module(margs):
   para_dict = {};
   errorStr = lrst_global.originalError;
   errorStr += get_common_param(margs, para_dict);

   if not errorStr == lrst_global.originalError:
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
input_files_group.add_argument("-P", "--inputPattern", type=str, default=None, help="The pattern of input files with *.");
input_files_group.add_argument("-p", "--downsample_percentage", type=float, default=1.0, help="The percentage of downsampling for quick run. Default: 1.0 without downsampling");

common_grp_param.add_argument("-g", "--log", type=str, default="log_output.log", help="Log file")
common_grp_param.add_argument("-G", "--Log_level", type=int, default=lrst_global.LOG_ERROR, help="Level for logging: ALL(0) < DEBUG(1) < INFO(2) < WARN(3) < ERROR(4) < FATAL(5) < OFF(6). Default: 4 (ERROR)")

common_grp_param.add_argument("-o", "--outputfolder", type=str, default="output_"+lrst_global.prg_name, help="The output folder.")
common_grp_param.add_argument("-t", "--thread", type=int, default=1, help="The number of threads used. Default: 1.")
common_grp_param.add_argument("-Q", "--outprefix", type=str, default="st_", help="The prefix of output. Default: `st_`.")
common_grp_param.add_argument("-s", "--seed", type=int, default=1, help="The number for random seed. Default: 1.")


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
    print(lrst_global.prg_name+" could not be run with lower version than python 2.7.")
else:
    if sys.version_info[1] < 6:
        print(lrst_global.prg_name+" could be run with python 2.7.")
    else:
        if len(sys.argv) < 2:
            parser.print_help()
        else:
            args = parser.parse_args()
            args.func(args)



