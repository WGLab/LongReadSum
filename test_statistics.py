"""
test_statistics.py:
Test expected values for output statistics using sample input files (FASTA, FAST5, FASTQ, BAM).
"""
import os
import lrst
import pytest


def test_fasta_n50():
    """
    Test that the calculated N50 value is correct for FASTA files.
    """

    # Set parameters
    input_para = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    input_para.output_folder = output_folder
    input_file = os.path.abspath(str("SampleData/fasta_trim1.fa"))  # Remote path
    # input_file = str("/home/perdomoj/github/LongReadSum/SampleData trimmed/SampleData/fasta_trim1.fa") # Local path
    input_para.add_input_file(input_file)
    input_para.out_prefix = str("fa_")

    # Run the FASTA statistics module
    fa_output = lrst.Output_FA()
    exit_code = lrst.callFASTAModule(input_para, fa_output)

    # Ensure the N50 value is correct
    n50_read_length = fa_output.long_read_info.n50_read_length
    assert (exit_code == 0) and (n50_read_length == 9726200)
