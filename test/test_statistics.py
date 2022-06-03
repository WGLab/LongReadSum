"""
test_statistics.py:
Test expected values for output statistics using sample input files (FASTA, FAST5, FASTQ, BAM).
"""
import os
import lrst


def test_fasta_n50():
    """
    Test that the calculated N50 value is correct for FASTA files.
    """
    # Define expected values
    expected_n50_read_length = 23758872

    # Set parameters
    input_para = lrst.Input_Para()
    input_para.threads = 1
    input_para.rdm_seed = 1
    input_para.downsample_percentage = 1.0
    input_para.other_flags = 0
    output_folder = os.path.abspath(str("output/"))
    input_para.output_folder = output_folder
    input_file = os.path.abspath(str("TestInputs/all_chr.hap1.cns.fa"))
    input_para.add_input_file(input_file)
    input_para.out_prefix = str("fa_")

    # Run the FASTA statistics module
    fa_output = lrst.Output_FA()
    lrst.generate_statistic_from_fa(input_para, fa_output)

    # Ensure the N50 value is correct
    n50_read_length = fa_output.long_read_info.n50_read_length
    assert n50_read_length == expected_n50_read_length


def test_two_fasta_n50():
    """
    Test that the calculated N50 value is correct when the input is two FASTA files.
    """
    # Define expected values
    expected_n50_read_length = 23544988

    # Set parameters
    input_para = lrst.Input_Para()
    input_para.threads = 1
    input_para.rdm_seed = 1
    input_para.downsample_percentage = 1.0
    input_para.other_flags = 0
    output_folder = os.path.abspath(str("output/"))
    input_para.output_folder = output_folder
    input_file1 = os.path.abspath(str("TestInputs/all_chr.hap1.cns.fa"))
    input_para.add_input_file(input_file1)
    input_file2 = os.path.abspath(str("TestInputs/all_chr.hap2.cns.fa"))
    input_para.add_input_file(input_file2)
    input_para.out_prefix = str("fa_")

    # Run the FASTA statistics module
    fa_output = lrst.Output_FA()
    lrst.generate_statistic_from_fa(input_para, fa_output)

    # Ensure the N50 value is correct
    n50_read_length = fa_output.long_read_info.n50_read_length
    assert n50_read_length == expected_n50_read_length
