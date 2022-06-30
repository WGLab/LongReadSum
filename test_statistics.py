"""
test_statistics.py:
Test expected values for output statistics using sample input files (FASTA, FAST5, FASTQ, BAM).
"""
import os
import lrst
import pytest

@pytest.fixture(scope='class')
def fasta_output():
    """
    Run the FASTA module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    input_file = os.path.abspath(str("SampleData/fasta_trim1.fa"))  # Remote path
    # input_file = str("/home/perdomoj/github/LongReadSum/SampleData trimmed/SampleData/fasta_trim1.fa") # Local path    default_parameters.add_input_file(input_file)
    default_parameters.add_input_file(input_file)
    default_parameters.out_prefix = str("fa_")

    # Run the FASTA statistics module
    output = lrst.Output_FA()
    exit_code = lrst.callFASTAModule(default_parameters, output)

    yield [exit_code, output]

class TestFASTA:
    def test_module_load(self, fasta_output):
        exit_code, output_statistics = fasta_output

        # Ensure the N50 value is correct
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert (exit_code == 0) and (n50_read_length == 9726200)

