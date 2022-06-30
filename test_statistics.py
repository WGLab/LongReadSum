"""
test_statistics.py:
Test expected values for output statistics using sample input files (FASTA, FAST5, FASTQ, BAM).
"""
import os
import lrst
import pytest


# Fixtures for module outputs
@pytest.fixture(scope='class')
def fasta_output():
    """
    Run the FASTA module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("fa_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/fasta_trim1.fa") # Local path
    else:
        input_file = os.path.abspath(str("SampleData/fasta_trim1.fa"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_FA()
    exit_code = lrst.callFASTAModule(default_parameters, output)

    yield [exit_code, output]

@pytest.fixture(scope='class')
def multiple_fasta_output():
    """
    Run the FASTA module for multiple inputs.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("faX2_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        # Local paths
        input_file1 = os.path.join(local_dir, "SampleData/fasta_trim1.fa")
        input_file2 = os.path.join(local_dir, "SampleData/fasta_trim2.fa")
    else:
        # Remote paths
        input_file1 = os.path.abspath(str("SampleData/fasta_trim1.fa"))
        input_file2 = os.path.abspath(str("SampleData/fasta_trim2.fa"))

    # Add input files
    default_parameters.add_input_file(input_file1)
    default_parameters.add_input_file(input_file2)

    # Run the FASTA statistics module
    output = lrst.Output_FA()
    exit_code = lrst.callFASTAModule(default_parameters, output)

    yield [exit_code, output]


@pytest.fixture(scope='class')
def fastq_output():
    """
    Run the FASTQ module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("fq_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/guppy.fastq") # Local path
    else:
        input_file = os.path.abspath(str("SampleData/guppy.fastq"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_FQ()
    exit_code = lrst.callFASTQModule(default_parameters, output)

    yield [exit_code, output]


@pytest.fixture(scope='class')
def bam_output():
    """
    Run the BAM module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("bam_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/guppy.bam") # Local path
    else:
        input_file = os.path.abspath(str("SampleData/guppy.bam"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output]


@pytest.fixture(scope='class')
def unmapped_bam_output():
    """
    Run the BAM module on unmapped inputs.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("bam_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/pacbio_unmapped_trim.bam") # Local path
    else:
        input_file = os.path.abspath(str("SampleData/pacbio_unmapped_trim.bam"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output]


# Filetype-specific unit test classes
class TestFASTA:
    """
    Unit tests for FASTA inputs.
    """
    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, fasta_output):
        exit_code = fasta_output[0]
        assert exit_code == 0

    # Unit tests
    @pytest.mark.dependency(depends=["TestFASTA::test_success"])
    def test_base_count(self, fasta_output):
        output_statistics = fasta_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 126055194

    @pytest.mark.dependency(depends=["TestFASTA::test_success"])
    def test_read_count(self, fasta_output):
        output_statistics = fasta_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 100

    @pytest.mark.dependency(depends=["TestFASTA::test_success"])
    def test_n50(self, fasta_output):
        output_statistics = fasta_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 9726200


class TestMultipleFASTA:
    """
    Unit tests for multiple FASTA inputs
    """
    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, multiple_fasta_output):
        exit_code = multiple_fasta_output[0]
        assert exit_code == 0

    # Unit tests
    @pytest.mark.dependency(depends=["TestMultipleFASTA::test_success"])
    def test_base_count(self, multiple_fasta_output):
        output_statistics = multiple_fasta_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 186794132

    @pytest.mark.dependency(depends=["TestMultipleFASTA::test_success"])
    def test_read_count(self, multiple_fasta_output):
        output_statistics = multiple_fasta_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 201

    @pytest.mark.dependency(depends=["TestMultipleFASTA::test_success"])
    def test_n50(self, multiple_fasta_output):
        output_statistics = multiple_fasta_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 15495669


class TestFASTQ:
    """
    Unit tests for FASTQ inputs.
    """
    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, fastq_output):
        exit_code = fastq_output[0]
        assert exit_code == 0

    # Unit tests
    @pytest.mark.dependency(depends=["TestFASTQ::test_success"])
    def test_base_count(self, fastq_output):
        output_statistics = fastq_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 340189

    @pytest.mark.dependency(depends=["TestFASTQ::test_success"])
    def test_read_count(self, fastq_output):
        output_statistics = fastq_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 50

    @pytest.mark.dependency(depends=["TestFASTQ::test_success"])
    def test_n50(self, fastq_output):
        output_statistics = fastq_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 8731


class TestBAM:
    """
    Unit tests for BAM inputs.
    """
    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, bam_output):
        exit_code = bam_output[0]
        assert exit_code == 0

    # Unit tests
    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_base_count(self, bam_output):
        output_statistics = bam_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 340189

    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_read_count(self, bam_output):
        output_statistics = bam_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 50

    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_n50(self, bam_output):
        output_statistics = bam_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 7415


class TestUnmappedBAM:
    """
    Unit tests for unmapped BAM inputs.
    """
    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, unmapped_bam_output):
        exit_code = unmapped_bam_output[0]
        assert exit_code == 0

    # Unit tests
    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_base_count(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 1297818

    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_read_count(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 95

    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_n50(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 21391
