"""
test_statistics.py:
Test expected values for output statistics using sample input files (FASTA, FAST5, FASTQ, BAM).
"""
import os
import pytest


# Add the parent directory to the path
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib import lrst


# FASTA tests
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
        input_file = os.path.join(local_dir, "SampleData/fasta_trim1.fa")  # Local path
    else:
        input_file = os.path.abspath(str("SampleData/fasta_trim1.fa"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_FA()
    exit_code = lrst.callFASTAModule(default_parameters, output)

    yield [exit_code, output]


class TestFASTA:
    """
    Tests for FASTA inputs.
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, fasta_output):
        exit_code = fasta_output[0]
        assert exit_code == 0

    # Tests
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


@pytest.fixture(scope='class')
def multiple_fasta_output():
    """
    Run the FASTA module for multiple inputs.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("fa_multi_")

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


class TestMultipleFASTA:
    """
    Tests for multiple FASTA inputs
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, multiple_fasta_output):
        exit_code = multiple_fasta_output[0]
        assert exit_code == 0

    # Tests
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


# FASTQ tests
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
        input_file = os.path.join(local_dir, "SampleData/guppy.fastq")  # Local path
    else:
        input_file = os.path.abspath(str("SampleData/guppy.fastq"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTQ statistics module
    output = lrst.Output_FQ()
    exit_code = lrst.callFASTQModule(default_parameters, output)

    yield [exit_code, output]


class TestFASTQ:
    """
    Tests for FASTQ inputs.
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, fastq_output):
        exit_code = fastq_output[0]
        assert exit_code == 0

    # Tests
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
        assert n50_read_length == 8733


# FAST5 tests
@pytest.fixture(scope='class')
def fast5_output():
    """
    Run the FAST5 module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("f5_")

    # Check if running remotely
    file_dir = ''
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        file_dir = os.path.join(local_dir, "SampleData/")  # Local path
    else:
        file_dir = os.path.abspath(str("SampleData/"))  # Remote path

    # # Get all FAST5 files
    # input_files = []
    # input_files.extend(glob.glob(file_dir + '*.fast5'))

    # Get the first 5 FAST5 files (Github Actions throws a seg.fault error if all 50 are used)
    file_names = \
        ["kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch101_read101553_strand.fast5",
         "kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch101_read110620_strand.fast5",
         "kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch101_read152950_strand.fast5",
         "kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch100_read6082_strand.fast5",
         "kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch101_read77163_strand.fast5"]

    # Add input files
    for file_name in file_names:
        input_file = os.path.join(file_dir, file_name)
        default_parameters.add_input_file(input_file)

    # Run the FAST5 statistics module
    output = lrst.Output_FAST5()
    exit_code = lrst.callFAST5Module(default_parameters, output)

    yield [exit_code, output]


class TestFAST5:
    """
    Tests for FAST5 inputs.
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, fast5_output):
        exit_code = fast5_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestFAST5::test_success"])
    def test_base_count(self, fast5_output):
        output_statistics = fast5_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 28581

    @pytest.mark.dependency(depends=["TestFAST5::test_success"])
    def test_read_count(self, fast5_output):
        output_statistics = fast5_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 5

    @pytest.mark.dependency(depends=["TestFAST5::test_success"])
    def test_n50(self, fast5_output):
        output_statistics = fast5_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 7050


# FAST5 signal tests
@pytest.fixture(scope='class')
def fast5s_output():
    """
    Run the FAST5 signal QC module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("f5s_")
    default_parameters.other_flags = 1  # 0 for normal QC, 1 for signal statistics output

    # Check if running remotely
    file_dir = ''
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        file_dir = os.path.join(local_dir, "SampleData/")  # Local path
    else:
        file_dir = os.path.abspath(str("SampleData/"))  # Remote path

    # Get the first 2 FAST5 files
    file_names = [
        "kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch100_read6082_strand.fast5",
        "kelvin_20160810_FN_MN17519_sequencing_run_160810_na12878_PCRSssI_92870_ch101_read110620_strand.fast5"
    ]

    # Add input files
    for file_name in file_names:
        input_file = os.path.join(file_dir, file_name)
        default_parameters.add_input_file(input_file)

    # Run the FAST5 statistics module
    output = lrst.Output_FAST5()
    exit_code = lrst.callFAST5Module(default_parameters, output)

    yield [exit_code, output]


class TestFAST5Signal:
    """
    Tests for FAST5 inputs with signal QC output.
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, fast5s_output):
        exit_code = fast5s_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestFAST5Signal::test_success"])
    def test_base_count(self, fast5s_output):
        output_statistics = fast5s_output[1]
        base_count = output_statistics.getTotalBaseCount()
        assert base_count == 9720

    @pytest.mark.dependency(depends=["TestFAST5Signal::test_success"])
    def test_read_count(self, fast5s_output):
        output_statistics = fast5s_output[1]
        read_count = output_statistics.getReadCount()
        assert read_count == 2

    @pytest.mark.dependency(depends=["TestFAST5Signal::test_success"])
    def test_window_length(self, fast5s_output):
        output_statistics = fast5s_output[1]
        first_read_data = output_statistics.getNthReadBaseSignals(0)
        window_length = len(first_read_data[0])
        assert window_length == 5


# BAM tests
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
        input_file = os.path.join(local_dir, "SampleData/guppy.bam")  # Local path
    else:
        input_file = os.path.abspath(str("SampleData/guppy.bam"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output]


class TestBAM:
    """
    Tests for BAM inputs.
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, bam_output):
        exit_code = bam_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_base_count(self, bam_output):
        output_statistics = bam_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 340189

    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_mapped_base_count(self, bam_output):
        output_statistics = bam_output[1]
        mapped_base_count = output_statistics.num_matched_bases
        assert mapped_base_count == 304323

    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_read_count(self, bam_output):
        output_statistics = bam_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 50

    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_mapped_read_count(self, bam_output):
        output_statistics = bam_output[1]
        mapped_read_count = output_statistics.num_primary_alignment
        assert mapped_read_count == 50

    @pytest.mark.dependency(depends=["TestBAM::test_success"])
    def test_n50(self, bam_output):
        output_statistics = bam_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 8733


@pytest.fixture(scope='class')
def unmapped_bam_output():
    """Run the BAM module on unmapped inputs."""
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("ubam_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/pacbio_unmapped_trim.bam")  # Local path
    else:
        input_file = os.path.abspath(str("SampleData/pacbio_unmapped_trim.bam"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the BAM statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output]


class TestUnmappedBAM:
    """Tests for unmapped BAM inputs."""

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, unmapped_bam_output):
        exit_code = unmapped_bam_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_base_count(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        base_count = output_statistics.long_read_info.total_num_bases
        assert base_count == 1297818

    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_mapped_base_count(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        mapped_base_count = output_statistics.num_matched_bases
        assert mapped_base_count == 0

    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_read_count(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        read_count = output_statistics.long_read_info.total_num_reads
        assert read_count == 95

    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_mapped_read_count(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        mapped_read_count = output_statistics.num_primary_alignment
        assert mapped_read_count == 0

    @pytest.mark.dependency(depends=["TestUnmappedBAM::test_success"])
    def test_n50(self, unmapped_bam_output):
        output_statistics = unmapped_bam_output[1]
        n50_read_length = output_statistics.long_read_info.n50_read_length
        assert n50_read_length == 22029


@pytest.fixture(scope='class')
def forward_base_mod_output():
    """Run the BAM module on a read aligned to the forward strand with base modifications."""
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("fwdmod_")
    default_parameters.mod_analysis = True
    default_parameters.base_mod_threshold = -1.0

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/forward_mod.bam")  # Local path
        ref_file = os.path.join(local_dir, "SampleData/chr11.fa")
    else:
        input_file = os.path.abspath(str("SampleData/forward_mod.bam"))  # Remote path
        ref_file = os.path.abspath(str("SampleData/chr11.fa"))

    # Add input files
    default_parameters.add_input_file(input_file)
    default_parameters.ref_genome = ref_file

    # Run the BAM statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output]


class TestForwardBaseModBAM:
    """Tests for BAM inputs with base modifications on the forward strand."""

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, forward_base_mod_output):
        exit_code = forward_base_mod_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestForwardBaseModBAM::test_success"])
    def test_modified_base_count(self, forward_base_mod_output):
        output_statistics = forward_base_mod_output[1]
        modified_base_count = output_statistics.sample_modified_base_count
        assert modified_base_count == 695

    @pytest.mark.dependency(depends=["TestForwardBaseModBAM::test_success"])
    def test_forward_modified_base_count(self, forward_base_mod_output):
        output_statistics = forward_base_mod_output[1]
        forward_modified_base_count = output_statistics.sample_modified_base_count_forward
        assert forward_modified_base_count == 695

    @pytest.mark.dependency(depends=["TestForwardBaseModBAM::test_success"])
    def test_reverse_modified_base_count(self, forward_base_mod_output):
        output_statistics = forward_base_mod_output[1]
        reverse_modified_base_count = output_statistics.sample_modified_base_count_reverse
        assert reverse_modified_base_count == 0

    @pytest.mark.dependency(depends=["TestForwardBaseModBAM::test_success"])
    def test_forward_cpg_modified_base_count(self, forward_base_mod_output):
        output_statistics = forward_base_mod_output[1]
        forward_cpg_modified_base_count = output_statistics.sample_cpg_forward_count
        assert forward_cpg_modified_base_count == 621

    @pytest.mark.dependency(depends=["TestForwardBaseModBAM::test_success"])
    def test_reverse_cpg_modified_base_count(self, forward_base_mod_output):
        output_statistics = forward_base_mod_output[1]
        reverse_cpg_modified_base_count = output_statistics.sample_cpg_reverse_count
        assert reverse_cpg_modified_base_count == 0


@pytest.fixture(scope='class')
def reverse_base_mod_output():
    """Run the BAM module on a read aligned to the reverse strand with base modifications."""
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("revmod_")
    default_parameters.mod_analysis = True
    default_parameters.base_mod_threshold = -1.0

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/reverse_mod.bam")
        ref_file = os.path.join(local_dir, "SampleData/chr11.fa")
    else:
        input_file = os.path.abspath(str("SampleData/reverse_mod.bam"))
        ref_file = os.path.abspath(str("SampleData/chr11.fa"))

    # Add input files
    default_parameters.add_input_file(input_file)
    default_parameters.ref_genome = ref_file

    # Run the BAM statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output]

class TestReverseBaseModBam:
    """Tests for BAM inputs with base modifications on the reverse strand."""

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, reverse_base_mod_output):
        exit_code = reverse_base_mod_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestReverseBaseModBam::test_success"])
    def test_modified_base_count(self, reverse_base_mod_output):
        output_statistics = reverse_base_mod_output[1]
        modified_base_count = output_statistics.sample_modified_base_count
        assert modified_base_count == 1094

    @pytest.mark.dependency(depends=["TestReverseBaseModBam::test_success"])
    def test_forward_modified_base_count(self, reverse_base_mod_output):
        output_statistics = reverse_base_mod_output[1]
        forward_modified_base_count = output_statistics.sample_modified_base_count_forward
        assert forward_modified_base_count == 0

    @pytest.mark.dependency(depends=["TestReverseBaseModBam::test_success"])
    def test_reverse_modified_base_count(self, reverse_base_mod_output):
        output_statistics = reverse_base_mod_output[1]
        reverse_modified_base_count = output_statistics.sample_modified_base_count_reverse
        assert reverse_modified_base_count == 1094

    @pytest.mark.dependency(depends=["TestReverseBaseModBam::test_success"])
    def test_forward_cpg_modified_base_count(self, reverse_base_mod_output):
        output_statistics = reverse_base_mod_output[1]
        forward_cpg_modified_base_count = output_statistics.sample_cpg_forward_count
        assert forward_cpg_modified_base_count == 0

    @pytest.mark.dependency(depends=["TestReverseBaseModBam::test_success"])
    def test_reverse_cpg_modified_base_count(self, reverse_base_mod_output):
        output_statistics = reverse_base_mod_output[1]
        reverse_cpg_modified_base_count = output_statistics.sample_cpg_reverse_count
        assert reverse_cpg_modified_base_count == 525


# sequencing_summary.txt tests
@pytest.fixture(scope='class')
def seqtxt_output():
    """
    Run the sequencing_summary.txt module.
    """
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("seqtxt_")

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/sequencing_summary.txt")  # Local path
    else:
        input_file = os.path.abspath(str("SampleData/sequencing_summary.txt"))  # Remote path

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the FASTA statistics module
    output = lrst.Output_SeqTxt()
    exit_code = lrst.callSeqTxtModule(default_parameters, output)

    yield [exit_code, output]


class TestSeqTxt:
    """
    Tests for sequencing_summary.txt inputs.
    """

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, seqtxt_output):
        exit_code = seqtxt_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestSeqTxt::test_success"])
    def test_base_count(self, seqtxt_output):
        output_statistics = seqtxt_output[1]
        base_count = output_statistics.all_long_read_info.long_read_info.total_num_bases
        assert base_count == 340189

    @pytest.mark.dependency(depends=["TestSeqTxt::test_success"])
    def test_passed_base_count(self, seqtxt_output):
        output_statistics = seqtxt_output[1]
        passed_base_count = output_statistics.passed_long_read_info.long_read_info.total_num_bases
        assert passed_base_count == 117697

    @pytest.mark.dependency(depends=["TestSeqTxt::test_success"])
    def test_read_count(self, seqtxt_output):
        output_statistics = seqtxt_output[1]
        read_count = output_statistics.all_long_read_info.long_read_info.total_num_reads
        assert read_count == 50

    @pytest.mark.dependency(depends=["TestSeqTxt::test_success"])
    def test_passed_read_count(self, seqtxt_output):
        output_statistics = seqtxt_output[1]
        passed_read_count = output_statistics.passed_long_read_info.long_read_info.total_num_reads
        assert passed_read_count == 17

    @pytest.mark.dependency(depends=["TestSeqTxt::test_success"])
    def test_n50(self, seqtxt_output):
        output_statistics = seqtxt_output[1]
        n50_read_length = output_statistics.all_long_read_info.long_read_info.n50_read_length
        assert n50_read_length == 8733

    @pytest.mark.dependency(depends=["TestSeqTxt::test_success"])
    def test_passed_n50(self, seqtxt_output):
        output_statistics = seqtxt_output[1]
        passed_n50_read_length = output_statistics.passed_long_read_info.long_read_info.n50_read_length
        assert passed_n50_read_length == 7050


@pytest.fixture(scope='class')
def rnaseq_bam_output():
    """Run the BAM module on RNASeq inputs."""
    # Set parameters
    default_parameters = lrst.Input_Para()
    output_folder = os.path.abspath(str("output/"))
    default_parameters.output_folder = output_folder
    default_parameters.out_prefix = str("rnaseq_")
    default_parameters.tin_sample_size = 100
    default_parameters.tin_min_coverage = 2

    # Check if running remotely
    local_dir = os.path.expanduser('~/github/LongReadSum')
    if os.getcwd() == local_dir:
        input_file = os.path.join(local_dir, "SampleData/GTEX-RNASeq-subset.bam")  # Local path
        default_parameters.gene_bed = os.path.join(local_dir, "SampleData/gencode.v46.basic.subset.bed")
    else:
        input_file = os.path.abspath(str("SampleData/GTEX-RNASeq-subset.bam"))  # Remote path
        default_parameters.gene_bed = os.path.abspath(str("SampleData/gencode.v46.basic.subset.bed"))

    # Add input files
    default_parameters.add_input_file(input_file)

    # Run the BAM statistics module
    output = lrst.Output_BAM()
    exit_code = lrst.callBAMModule(default_parameters, output)

    yield [exit_code, output, input_file]


class TestRNASeqBAM:
    """Tests for RNASeq BAM inputs."""

    # Ensure the module ran successfully
    @pytest.mark.dependency()
    def test_success(self, rnaseq_bam_output):
        exit_code = rnaseq_bam_output[0]
        assert exit_code == 0

    # Tests
    @pytest.mark.dependency(depends=["TestRNASeqBAM::test_success"])
    def test_tin_count(self, rnaseq_bam_output):
        output_statistics = rnaseq_bam_output[1]
        input_file = rnaseq_bam_output[2]
        tin_count = output_statistics.getTINCount(input_file)
        assert tin_count == 9

    @pytest.mark.dependency(depends=["TestRNASeqBAM::test_success"])
    def test_tin_mean(self, rnaseq_bam_output):
        output_statistics = rnaseq_bam_output[1]
        input_file = rnaseq_bam_output[2]
        tin_mean = output_statistics.getTINMean(input_file)
        assert round(tin_mean, 1) == 60.2

    @pytest.mark.dependency(depends=["TestRNASeqBAM::test_success"])
    def test_tin_median(self, rnaseq_bam_output):
        output_statistics = rnaseq_bam_output[1]
        input_file = rnaseq_bam_output[2]
        tin_median = output_statistics.getTINMedian(input_file)
        assert round(tin_median, 1) == 70.0

    @pytest.mark.dependency(depends=["TestRNASeqBAM::test_success"])
    def test_tin_stddev(self, rnaseq_bam_output):
        output_statistics = rnaseq_bam_output[1]
        input_file = rnaseq_bam_output[2]
        tin_stddev = output_statistics.getTINStdDev(input_file)
        assert round(tin_stddev, 1) == 31.0
        