# LongReadSum: A fast and flexible QC tool for long read sequencing data

![image](https://github.com/user-attachments/assets/ae3d2b6a-46da-476a-8fe5-7e0a34751585)

[![build tests](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml/badge.svg)](https://github.com/WGLab/LongReadSum/actions/workflows/build-test.yml)

LongReadSum supports FASTA, FASTQ, BAM, FAST5, and sequencing_summary.txt file formats for quick generation of QC data in HTML and text format.

# README Contents
- [Installation using Anaconda (recommended)](#installation-using-anaconda)
- [Installation using Docker](#installation-using-anaconda)
- [Building from source](#building-from-source)
- General usage for common filetypes:
  - [Common parameters](#common-parameters)
  - [WGS BAM](#wgs-bam)
  - [BAM with base modifications/methylation](#bam-with-base-modifications)
  - [RRMS BAM](#rrms-bam)
  - [PacBio unaligned BAM](#pacbio-unaligned-bam)
  - [RNA-Seq BAM (TIN values)](#rna-seq-bam)
  - [ONT POD5](#ont-pod5)
  - [ONT FAST5](#ont-fast5)
    - [Signal QC](#signal-qc)
    - [Sequence QC](#sequence-qc)
  - [Basecall summary (ONT sequencing_summary.txt)](#basecall-summary)
  - [FASTQ](#fastq)
  - [FASTA](#fasta)
- [Revision history](#revision-history)
- [Getting help](#getting-help)
- [Citing LongReadSum](#citing-longreadsum)

# Installation using Anaconda
First, install [Anaconda](https://www.anaconda.com/).

Next, create a new environment. This installation has been tested with Python 3.10, Linux 64-bit.

```
conda create -n longreadsum python=3.9
conda activate longreadsum
```

LongReadSum and its dependencies can then be installed using the following command:

```
conda install -c wglab -c conda-forge -c jannessp -c bioconda longreadsum=1.4.0
```

# Installation using Docker
First, install [Docker](https://docs.docker.com/engine/install/).
Pull the latest image from Docker hub, which contains the latest longreadsum release and its dependencies.

```
docker pull genomicslab/longreadsum
```

## Running

On Unix/Linux:
```
docker run -v C:/Users/.../DataDirectory:/mnt/ -it genomicslab/longreadsum bam -i /mnt/input.bam -o /mnt/output
```
Note that the `-v` command is required for Docker to find the input file. Use a directory under `C:/Users/` to ensure volume files are mounted correctly. In the above example, the local directory `C:/Users/.../DataDirectory` containing the input file `input.bam` is mapped to a directory `/mnt/` in the Docker container. Thus, the input file and output directory arguments are relative to the `/mnt/` directory, but the output files will also be saved locally in `C:/Users/.../DataDirectory` under the specified subdirectory `output`.


# Building from source
To get the latest updates in longreadsum, you can build from source.
First install [Anaconda](https://www.anaconda.com/). Then follow the instructions below to install LongReadSum and its dependencies:

```
# Pull the latest updates
git clone https://github.com/WGLab/LongReadSum
cd LongReadSum

# Create the longreadsum environment, install dependencies, and activate
conda env create -f environment.yml
conda activate longreadsum

# Build the program
make
```


## Running
Activate the conda environment and then run with arguments:
```
conda activate longreadsum
longreadsum <FILETYPE> [arguments]
```

# General Usage

Specify the filetype followed by parameters:
```
longreadsum <FILETYPE> -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# Common parameters

To see all parameters for a filetype, run:

```longreadsum <FILETYPE> --help```

This section describes parameters common to all filetypes:

| Parameter	| Description | Default |
| --- | --- | --- |
| -i, --input | A single input filepath
| -I, --inputs | Multiple comma-separated input filepaths
| -P, --pattern | Use pattern matching (*) to specify multiple input files. Enclose the pattern in double quotes.
| -g, --log | Log file path | log_output.log
| -G, --log-level |Logging level (1: DEBUG, 2: INFO, 3: WARNING, 4: ERROR, 5: CRITICAL) | 2
| -o, --outputfolder | Output directory | output_longreadsum
| -t, --threads | The number of threads used | 1
| -Q, --outprefix | Output file prefix |  QC_

# WGS BAM

This section describes how to generate QC reports for BAM files from whole-genome sequencing
(WGS) with alignments to a linear reference genome such as GRCh38 (data shown is HG002 sequenced with ONT Kit V14
Promethion R10.4.1 from https://labs.epi2me.io/askenazi-kit14-2022-12/)

![image](https://github.com/user-attachments/assets/166f6d04-26ca-4469-be2c-ce466597a68a)

![image](https://github.com/user-attachments/assets/7d83e55c-85a2-48a8-b9a7-92b671de758f)

![image](https://github.com/user-attachments/assets/d303f5a9-8e1b-425e-b0b0-f46f263a3f9f)

![image](https://github.com/user-attachments/assets/f74f985a-3c3d-4b00-bf98-d59b128d8722)

## General usage
```
longreadsum bam -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# BAM with base modifications

This section describes how to generate QC reports for BAM files with MM, ML base modification tags (data shown is HG002 sequenced with ONT
MinION R9.4.1 from https://labs.epi2me.io/gm24385-5mc/)

![image](https://github.com/user-attachments/assets/5d97a949-842a-4f41-bfc5-81e9f30c57bc)


## Parameters
| Parameter	| Description | Default |
| --- | --- | --- |
| --mod | Run base modification analysis on the BAM file | False
| --modprob | Base modification filtering threshold. Above/below this value, the base is considered modified/unmodified. | 0.8
| --ref | The reference genome FASTA file to use for identifying CpG sites (optional)

## General usage
```
longreadsum bam -i $INPUT_FILE -o $OUTPUT_DIRECTORY --mod --modprob 0.8 --ref $REF_GENOME
```

# RRMS BAM

This section describes describes how to generate QC reports for ONT RRMS BAM files and associated CSVs (data shown is HG002 RRMS using ONT
R9.4.1).

### Accepted reads:
![image](https://github.com/user-attachments/assets/c0e69e53-0a1e-432d-ad4c-9edfac764514)

![image](https://github.com/user-attachments/assets/105a47ff-7bd8-436e-9d3d-1b112b94fb5e)


### Rejected reads:
![image](https://github.com/user-attachments/assets/7c213975-ec6b-4476-81c9-8853c664b653)

![image](https://github.com/user-attachments/assets/604ca74a-516b-48a7-8b02-931e27255bd8)


## Parameters
| Parameter	| Description | Default |
| --- | --- | --- |
| -c, --csv | CSV file containing read IDs to extract from the BAM file*

The CSV file should contain a `read_id` column with the read IDs in the BAM
file, and a `decision` column with the accepted/rejected status of the read.
Accepted reads will have `stop_receiving` in the `decision` column, while rejected
reads will have `unblock`:

```
batch_time,read_number,channel,num_samples,read_id,sequence_length,decision
1675186897.6034577,93,4,4011,f943c811-3f97-4971-8aed-bb9f36ffb8d1,361,unblock
1675186897.7544408,80,68,4025,fab0c19d-8085-454c-bfb7-c375bbe237a1,462,unblock
1675186897.7544408,93,127,4028,5285e0ba-86c0-4b5d-ba27-5783acad6105,438,unblock
1675186897.7544408,103,156,4023,65d8befa-eec0-4496-bf2b-aa1a84e6dc5e,362,stop_receiving
...
```

## General usage
```
longreadsum rrms -i $INPUT_FILE -o $OUTPUT_DIRECTORY -c $RRMS_CSV
```

# RNA-Seq BAM

This section describes how to generate QC reports for TIN (transcript integrity
number) scores from RNA-Seq BAM files (data shown is Adult GTEx v9 long-read RNA-seq data sequenced with ONT
cDNA-PCR protocol from https://www.gtexportal.org/home/downloads/adult-gtex/long_read_data).

## Outputs
A TSV file with scores for each transcript:

```
geneID	chrom	tx_start	tx_end	TIN
ENST00000456328.2	chr1	11868	14409	2.69449577083296
ENST00000450305.2	chr1	12009	13670	0.00000000000000
ENST00000488147.2	chr1	14695	24886	94.06518975035769
ENST00000619216.1	chr1	17368	17436	0.00000000000000
ENST00000473358.1	chr1	29553	31097	0.00000000000000
...
```

An TSV file with TIN score summary statistics:

```
Bam_file	TIN(mean)	TIN(median)	TIN(stddev)
/mnt/isilon/wang_lab/perdomoj/data/GTEX/GTEX-14BMU-0526-SM-5CA2F_rep.FAK93376.bam	67.06832655372376	74.24996965188242	26.03788585287367
```

A summary table in the HTML report:

![image](https://github.com/user-attachments/assets/400bcd68-05fc-4f08-8b70-b981cd9dc994)

## Parameters
| Parameter	| Description | Default |
| --- | --- | --- |
| --genebed | Gene BED12 file required for calculating TIN scores
| --sample-size | Sample size for TIN calculation | 100
| --min-coverage | Minimum coverage for TIN calculation | 10

## General usage
```
longreadsum bam -i $INPUT_FILE -o $OUTPUT_DIRECTORY --genebed $BED_FILE --min-coverage <COVERAGE> --sample-size <SIZE>
```

Download an example HTML report [here]() (data is Adult GTEx v9 long-read RNA-seq data sequenced with ONT
cDNA-PCR protocol from https://www.gtexportal.org/home/downloads/adult-gtex/long_read_data)

# PacBio unaligned BAM

This section describes how to generate QC reports for PacBio BAM files without alignments (data shown is HG002 sequenced with PacBio
Revio HiFi long reads obtained from https://www.pacb.com/connect/datasets/#WGS-datasets).

![image](https://github.com/user-attachments/assets/76374274-3671-49d2-984f-0208e9d8e3e7)

![image](https://github.com/user-attachments/assets/15112738-b6cd-4d1d-b0c4-e0bb31464374)

![image](https://github.com/user-attachments/assets/e3935f58-eb53-4f9d-b4b5-7287fcdc3252)

![image](https://github.com/user-attachments/assets/8b17c9e2-8932-45b3-a673-b5b35ae994e6)

## General usage
```
longreadsum bam -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# ONT POD5

This section describes how to generate QC reports for ONT POD5 (signal) files and their corresponding basecalled BAM files (data shown is HG002 using ONT
R10.4.1 and LSK114 downloaded from the tutorial https://github.com/epi2me-labs/wf-basecalling).

![image](https://github.com/user-attachments/assets/62c3c810-5c1a-4124-816b-74245af8b57c)


## Parameters
> [!NOTE]
> The interactive signal-base correspondence plots in the HTML report use a
lot of memory (RAM) which can make your web browser slow. Thus by default, we
randomly sample only a few reads, and the user can specify a list of read IDs as
well (e.g. from a specific region of interest).

| Parameter	| Description | Default |
| --- | --- | --- |
| -b, --basecalls | The basecalled BAM file to use for signal extraction
| -r, --read_ids | A comma-separated list of read IDs to extract from the file
| -R, --read-count | Set the number of reads to randomly sample from the file | 3

## General usage
```
# Individual file:
longreadsum pod5 -i $INPUT_FILE -o $OUTPUT_DIRECTORY --basecalls $INPUT_BAM [--read-count <COUNT> | --read-ids <IDS>]

# Directory:
longreadsum pod5 -P "$INPUT_DIRECTORY/*.fast5" -o $OUTPUT_DIRECTORY --basecalls $INPUT_BAM [--read-count <COUNT> | --read-ids <IDS>]
```

# ONT FAST5

## Signal QC

This section describes how to generate QC reports for generating a signal and basecalling QC
report from ONT FAST5 files with signal and basecall information (data shown is HG002 sequenced with ONT MinION R9.4.1 from https://labs.epi2me.io/gm24385-5mc/)

![image](https://github.com/user-attachments/assets/b2a56529-55e3-4678-918a-e5787881643e)


## Parameters
> [!NOTE]
> The interactive signal-base correspondence plots in the HTML report use a
lot of memory (RAM) which can make your web browser slow. Thus by default, we
randomly sample only a few reads, and the user can specify a list of read IDs as
well (e.g. from a specific region of interest).

| Parameter	| Description | Default |
| --- | --- | --- |
| -r, --read_ids | A comma-separated list of read IDs to extract from the file
| -R, --read-count | Set the number of reads to randomly sample from the file | 3

## General usage
```
# Individual file:
longreadsum f5s -i $INPUT_FILE -o $OUTPUT_DIRECTORY [--read-count <COUNT> | --read-ids <IDS>]

# Directory:
longreadsum f5s -P "$INPUT_DIRECTORY/*.fast5" -o $OUTPUT_DIRECTORY [--read-count <COUNT> | --read-ids <IDS>]
```

## Sequence QC

This section describes how to generate QC reports for sequence data from ONT FAST5 files (data shown is HG002 sequenced with ONT MinION R9.4.1 from https://labs.epi2me.io/gm24385-5mc/)

![image](https://github.com/user-attachments/assets/97876343-cd34-4bfe-9612-7f6b14a2be0d)

![image](https://github.com/user-attachments/assets/be8415a5-63ee-403b-931a-66d79a3b28a5)

![image](https://github.com/user-attachments/assets/2d822263-d9a4-470f-aa6d-5bccf570edac)

![image](https://github.com/user-attachments/assets/0fd1ee15-1e2f-492a-b072-bfaf768a448a)

![image](https://github.com/user-attachments/assets/aa5b37c8-7c83-418b-b5ea-6e53d557fb93)


## General usage
```
longreadsum f5 -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# Basecall summary

This section describes how to generate QC reports for ONT basecall summary (sequencing_summary.txt) files (data shown is HG002 sequenced with ONT
PromethION R10.4 from https://labs.epi2me.io/gm24385_q20_2021.10/, filename `gm24385_q20_2021.10/analysis/20210805_1713_5C_PAH79257_0e41e938/guppy_5.0.15_sup/sequencing_summary.txt`)

![image](https://github.com/user-attachments/assets/ad094b0a-7878-4937-840c-ad0d7c09335b)

![image](https://github.com/user-attachments/assets/5e2417e8-74b8-4f39-8c3d-d6481749711d)

![image](https://github.com/user-attachments/assets/f25841bf-8129-41bc-a90f-0196ca14159f)


## General usage
```
longreadsum seqtxt -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# FASTQ

This section describes how to generate QC reports for FASTQ files (data shown is HG002 ONT 2D from GIAB
 [FTP index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_HG002_Cornell_Oxford_Nanopore_fasta_fastq_10132015.HG002))

![image](https://github.com/user-attachments/assets/09d80928-bae7-4c64-a21c-8ef21fe9ab60)

![image](https://github.com/user-attachments/assets/d5ee0aa5-9127-447f-b96b-26f3fad7a963)

![image](https://github.com/user-attachments/assets/cea3c23c-44a3-4313-9d31-4c8559073b22)

![image](https://github.com/user-attachments/assets/acb199f6-4529-43ce-9212-f938128b0706)

![image](https://github.com/user-attachments/assets/47e395fc-b33d-45d0-b3ac-c658a84f62cb)


## General usage
```
longreadsum fq -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# FASTA

This section describes how to generate QC reports for FASTA files (data shown is HG002 ONT 2D from GIAB
 [FTP index](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_indexes/AshkenazimTrio/sequence.index.AJtrio_HG002_Cornell_Oxford_Nanopore_fasta_fastq_10132015.HG002)).

![image](https://github.com/user-attachments/assets/d4862e6d-435e-4317-b331-4af0428a6419)

![image](https://github.com/user-attachments/assets/af3b736b-beb6-44e4-a0d8-df736c288389)

![image](https://github.com/user-attachments/assets/cd4e1f59-0c34-41a6-91ea-08381bdc906a)

## General usage
```
longreadsum fa -i $INPUT_FILE -o $OUTPUT_DIRECTORY
```

# Revision history
For release history, please visit [here](https://github.com/WGLab/LongReadSum/releases). 

# Getting help
Please refer to the [LongReadSum issue pages](https://github.com/WGLab/LongReadSum/issues) for posting your issues. We will also respond your questions quickly. Your comments are criticl to improve our tool and will benefit other users.

# Citing LongReadSum
### Please cite the preprint below if you use our tool:

Perdomo, J. E., Ahsan, M. U., Liu, Q., Fang, L. & Wang, K. LongReadSum: A fast and flexible quality control and signal summarization tool for long-read sequencing data. bioRxiv, 2024.2008.2005.606643, [doi:10.1101/2024.08.05.606643](https://doi.org/doi:10.1101/2024.08.05.606643) (2024).
