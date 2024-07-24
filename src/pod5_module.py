import os
import numpy as np
import logging
import pod5 as p5

# Import the longreadsum module
if __name__ == 'src.pod5_module':
    from lib import lrst  # For debugging
else:
    import lrst

def generate_pod5_qc(input_data: dict) -> dict:
    """Generate QC for POD5"""
    logging.info("Generating QC for POD5")

    # Get the list of read IDs to process (if specified)
    read_id_list = []
    if (input_data['read_ids']):
        # Parse the comma-separated list of read IDs
        read_id_list = input_data['read_ids'].split(',')
        logging.info("Processing %d read IDs", len(read_id_list))

    # Loop through each input file and get the QC data across files
    # Dictionary to store read signal QC (key: read_id, [signal, mean_signal, median_signal, std_signal, skewness, kurtosis])
    read_signal_dict = {}
    for input_file in input_data['input_files']:
        logging.info("File name: %s", input_file)

        # Iterate through each read in the file
        with p5.Reader(input_file) as reader:
            for read in reader:
                # Check if the read is in the list of reads to process (if
                # provided)
                # Convert the UUID to a string
                read_id = str(read.read_id)
                if read_id_list and read_id not in read_id_list:
                    logging.info("Skipping read ID: %s", read_id)
                    continue
                # logging.info("Processing read ID: %s", read_id)

                # Get the basecall signals
                read_signal = read.signal

                read_signal_dict[read_id] = {}
                # logging.info("Adding read signals")
                read_signal_dict[read_id]['signal'] = read_signal

                # Calculate QC metrics using numpy
                # Mean signal
                mean_signal = np.mean(read_signal)
                read_signal_dict[read_id]['mean'] = mean_signal

                # Median signal
                median_signal = np.median(read_signal)
                read_signal_dict[read_id]['median'] = median_signal

                # Standard deviation
                std_signal = np.std(read_signal)
                read_signal_dict[read_id]['std'] = std_signal

                # Pearson skewness coefficient
                moment3 = np.mean((read_signal - mean_signal) ** 3)
                moment2 = np.mean((read_signal - mean_signal) ** 2)
                skewness = moment3 / (moment2 ** 1.5)
                read_signal_dict[read_id]['skewness'] = skewness

                # Kurtosis
                moment4 = np.mean((read_signal - mean_signal) ** 4)
                kurtosis = moment4 / (moment2 ** 2)
                read_signal_dict[read_id]['kurtosis'] = kurtosis

                
    return read_signal_dict
