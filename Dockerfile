# Use the miniconda container
FROM continuumio/miniconda3

# Copy the project directory
COPY . /app/longreadsum
WORKDIR /app/longreadsum

# Install build tools
RUN apt-get update && apt-get install build-essential -y

# Install VBZ compression
RUN wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz
RUN tar -xf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz

# Create the environment
RUN conda env create -f environment.yml

# Activate the environment
RUN echo "conda activate longreadsum" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Ensure the correct environment is being used
RUN export PATH="/opt/conda/envs/longreadsum/bin/python"
RUN which python

# Build LongReadSum
RUN make

# Set up the HDF5 plugin path
ENV HDF5_PLUGIN_PATH="/longreadsum/lib/"

# The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "longreadsum", "python", "/app/longreadsum"]
