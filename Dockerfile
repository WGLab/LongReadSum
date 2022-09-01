# Use the miniconda container
FROM continuumio/miniconda3

# Copy the project directory
COPY . /longreadsum
WORKDIR /longreadsum
RUN ls

# Install build tools
RUN apt-get update && apt-get install build-essential -y

# Create the environment
RUN conda env create -f environment.yml

# Activate the environment
RUN echo "conda activate lrst_py39" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Ensure the correct environment is being used
RUN export PATH="/opt/conda/envs/lrst_py39/bin/python"
RUN which python

# Build LongReadSum
RUN make

# The code to run when container is started:
#NTRYPOINT ["python", "longreadsum"]
WORKDIR /
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "lrst_py39", "python", "longreadsum"]

