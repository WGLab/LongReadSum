# Use the miniconda container
FROM continuumio/miniconda3:main

WORKDIR /app

RUN apt-get update
RUN conda update conda

# Install LongReadSum
RUN conda config --add channels wglab
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda config --add channels jannessp
RUN conda create -n longreadsum python=3.9
RUN echo "conda activate longreadsum" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
RUN conda install -n longreadsum -c wglab -c conda-forge -c jannessp -c bioconda longreadsum=1.5.0 && conda clean -afy

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "longreadsum", "longreadsum"]
