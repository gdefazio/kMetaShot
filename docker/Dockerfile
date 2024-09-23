FROM continuumio/miniconda3

# Set working directory 
WORKDIR /kmetashot

# Create Conda environment from the YAML file
COPY environment.yml .
RUN conda env create -f environment.yml

# Override default shell and use bash
SHELL ["conda", "run", "-n", "env", "/bin/bash", "-c"]

# Activate Conda environment and check if it is working properly
RUN echo "Making sure kmetashot is installed correctly..."
RUN kMetaShot_classifier_NV.py --help

ENTRYPOINT ["conda", "run", "--name", "env"]
