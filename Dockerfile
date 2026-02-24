# Use the official Miniconda3 lightweight image
FROM continuumio/miniconda3:latest

# Provide metadata about the image
LABEL maintainer="Pablo Vargas"
LABEL description="Docker image for SIREN: Suite for Intelligent RNAi Design"

# Install RNAhybrid from the Bioconda channel
RUN conda install -c conda-forge -c bioconda rnahybrid -y && \
    conda clean -a -y

# Install SIREN directly from PyPI (version 0.1.9)
RUN pip install --no-cache-dir siren-rnai==0.1.9

# Set the working directory inside the container
WORKDIR /data

# Make the container run like a native executable
ENTRYPOINT ["siren-rnai"]

# If no arguments are passed, display the help menu
CMD ["--help"]
