#!/bin/bash -l

# Script to set up Nextflow and a conda environment for running the WGS and MitoEdit Nextflow pipelines on
# Virtual Machine instances

# Download Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_25.7.0-2-Linux-x86_64.sh
# Install Miniconda with bash set to auto accept/yes
bash Miniconda3-py310_25.7.0-2-Linux-x86_64.sh -b -p $HOME/miniconda3

# Cleanup/remove installer
rm -rf Miniconda3-py310_25.7.0-2-Linux-x86_64.sh

# Add channels
conda config --add channels bioconda
conda config --add channels conda-forge

# Accept the Anaconda Terms of Service for channels
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# Install mamba and update conda environment
conda install mamba -n base -y
mamba update --all -n base -y

# Make 'bin' directory in home for storing some tools
mkdir -p $HOME/bin

# Download Java 21 and Nextflow for running Nextflow pipelines
wget https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.8%2B9/OpenJDK21U-jre_x64_linux_hotspot_21.0.8_9.tar.gz
tar -xvf OpenJDK21U-jre_x64_linux_hotspot_21.0.8_9.tar.gz
mv jdk-21.0.8+9-jre/ bin/
rm -rf OpenJDK21U-jre_x64_linux_hotspot_21.0.8_9.tar.gz

# Set JAVA_HOME and PATH for Java 21
export JAVA_HOME=${HOME}/bin/jdk-21.0.8+9-jre
export PATH="$JAVA_HOME/bin:$PATH"

# Now get Nextflow
wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v25.04.7/nextflow | bash
mv nextflow ${HOME}/bin/
chmod +x ${HOME}/bin/nextflow

# Now put bin on PATH
export PATH="${HOME}/bin:$PATH"

# Now install conda environments necessary for running Mitochondrial Variant Calling Nextflow pipeline
mamba env create -f nf_integrated/Mito_WGS_Nextflow/conda_environment/nf_mitvar.yml
mamba env create -f nf_integrated/Mito_WGS_Nextflow/conda_environment/annot_python2.yml

# Install mitoedit in its own conda environment
mamba env create -f nf_integrated/nf_mte/mte.yml

# Finally install streamlit and igv-reports packages for front-end UI and generating IGV report vis
mamba env create -f igv/ends.yml

# Create empty file for bystander input
touch NOFILE

# Uncompress WGS annovar data
gunzip -v nf_integrated/Mito_WGS_Nextflow/ANNOVAR/annovar_humandb/*.gz

# # Get mitochondrial reference
# mamba install -n base -y entrez-direct
# esearch -db nucleotide -query "NC_012920.1" | efetch -format fasta > NC_012920.1.fasta
