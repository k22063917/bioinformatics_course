#!/usr/bin/env bash

# Make a directory that I would conduct my assignment
mkdir K22063917_Adv_Bioinformatics_Assignment

# Getting and installing Conda
# Download Miniconda3 latest version on Open Stack
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

#Then run this bash script to install miniconda3
bash Miniconda3-latest-Linux-x86_64.sh

# To reload our current shell session
source ~/.bashrc

# Create a new conda environment
# to install new packages within the new-env
# environment location: /home/ubuntu/miniconda3/envs/new-env
conda create -n new-env

#to enter the new-env
conda activate new-env

#Then install the packages required
conda install trimmomatic
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install fastqc


#vcffilter already available within miniconda, no need to install

