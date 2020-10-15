FROM ubuntu:20.04
FROM python:latest

LABEL author="Greydon Gilmore <greydon.gilmore@gmail.com>"

ENV DEBIAN_FRONTEND noninteractive

############# SYSTEM LEVEL INSTALLS #############
# install basic ubuntu utilities
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    apt-utils build-essential wget \
    python dialog git g++ gcc unzip \
    curl libjpeg62 libopenjp2-7 libopenjp2-7-dev \
    pkg-config cmake && \
	apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y

# Python dependencies
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install snakemake mne-bids numpy scipy mne \
	dicom2nifti pydicom dcmstack pandas heudiconv


################### dcm2niix ###################
# Install Dependencies
RUN apt-get update && apt-get upgrade -y && \
	apt-get install -y pigz && \
	apt-get clean -y && apt-get autoclean -y && apt-get autoremove -y

# Get dcm2niix from github and compile
RUN cd /tmp && \
	curl -fLO https://github.com/rordenlab/dcm2niix/releases/latest/download/dcm2niix_lnx.zip && \
	unzip dcm2niix_lnx.zip && \
	mv ./dcm2niix /usr/local/bin


############### NVM and Node ###############
ENV NODE_VERSION=12.8.0
RUN curl -o- https://raw.githubusercontent.com/creationix/nvm/v0.34.0/install.sh | bash
ENV NVM_DIR=/root/.nvm
RUN . "$NVM_DIR/nvm.sh" && nvm install ${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm use v${NODE_VERSION}
RUN . "$NVM_DIR/nvm.sh" && nvm alias default v${NODE_VERSION}
ENV PATH="/root/.nvm/versions/node/v${NODE_VERSION}/bin/:${PATH}"
RUN node --version
RUN npm --version


############### bids-validator ###############
# bids-validator
RUN npm i -g bids-validator


############# KEEP BELOW SYSTEM LEVEL INSTALLS #############
# setup working directories
#WORKDIR /d2b-clinical

# copy over data files
#COPY ./data /data
#COPY ./workflow /d2b-clinical/workflow
#COPY ./config/config_docker.yaml /d2b-clinical/config/config.yaml
