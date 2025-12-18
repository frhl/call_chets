# Use an official base image, e.g., Ubuntu
FROM ubuntu:20.04

# Set maintainer label
LABEL maintainer="flassen@well.ox.ac.uk"

# Define build arguments for git information
ARG GIT_COMMIT=unknown
ARG GIT_DATE=unknown
# For DockerHub automated builds, these can be populated from build hooks

# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    make \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libhts-dev \
    git

# Clone htslib
WORKDIR /usr/src
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git

# Install htslib

# Install htslib
WORKDIR /usr/src/htslib
RUN make
RUN make install
WORKDIR /usr/src

# Fix lhts paths
RUN ldconfig

# copy scripts
WORKDIR /usr/src/app
COPY makefile makefile
COPY src src
COPY scripts scripts
RUN mkdir bin

# Set execute permission for make_version.sh
RUN chmod +x scripts/make_version.sh

# Pass git information to the version script
ENV GIT_COMMIT=$GIT_COMMIT
ENV GIT_DATE=$GIT_DATE

# Build and Install
# Pass GIT_COMMIT and GIT_DATE variables to make
RUN make GIT_COMMIT=$GIT_COMMIT GIT_DATE=$GIT_DATE
RUN make install

# No need to manually move binaries or set PATH if installing to /usr/local/bin (default)
# But strictly, /usr/local/bin is usually in PATH.

# Set default command
CMD ["interpret_phase", "--help"]
