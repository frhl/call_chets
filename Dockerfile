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

# Create version file (will be embedded in binaries through makef# Build
RUN make
# No need to move binaries as make now outputs to bin/ and Dockerfile WORKDIR is /usr/src/app
# However, if we want them in path:
ENV PATH="/usr/src/app/bin:${PATH}"
RUN mv bin/transform /usr/local/bin/.
RUN mv bin/count_by_gene /usr/local/bin/.
RUN mv bin/recode_vcf /usr/local/bin/.


# Set default command to R when the container starts
#CMD ["bash"]

