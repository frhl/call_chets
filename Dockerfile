# Use an official base image, e.g., Ubuntu
FROM ubuntu:20.04

# Set maintainer label
LABEL maintainer="flassen@well.ox.ac.uk"

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
WORKDIR htslib
RUN make
RUN make install
WORKDIR ..

# copy scripts
WORKDIR app
COPY makefile makefile
COPY get_non_ref_sites.cpp get_non_ref_sites.cpp
COPY call_chets.cpp call_chets.cpp
COPY encode_vcf.cpp encode_vcf.cpp
COPY .version .version
RUN make

# move to folder in PATH
RUN mv call_chets /usr/local/bin/.
RUN mv get_non_ref_sites /usr/local/bin/.
RUN mv encode_vcf /usr/local/bin/.

# Set default command to R when the container starts
#CMD ["bash"]


