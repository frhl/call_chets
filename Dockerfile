# Use an official base image, e.g., Ubuntu
FROM ubuntu:20.04

# Set maintainer label
LABEL maintainer="flassen@well.ox.ac.uk"

# Update and install system dependencies
RUN apt-get update && apt-get install -y \
    g++ \
    zlib1g-dev \
    git

# Clone htslib and install
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git
WORKDIR htslib

# copy scripts 
WORKDIR ../app
COPY makefile makefile
COPY get_non_ref_sites.cpp get_non_ref_sites.cpp
COPY call_chets.cpp call_chets.cpp
COPY encode_vcf.cpp enxode_vcf.cpp
RUN make


# Set default command to R when the container starts
CMD ["bash"]

