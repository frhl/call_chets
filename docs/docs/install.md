## Installation

### Quick start with Docker

The easiest way to run **arcade** is via Docker:

```bash
# Pull the image
docker pull fhlassen/arcade:latest

# Run a tool (mount current directory to /data)
docker run -v $PWD:/data fhlassen/arcade interpret_phase --help
```

All tools are available in the Docker image. To run commands on local files, mount
the relevant directory:

```bash
docker run -v $PWD:/data fhlassen/arcade \
    interpret_phase \
        --geno /data/genotypes.txt.gz \
        --gene-map /data/gene_map.txt \
        > results.txt
```

### Building from source

#### Pre-requisites

- **g++** with C++11 support
- **htslib** (>= 1.19 recommended)
- **zlib**
- **make**

On **macOS** (Homebrew):
```bash
brew install htslib
```

On **Ubuntu/Debian**:
```bash
sudo apt-get install libhts-dev zlib1g-dev
```

#### Compile and install

```bash
git clone https://github.com/frhl/call_chets.git
cd call_chets
make
sudo make install
```

This installs the following executables to `/usr/local/bin`:

| Executable | Description |
|------------|-------------|
| `interpret_phase` | Compound heterozygote caller |
| `make_pseudo_vcf` | Pseudo-variant VCF generator |
| `recode` | Genotype recoder/orthogonalizer |
| `filter_pp` | Posterior probability filter |
| `count_by_gene` | Per-gene genotype counter |

To install to a custom location:
```bash
make install PREFIX=/path/to/install
```

#### Verify installation

```bash
interpret_phase --help
recode --help
```

### External dependencies

**arcade** requires [bcftools](https://samtools.github.io/bcftools/) for the genotype extraction step of the gene-level pipeline. bcftools is not bundled with **arcade** but is available via most package managers:

```bash
# macOS
brew install bcftools

# Ubuntu/Debian
sudo apt-get install bcftools
```
