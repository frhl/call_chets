# Version determination
VERSION_FILE := $(shell cat .version 2>/dev/null)

# Try to get git info regardless of .version file
GIT_COMMIT ?= $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GIT_DATE ?= $(shell git log -1 --format=%cd --date=short 2>/dev/null || date +%Y-%m-%d)

ifneq ($(VERSION_FILE),)
    GIT_VERSION := $(VERSION_FILE)
else
    GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags 2>/dev/null || date +%Y-%m-%d)
endif

# Compiler
CXX := g++

# Compiler flags
CXXFLAGS := -O2 -std=c++11 -Wall -Wextra -Wpedantic -Wshadow -DVERSION=\"$(GIT_VERSION)\" -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DGIT_DATE=\"$(GIT_DATE)\"

# Include directories
INCLUDES := -I./src

# Library paths and libraries
# Add Homebrew paths on macOS if they exist
ifeq ($(shell uname), Darwin)
    ifneq (,$(wildcard /opt/homebrew/include))
        INCLUDES += -I/opt/homebrew/include
        LIBS_PATH += -L/opt/homebrew/lib
    endif
    ifneq (,$(wildcard /usr/local/include))
        INCLUDES += -I/usr/local/include
        LIBS_PATH += -L/usr/local/lib
    endif
endif

LIBS := $(LIBS_PATH) -lz -lhts

# Output directory
BIN_DIR := bin

# Target executables
TARGETS := bin/interpret_phase bin/make_pseudo_vcf bin/recode bin/filter_pp bin/count_by_gene

# Legacy symlinks (for backward compatibility)
LEGACY_LINKS := bin/call_chets bin/encode_vcf bin/transform bin/filter_vcf_by_pp bin/orthogonalize

# Source files
SRC_DIR := src
SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(SRCS:.cpp=.o)

# Default target
all: $(TARGETS) legacy_links

# Link targets
bin/interpret_phase: src/interpret_phase.o src/ChetCaller.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/make_pseudo_vcf: src/make_pseudo_vcf.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/recode: src/recode.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/filter_pp: src/filter_pp.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/count_by_gene: src/count_by_gene.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Compile source files
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create legacy symlinks
legacy_links: $(TARGETS)
	@ln -sf interpret_phase bin/call_chets
	@ln -sf make_pseudo_vcf bin/encode_vcf
	@ln -sf recode bin/transform
	@ln -sf recode bin/orthogonalize
	@ln -sf filter_pp bin/filter_vcf_by_pp

# Clean
clean:
	rm -f $(TARGETS) $(LEGACY_LINKS) src/*.o

# Install
PREFIX ?= /usr/local
install: all
	@mkdir -p $(PREFIX)/bin
	@cp $(TARGETS) $(PREFIX)/bin/
	@ln -sf interpret_phase $(PREFIX)/bin/call_chets
	@ln -sf make_pseudo_vcf $(PREFIX)/bin/encode_vcf
	@ln -sf recode $(PREFIX)/bin/transform
	@ln -sf recode $(PREFIX)/bin/orthogonalize
	@ln -sf filter_pp $(PREFIX)/bin/filter_vcf_by_pp

.PHONY: all clean legacy_links install
