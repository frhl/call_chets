# Version determination: Try .version file first (for Docker), then git, then date
VERSION_FILE := $(shell cat .version 2>/dev/null)
ifneq ($(VERSION_FILE),)
    GIT_VERSION := $(VERSION_FILE)
else
    GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags 2>/dev/null || date +%Y-%m-%d)
endif

# Compiler
CXX := g++

# Compiler flags
CXXFLAGS := -O2 -std=c++11 -Wall -DVERSION=\"$(GIT_VERSION)\"

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

# Output binaries
TARGET_CALL_CHETS := $(BIN_DIR)/call_chets
TARGET_ENCODE := $(BIN_DIR)/encode_vcf
TARGET_RECODE := $(BIN_DIR)/recode_vcf
TARGET_TRANSFORM := $(BIN_DIR)/transform
# Target executables
TARGETS := bin/interpret_phase bin/make_pseudo_vcf bin/orthogonalize bin/filter_pp bin/count_by_gene bin/encode_vcf_by_group

# Legacy symlinks (for backward compatibility)
LEGACY_LINKS := bin/call_chets bin/encode_vcf bin/transform bin/filter_vcf_by_pp

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

bin/orthogonalize: src/orthogonalize.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/filter_pp: src/filter_pp.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/count_by_gene: src/count_by_gene.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

bin/encode_vcf_by_group: src/encode_vcf_by_group.o
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LIBS)

# Compile source files
src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Create legacy symlinks
legacy_links: $(TARGETS)
	@ln -sf interpret_phase bin/call_chets
	@ln -sf make_pseudo_vcf bin/encode_vcf
	@ln -sf orthogonalize bin/transform
	@ln -sf filter_pp bin/filter_vcf_by_pp

# Clean
clean:
	rm -f $(TARGETS) $(LEGACY_LINKS) src/*.o

# Install
PREFIX ?= /usr/local
install: all
	@mkdir -p $(PREFIX)/bin
	@cp $(TARGETS) $(PREFIX)/bin/
	@# Copy symlinks if desired, or let user rely on new names

.PHONY: all clean legacy_links install
