CC = g++
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lhts -lz

# Get version info from environment or git
GIT_COMMIT := $(shell if [ -n "$$GIT_COMMIT" ]; then echo "$$GIT_COMMIT"; else git rev-parse --short HEAD 2>/dev/null || echo "unknown"; fi)
GIT_DATE := $(shell if [ -n "$$GIT_DATE" ]; then echo "$$GIT_DATE"; else git log -1 --format=%cd --date=short 2>/dev/null || date +"%Y-%m-%d"; fi)
VERSION := $(shell cat .version 2>/dev/null || echo "0.3.0")

# Add git info to compilation flags
CFLAGS += -DGIT_COMMIT=\"$(GIT_COMMIT)\" -DGIT_DATE=\"$(GIT_DATE)\" -DVERSION=\"$(VERSION)\"

TARGET1 = filter_vcf_by_pp
TARGET2 = call_chets
TARGET3 = encode_vcf
TARGET4 = transform
TARGET5 = count_by_gene
TARGET6 = encode_vcf_by_group
TARGET7 = recode_vcf

SOURCE1 = filter_vcf_by_pp.cpp
SOURCE2 = call_chets.cpp
SOURCE3 = encode_vcf.cpp
SOURCE4 = transform.cpp
SOURCE5 = count_by_gene.cpp
SOURCE6 = encode_vcf_by_group.cpp
SOURCE7 = recode_vcf.cpp

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)

$(TARGET1): $(SOURCE1)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE1) -o $(TARGET1) $(LIBS)

$(TARGET2): $(SOURCE2)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE2) -o $(TARGET2) $(LIBS)

$(TARGET3): $(SOURCE3)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE3) -o $(TARGET3) $(LIBS)

$(TARGET4): $(SOURCE4)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE4) -o $(TARGET4) $(LIBS)

$(TARGET5): $(SOURCE5)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE5) -o $(TARGET5) $(LIBS)

$(TARGET6): $(SOURCE6)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE6) -o $(TARGET6) $(LIBS)

$(TARGET7): $(SOURCE7)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE7) -o $(TARGET7) $(LIBS)





clean:
		rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)

