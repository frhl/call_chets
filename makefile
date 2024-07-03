CC = g++
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lhts -lz
TARGET1 = filter_vcf_by_pp
TARGET2 = call_chets
TARGET3 = encode_vcf
TARGET4 = transform
TARGET5 = count_by_gene
TARGET6 = encode_vcf_by_group
SOURCE1 = filter_vcf_by_pp.cpp
SOURCE2 = call_chets.cpp
SOURCE3 = encode_vcf.cpp
SOURCE4 = transform.cpp
SOURCE5 = count_by_gene.cpp
SOURCE6 = encode_vcf_by_group.cpp

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6)

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




clean:
		rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6)

