CC = gcc
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lhts
TARGET = get_phased_sites
SOURCE = get_phased_sites.cpp

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET)


