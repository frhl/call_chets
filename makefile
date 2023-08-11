CC = g++
CFLAGS = -Wall -O2
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lhts
TARGET1 = get_phased_sites
TARGET2 = call_chets
TARGET3 = long_to_wide
SOURCE1 = get_phased_sites.cpp
SOURCE2 = call_chets.cpp
SOURCE3 = long_to_wide.cpp

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(SOURCE1)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE1) -o $(TARGET1) $(LIBS)

$(TARGET2): $(SOURCE2)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE2) -o $(TARGET2) $(LIBS)

$(TARGET3): $(SOURCE3)
		$(CC) $(CFLAGS) $(INCLUDES) $(SOURCE3) -o $(TARGET3) $(LIBS)

clean:
		rm -f $(TARGET1) $(TARGET2) $(TARGET3)

