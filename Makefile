CC=g++
CFLAGS=-c -Wall -O3 
LFLAGS = -O3 
all: peak_tool.o
	$(CC) $(LFLAGS) peak_tool.o -o peak_tool

peak_tool.o: peak_tool.cpp
	$(CC) $(CFLAGS)  peak_tool.cpp
clean:
	rm -rf *o peak_tool