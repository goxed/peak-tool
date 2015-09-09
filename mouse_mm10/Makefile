CC=g++
CFLAGS=-c -Wall -O3  -std=c++11
LFLAGS = -O3 
all: peak_tool.o
	$(CC) $(LFLAGS) peak_tool_multi.o -o peak_tool_multi

peak_tool.o: peak_tool_multi.cpp
	$(CC) $(CFLAGS)  peak_tool_multi.cpp
clean:
	rm -rf *o peak_tool_multi
