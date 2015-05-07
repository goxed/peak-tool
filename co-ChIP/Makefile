CC=g++
CFLAGS=-c -Wall -O3 -march=native -mtune=native -flto
LFLAGS = -O3 -march=native -mtune=native -Wall -flto
all: peak_tool.o
	$(CC) $(LFLAGS) peak_tool.o -o peak_tool_cochip

peak_tool.o: peak_tool.cpp
	$(CC) $(CFLAGS)  peak_tool.cpp
clean:
	rm -rf *o peak_tool_cochip
