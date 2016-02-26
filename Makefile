CC=g++
CFLAGS=-c -Wall -std=c++0x -g -fopenmp
LDFLAGS=-fopenmp
SOURCES=box.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=box

all: $(SOURCES) $(EXECUTABLE)
    
linear: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)

