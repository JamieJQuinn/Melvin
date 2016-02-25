CC=g++
CFLAGS=-c -std=c++11 -g
LDFLAGS=
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

