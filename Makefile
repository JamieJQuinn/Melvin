CC=g++
CFLAGS=-c -std=c++11
LDFLAGS=
LIN_SOURCES=boxLinear.cpp 
LIN_OBJECTS=$(LIN_SOURCES:.cpp=.o)
LIN_EXECUTABLE=boxLinear

all: $(LIN_SOURCES) $(LIN_EXECUTABLE)
    
linear: $(LIN_EXECUTABLE)

$(LIN_EXECUTABLE): $(LIN_OBJECTS) 
	$(CC) $(LDFLAGS) $(LIN_OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(LIN_EXECUTABLE)

