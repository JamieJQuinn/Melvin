CC=g++
CFLAGS=-c -std=c++0x -fopenmp
DEBUG_CFLAGS= -DDEBUG -g -pg -Wall
BUILD_CFLAGS= -DNDEBUG -O2
LDFLAGS=-fopenmp
DEBUG_LDFLAGS=-pg
SOURCES=box.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=box

all: $(SOURCES) debug

test: CFLAGS += -DNDEBUG -pg
test: LDFLAGS += -pg 
test: $(EXECUTABLE)

build: CFLAGS += $(BUILD_CFLAGS)
build: $(EXECUTABLE)

debug: CFLAGS += $(DEBUG_CFLAGS) 
debug: LDFLAGS += $(DEBUG_LDFLAGS) 
debug: $(EXECUTABLE)

linear: CFLAGS += -DLINEAR
linear: build

linearDebug: CFLAGS += -DLINEAR
linearDebug: debug

nonlinear: CFLAGS += -DNONLINEAR
nonlinear: build

nonlinearDebug: CFLAGS += -DNONLINEAR
nonlinearDebug: debug

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)

