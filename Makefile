CC=g++
CFLAGS=-c -std=c++0x
LDFLAGS=
SRC_DIR=src
BUILD_DIR=build
INCLUDE_DIR=include

SOURCES=$(wildcard $(SRC_DIR)/*.cpp)
OBJECTS=$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
EXECUTABLE=exe

all: $(BUILD_DIR) build

$(BUILD_DIR)/$(EXECUTABLE): $(OBJECTS)
			$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
			$(CC) $(CFLAGS) $< -o $@

$(BUILD_DIR):
			mkdir -p $@

.PHONY: clean
clean:
			rm -rf $(BUILD_DIR)

build: CFLAGS += -DNDEBUG -O2 -fopenmp
build: LDFLAGS += -fopenmp
build: $(BUILD_DIR)/$(EXECUTABLE)

debug: CFLAGS += -DDEBUG -g -pg -Wall
debug: LDFLAGS += -pg
debug: $(BUILD_DIR)/$(EXECUTABLE)

ddcLinear: CFLAGS += -DLINEAR -DDDC
ddcLinear: build

linear: CFLAGS += -DLINEAR
linear: build

linearDebug: CFLAGS += -DLINEAR
linearDebug: debug

nonlinear: CFLAGS += -DNONLINEAR
nonlinear: build

nonlinearDebug: CFLAGS += -DNONLINEAR
nonlinearDebug: debug
