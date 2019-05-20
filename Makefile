CC=g++
CFLAGS=-c -std=c++0x -I$(INCLUDE_DIR)
LDFLAGS=
SRC_DIR=src
BUILD_DIR=build
INCLUDE_DIR=include

#SOURCES=src/main.cpp src/sim.cpp src/numerical_methods.cpp src/thomas_algorithm.cpp src/constants.cpp src/variable.cpp src/double_diffusive_sim.cpp
SOURCES=$(wildcard $(SRC_DIR)/*.cpp)
OBJECTS=$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
EXECUTABLE=exe

.PHONY: all
all: release

$(BUILD_DIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

$(BUILD_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

.PHONY: test
test:
	test/test_all_linear.sh
	test/test_all_nonlinear.sh

profile: CFLAGS += -pg
profile: LDFLAGS += -pg
profile: release

release: CFLAGS += -DNDEBUG -O2 -fopenmp
release: LDFLAGS += -fopenmp
release: $(BUILD_DIR) $(BUILD_DIR)/$(EXECUTABLE)

debug: CFLAGS += -DDEBUG -g -pg -Wall
debug: LDFLAGS += -pg
debug: $(BUILD_DIR) $(BUILD_DIR)/$(EXECUTABLE)
