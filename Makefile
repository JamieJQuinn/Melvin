CC=g++
CFLAGS=-c -I$(INCLUDE_DIR) --std=c++14
LDFLAGS=
SRC_DIR=src
BUILD_DIR=build
INCLUDE_DIR=include

#SOURCES=src/main.cpp src/sim.cpp src/numerical_methods.cpp src/thomas_algorithm.cpp src/constants.cpp src/variable.cpp src/double_diffusive_sim.cpp
SOURCES=$(wildcard $(SRC_DIR)/*.cpp)
OBJECTS=$(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
EXECUTABLE=exe

GPU_SOURCES=$(wildcard $(SRC_DIR)/*.cu)
GPU_OBJECTS=$(patsubst $(SRC_DIR)/%.cu,$(BUILD_DIR)/%.o,$(GPU_SOURCES))

TEST_DIR=unit_test
TEST_SOURCES=$(wildcard $(TEST_DIR)/*.cpp)
TEST_OBJECTS=$(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(TEST_SOURCES))
GPU_TEST_SOURCES=$(wildcard $(TEST_DIR)/*.cu)
GPU_TEST_OBJECTS=$(patsubst $(TEST_DIR)/%.cu,$(BUILD_DIR)/%.o,$(GPU_TEST_SOURCES))

TEST_EXECUTABLE=test_exe
GPU_TEST_EXECUTABLE=test_exe_gpu

.PHONY: all
all: release

$(BUILD_DIR)/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	nvcc $(CFLAGS) $< -o $@

$(BUILD_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)

.PHONY: full-test
full-test:
	test/test_all_linear.sh
	test/test_all_nonlinear.sh

.PHONY: profile
profile: CFLAGS += -pg
profile: LDFLAGS += -pg
profile: release

.PHONY: release
release: CFLAGS += -DNDEBUG -O2 -fopenmp
release: LDFLAGS += -fopenmp
release: $(BUILD_DIR) $(BUILD_DIR)/$(EXECUTABLE)

.PHONY: debug
debug: CFLAGS += -DDEBUG -g -pg -Wall
debug: LDFLAGS += -pg
debug: $(BUILD_DIR) $(BUILD_DIR)/$(EXECUTABLE)

.PHONY: gpu
gpu: CFLAGS += -DCUDA
gpu: OBJECTS += $(GPU_OBJECTS)
gpu: CC = nvcc
gpu: $(BUILD_DIR) $(GPU_OBJECTS) $(BUILD_DIR)/$(EXECUTABLE)

.PHONY: gpu-test
gpu-test: CFLAGS += -DCUDA -pg
gpu-test: LDFLAGS += -pg
gpu-test: CC = nvcc
gpu-test: $(BUILD_DIR) $(BUILD_DIR)/$(GPU_TEST_EXECUTABLE)
	python3 tools/make_initial_conditions.py --output $(BUILD_DIR)/ICn1nZ101nN51 --salt_fingering --n_modes 64 --n_gridpoints 128 --modes $(shell seq 1 63)
	cd $(BUILD_DIR); ../test/print_test_constants.sh
	cd $(BUILD_DIR); ./$(GPU_TEST_EXECUTABLE)

.PHONY: test
test: CFLAGS += -DNDEBUG -O2 -fopenmp -pg
test: LDFLAGS += -fopenmp -pg
test: all $(BUILD_DIR)/$(TEST_EXECUTABLE)
	cd $(BUILD_DIR); ./$(TEST_EXECUTABLE)

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CC) $(CFLAGS) $< -o $@

$(BUILD_DIR)/%.o: $(TEST_DIR)/%.cu
	nvcc $(CFLAGS) $< -o $@

$(BUILD_DIR)/$(GPU_TEST_EXECUTABLE): $(OBJECTS) $(GPU_OBJECTS) $(GPU_TEST_OBJECTS) $(TEST_OBJECTS)
	$(CC) $(LDFLAGS) $(filter-out $(BUILD_DIR)/main.o, $(OBJECTS)) $(GPU_OBJECTS) $(GPU_TEST_OBJECTS) $(TEST_OBJECTS) -o $@

$(BUILD_DIR)/$(TEST_EXECUTABLE): $(OBJECTS) $(TEST_OBJECTS)
	$(CC) $(LDFLAGS) $(filter-out $(BUILD_DIR)/main.o, $(OBJECTS)) $(TEST_OBJECTS) -o $@
