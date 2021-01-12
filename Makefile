NVCC=nvcc
CXX=clang++
CXXFLAGS=-std=c++14

BIN=alignSequence
TEST_BIN=runTests

all : $(BIN) $(TEST_BIN)

$(BIN) : mainDriver.cu utils.hpp SequenceAlignment.hpp alignSequenceCPU.cpp
		$(NVCC) $(CXXFLAGS) mainDriver.cu -o $(BIN)

$(TEST_BIN) : test/tests.cpp utils.hpp SequenceAlignment.hpp alignSequenceCPU.cpp
		$(NVCC) $(CXXFLAGS) test/tests.cpp -o $(TEST_BIN)

clean :
		rm *.o $(BIN) $(TEST_BIN)
