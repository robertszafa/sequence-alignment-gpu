NVCC=nvcc
CXX=clang++
CXXFLAGS=-std=c++14 -m64 #-O3

BIN=alignSequence
TEST_BIN=runTests

all : $(BIN) $(TEST_BIN)

$(BIN) : mainDriver.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CXXFLAGS) mainDriver.cu -o $(BIN)

$(TEST_BIN) : test/tests.cpp utilities.cpp alignSequenceCPU.cpp SequenceAlignment.hpp
		$(NVCC) $(CXXFLAGS) test/tests.cpp -o $(TEST_BIN)

clean :
		rm *.o $(BIN) $(TEST_BIN)
