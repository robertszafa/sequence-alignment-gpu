# nvcc --std=c++11 alignSequence.cu -o alignSequence
NVCC=nvcc
CXX=clang++
CXXFLAGS=-std=c++11

BIN=alignSequence
TEST_BIN=runTests

all : alignSequence test

alignSequence : alignSequence.cu
		$(NVCC) $(CXXFLAGS) alignSequence.cu -o $(BIN)

test : test/tests.cpp
		$(NVCC) $(CXXFLAGS) test/tests.cpp -o $(TEST_BIN)


clean :
		rm $(BIN)
