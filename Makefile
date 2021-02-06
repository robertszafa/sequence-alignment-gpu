NVCC=nvcc
CXX=clang++
CXXFLAGS=-std=c++14 -m64 -arch=sm_30 -O3

BIN=alignSequence
TEST=runTests
BENCHMARK=benchmark

all : $(BIN) $(TEST) $(BENCHMARK)

$(BIN) : mainDriver.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CXXFLAGS) $(FLAGS) mainDriver.cu -o $(BIN)

$(TEST) : test/tests.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CXXFLAGS) test/tests.cu -o $(TEST)

$(BENCHMARK) : alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp benchmark.cu
		$(NVCC) $(CXXFLAGS) benchmark.cu -o $(BENCHMARK)

clean :
		rm *.o $(BIN) $(TEST_BIN) $(BENCHMARK)
