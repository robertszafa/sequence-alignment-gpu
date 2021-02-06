NVCC=nvcc
CUDA_FLAGS=-std=c++14 -m64 -arch=sm_30 -O3

BIN=alignSequence
TEST_BIN=runTests
BENCHMARK_BIN=runBenchmarks

all : $(BIN) $(TEST_BIN) $(BENCHMARK_BIN)

$(BIN) : mainDriver.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CUDA_FLAGS) $(FLAGS) mainDriver.cu -o $(BIN)

$(TEST_BIN) : test/tests.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CUDA_FLAGS) test/tests.cu -o $(TEST_BIN)

$(BENCHMARK_BIN) : alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp test/benchmarks.cu
		$(NVCC) $(CUDA_FLAGS) test/benchmarks.cu -o $(BENCHMARK_BIN)

clean :
		rm -f *.o $(BIN) $(TEST_BIN) $(BENCHMARK_BIN)
