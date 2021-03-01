NVCC=nvcc
CUDA_FLAGS=-std=c++14 -m64 -arch=sm_30 -O3

BIN=alignSequence
TEST_BIN=test
BENCHMARK_BIN=benchmark

all : $(BIN) $(TEST_BIN) $(BENCHMARK_BIN)

$(BIN) : mainDriver.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CUDA_FLAGS) $(FLAGS) mainDriver.cu -o $(BIN)

$(TEST_BIN) : tests/tests.cu utilities.cpp alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp
		$(NVCC) $(CUDA_FLAGS) tests/tests.cu -o $(TEST_BIN)

$(BENCHMARK_BIN) : alignSequenceCPU.cpp alignSequenceGPU.cu SequenceAlignment.hpp tests/benchmarks.cu tests/old_alignSequenceGPU.cu
		$(NVCC) $(CUDA_FLAGS) tests/benchmarks.cu -o $(BENCHMARK_BIN)

clean :
		rm -f *.o $(BIN) $(TEST_BIN) $(BENCHMARK_BIN)
