NVCC=nvcc
CXX=clang++
CXXFLAGS=-std=c++14

BIN=alignSequence
TEST_BIN=runTests

all : $(BIN) $(TEST_BIN)

$(BIN) : mainDriver.o alignSequenceCPU.o
		$(NVCC) $(CXXFLAGS) alignSequenceCPU.o mainDriver.o -o $(BIN)

$(TEST_BIN) : runTests.o alignSequenceCPU.o
		$(NVCC) $(CXXFLAGS) runTests.o alignSequenceCPU.o -o $(TEST_BIN)


mainDriver.o : mainDriver.cu SequenceAlignment.hpp utilities.hpp
	$(NVCC) $(CXXFLAGS) -c mainDriver.cu -o mainDriver.o

runTests.o : test/tests.cpp utilities.hpp SequenceAlignment.hpp alignSequenceCPU.cpp
	$(CXX) $(CXXFLAGS) -c test/tests.cpp -o runTests.o

alignSequenceCPU.o : SequenceAlignment.hpp alignSequenceCPU.cpp
	$(CXX) $(CXXFLAGS) -c alignSequenceCPU.cpp -o alignSequenceCPU.o

clean :
		rm *.o $(BIN) $(TEST_BIN)
