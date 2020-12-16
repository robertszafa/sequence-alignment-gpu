# nvcc --std=c++11 alignSequence.cu -o alignSequence
NVCC=nvcc
NVCCFLAGS=-std=c++11
BIN=alignSequence

all : alignSequence

alignSequence : alignSequence.cu
		$(NVCC) $(NVCCFLAGS) alignSequence.cu -o $(BIN)

clean :
		rm $(BIN)
