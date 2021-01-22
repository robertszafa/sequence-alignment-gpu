#include "SequenceAlignment.hpp"

#include <cuda.h>


__global__ void alignSequenceGlobalCUDA(const char *textBytes, const uint64_t textNumBytes,
                                        const char *patternBytes, const uint64_t patternNumBytes,
                                        const char *alphabet, const int alphabetSize,
                                        const short *scoreMatrix, const short gap_penalty,
                                        int *thisScores, int *prevScores)
{
}


void SequenceAlignment::alignSequenceGlobalGPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{
    const unsigned int numCols = request.textNumBytes + 1;
    const unsigned int numRows = request.patternNumBytes + 1;

    int *d_thisScores, *d_prevScores;
    short *d_scoreMatrix;
    char *d_textBytes, *d_patternBytes, *d_alphabet;

    if (cudaMalloc(&d_thisScores, sizeof(int) * numCols) != cudaSuccess ||
        cudaMalloc(&d_prevScores, sizeof(int) * numCols) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(short) * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_alphabet, request.alphabetSize) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cudaFree(d_thisScores);
        cudaFree(d_prevScores);
        cudaFree(d_scoreMatrix);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return;
    }

    if (cudaMemcpy(d_textBytes, &request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_patternBytes, &request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_scoreMatrix, &request.scoreMatrix, (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_alphabet, &request.alphabet, request.alphabetSize, cudaMemcpyHostToDevice) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cudaFree(d_thisScores);
        cudaFree(d_prevScores);
        cudaFree(d_scoreMatrix);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return;
    }

    alignSequenceGlobalCUDA<<<1, numCols>>>(d_textBytes, request.textNumBytes,
                                            d_patternBytes, request.patternNumBytes,
                                            d_alphabet, request.alphabetSize,
                                            d_scoreMatrix, request.gapPenalty,
                                            d_thisScores, d_prevScores);

    if (cudaMemcpy(&(d_thisScores + numCols - 1), response->score, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
        std::cout << "Could not copy back to host memory" << std::endl;

    cudaFree(d_thisScores);
    cudaFree(d_prevScores);
    cudaFree(d_scoreMatrix);
    cudaFree(d_textBytes);
    cudaFree(d_patternBytes);

    std::cout << "# Score: \t" << response.score << "\n";

}
