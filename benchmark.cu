#include "SequenceAlignment.hpp"

#include <iostream>
#include <utility>

#include <chrono>


unsigned int wrapperGPU_fillMatrixNW(char *M, const unsigned int numRows, const unsigned int numCols,
                                     const SequenceAlignment::Request &request, bool diagonal = false)
{
    int *d_finalScore;
    short *d_scoreMatrix;
    char *d_textBytes, *d_patternBytes, *d_M;

    /** Allocate and transfer memory */
    if (cudaMalloc(&d_finalScore, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(short) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess)
    {
        cudaFree(d_finalScore);
        cudaFree(d_scoreMatrix);
        cudaFree(d_M);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return 0;
    }
    if (cudaMemcpy(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_scoreMatrix, request.scoreMatrix, sizeof(short) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice) != cudaSuccess)
    {
        cudaFree(d_finalScore);
        cudaFree(d_scoreMatrix);
        cudaFree(d_M);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return 0;
    }
    /** End Allocate and transfer memory */

    auto begin = std::chrono::steady_clock::now();

    if (diagonal)
    {
        const unsigned int sharedMemSize = 3 * numRows * sizeof(int);
        cuda_fillMatrixDiagonalNW<<<1, numRows, sharedMemSize>>>(d_textBytes, request.textNumBytes,
                                                            d_patternBytes, request.patternNumBytes,
                                                            d_scoreMatrix, request.alphabetSize,
                                                            request.gapPenalty, numRows, numCols,
                                                            d_M, d_finalScore);

    }
    else
    {
        const unsigned int sharedMemSize = 2 * numCols * sizeof(int);

        cuda_fillMatrixNW<<<1, numCols, sharedMemSize>>>(d_textBytes, request.textNumBytes,
                                                    d_patternBytes, request.patternNumBytes,
                                                    d_scoreMatrix, request.alphabetSize,
                                                    request.gapPenalty, numRows, numCols,
                                                    d_M, d_finalScore);
    }

    if (cudaMemcpy(M, d_M, numRows*numCols, cudaMemcpyDeviceToHost) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        return 0;
    }

    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

    cudaFree(d_finalScore);
    cudaFree(d_scoreMatrix);
    cudaFree(d_M);
    cudaFree(d_textBytes);
    cudaFree(d_patternBytes);

    return (unsigned int) duration;
}



int main(int argc, const char *argv[])
{
    const int NUM_REPEATS = 20;

    std::vector<std::pair<int, int>> benchmarkSizes =
    {
        std::make_pair(1024, 1024),
        std::make_pair(1024*2, 1024),
        std::make_pair(1024*4, 1024),
        // std::make_pair(1024, 1024*8),
        // std::make_pair(1024, 1024*16),
        // std::make_pair(1024, 1024*32),
        // std::make_pair(1024, 1024*64),
    };


    for (const auto &sizePair : benchmarkSizes)
    {
        auto numRows = sizePair.first;
        auto numCols = sizePair.second;

        SequenceAlignment::Request request;
        request.sequenceType = SequenceAlignment::programArgs::PROTEIN;
        request.alignmentType = SequenceAlignment::programArgs::GLOBAL;
        request.alphabet = SequenceAlignment::PROTEIN_ALPHABET;
        request.alphabetSize = SequenceAlignment::NUM_PROTEIN_CHARS;
        request.gapPenalty = 5;
        request.textNumBytes = numCols - 1;
        request.patternNumBytes = numRows - 1;
        request.textBytes = new char[request.textNumBytes];
        request.patternBytes = new char[request.patternNumBytes];
        std::fill_n(request.textBytes, request.textNumBytes, 'A');
        std::fill_n(request.patternBytes, request.patternNumBytes, 'A');
        char *M = &(std::vector<char>(numRows * numCols))[0]; // Automatically delete M.

        int totalTimeCPU = 0;
        for (int i=0; i<NUM_REPEATS; ++i)
        {
            auto begin = std::chrono::steady_clock::now();

            fillMatrixNW(M, numRows, numCols, request);
            auto end = std::chrono::steady_clock::now();
            totalTimeCPU += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
        }

        int totalTimeGPU = 0;
        for (int i=0; i<NUM_REPEATS; ++i)
        {
            totalTimeGPU += wrapperGPU_fillMatrixNW(M, numRows, numCols, request);
        }

        int totalTimeDiagonalGPU = 0;
        for (int i=0; i<NUM_REPEATS; ++i)
        {
            totalTimeDiagonalGPU += wrapperGPU_fillMatrixNW(M, numCols, numRows, request, true);
        }

        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";
        std::cout << "CPU = " << (totalTimeCPU/NUM_REPEATS) << " ms\n";
        std::cout << "GPU = " << (totalTimeGPU/NUM_REPEATS) << " ms\n";
        std::cout << "GPU (diagonal) = " << (totalTimeDiagonalGPU/NUM_REPEATS) << " ms\n";
    }

    return 0;
}

