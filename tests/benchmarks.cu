#define BENCHMARK
/// Number of times an experiment is ran. Then the average is taken.
#define NUM_REPEATS 10
#define NUM_WARMUPS 5

#include "../SequenceAlignment.hpp"
#include "old_alignSequenceGPU.cu"

#include <iostream>
#include <utility>
#include <chrono>

#include <stdlib.h>


void fillDummyRequest(SequenceAlignment::Request &request, const uint64_t numRows, const uint64_t numCols)
{
    request.sequenceType = SequenceAlignment::programArgs::PROTEIN;
    request.alignmentType = SequenceAlignment::programArgs::GLOBAL;
    request.alphabet = SequenceAlignment::PROTEIN_ALPHABET;
    request.alphabetSize = SequenceAlignment::NUM_PROTEIN_CHARS;
    request.gapPenalty = 5;
    request.textNumBytes = numCols - 1;
    request.patternNumBytes = numRows - 1;
    request.textBytes = new char[request.textNumBytes];
    request.patternBytes = new char[request.patternNumBytes];

    // Fill with random letters from alphabet, translated to indices.
    for (int i=0; i<request.textNumBytes; ++i)
        request.textBytes[i] = rand() % (request.alphabetSize-1);
    for (int i=0; i<request.patternNumBytes; ++i)
        request.patternBytes[i] = rand() % (request.alphabetSize-1);

    parseScoreMatrixFile(SequenceAlignment::DEFAULT_PROTEIN_SCORE_MATRIX_FILE,
                         request.alphabetSize, request.scoreMatrix);
}


void benchmarkFillMatrixNW(bool cpu, bool gpu, bool old_gpu)
{
    std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(1024, 1024*2),
        std::make_pair(1024, 1024*4),
        std::make_pair(1024, 1024*8),
        // std::make_pair(1024, 1024*16),
        // std::make_pair(4000, 16000),
        // std::make_pair(1000, 56000),
        // std::make_pair(14000, 16000),
        // std::make_pair(200000, 220000),
    };

    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        fillDummyRequest(request, numRows, numCols);
        char *M = new char[numRows * numCols];

        int totalTimeCPU = 0;
        if (cpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                fillMatrixNW(M, numRows, numCols, request);
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                auto begin = std::chrono::steady_clock::now();
                fillMatrixNW(M, numRows, numCols, request);
                auto end = std::chrono::steady_clock::now();
                totalTimeCPU += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
            }
            std::cout << "CPU = " << (totalTimeCPU/NUM_REPEATS) << " ms\n";
        }

        int totalTimeGPU = 0;
        if (gpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                SequenceAlignment::alignSequenceGPU(request, &response);
            for (int i=0; i<NUM_REPEATS; ++i)
                totalTimeGPU += SequenceAlignment::alignSequenceGPU(request, &response);
            std::cout << "GPU = " << (totalTimeGPU/NUM_REPEATS) << " ms\n";
        }

        int totalTimeGPU_horizontal = 0;
        if (old_gpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                old_alignSequenceGPU(request, &response);
            for (int i=0; i<NUM_REPEATS; ++i)
                totalTimeGPU_horizontal += old_alignSequenceGPU(request, &response);
            std::cout << "GPU (horizontal) = " << (totalTimeGPU_horizontal/NUM_REPEATS) << " ms\n";
        }

        std::cout << "GPU Speedup = " << (double(totalTimeCPU)/double(totalTimeGPU)) << "\n";

        delete [] M;
    }
}



int main(int argc, const char *argv[])
{
    benchmarkFillMatrixNW(true, true, false);

    return 0;
}

