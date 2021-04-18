#define BENCHMARK

/// Number of experiment repetitions.
#define NUM_REPEATS 5
/// Number of discarderded runs for warm up, bring data in mem, etc.
#define NUM_WARMUPS 5

#include <iostream>
#include <utility>
// #include <chrono>
// chrono is not available on barkla.
#include <omp.h>
#include <stdlib.h>

#include "../SequenceAlignment.hpp"
#include "old_alignSequenceGPU.cu"


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
        std::make_pair(1024, 1024),
        std::make_pair(1024*4, 1024*4),
        std::make_pair(1024*8, 1024*8),
        std::make_pair(1024*16, 1024*16),
        std::make_pair(1024*32, 1024*32),
        std::make_pair(1024*64, 1024*64),
    };

    std::cout << "\nGlobal alignment benchmark:\n";
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
                auto begin = omp_get_wtime();
                fillMatrixNW(M, numRows, numCols, request);
                auto end = omp_get_wtime();
                totalTimeCPU += (end - begin) * 1000.0;
            }
            std::cout << "CPU = " << (totalTimeCPU/NUM_REPEATS) << " ms\n"
                      << "MCUPS: " << int((numRows * numCols) / ((totalTimeCPU/NUM_REPEATS) * 1000.0)) << "\n\n"; // /1000 for seconds * 1000000 for MCUPS
        }

        int totalTimeGPU = 0;
        if (gpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                SequenceAlignment::alignSequenceGPU(request, &response);
            for (int i=0; i<NUM_REPEATS; ++i)
                totalTimeGPU += SequenceAlignment::alignSequenceGPU(request, &response);
            std::cout << "GPU = " << (totalTimeGPU/NUM_REPEATS) << " ms\n"
                      << "MCUPS: " << int((numRows * numCols) / ((totalTimeGPU/NUM_REPEATS) * 1000.0)) << "\n\n";
        }

        int totalTimeGPU_horizontal = 0;
        if (old_gpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                old_alignSequenceGPU(request, &response);
            for (int i=0; i<NUM_REPEATS; ++i)
                totalTimeGPU_horizontal += old_alignSequenceGPU(request, &response);
            std::cout << "GPU (horizontal) = " << (totalTimeGPU_horizontal/NUM_REPEATS) << " ms\n"
                      << "MCUPS: " << int((numRows * numCols) / ((totalTimeGPU_horizontal/NUM_REPEATS) * 1000.0)) << "\n\n";
        }

        std::cout << "GPU Speedup = " << (double(totalTimeCPU)/double(totalTimeGPU)) << "\n";

        delete [] M;
    }
}


void benchmarkFillMatrixSW(bool cpu, bool gpu)
{
    std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(512*2, 1024*2),
        std::make_pair(512*4, 1024*4),
        std::make_pair(512*8, 1024*8),
        std::make_pair(512*16, 1024*16),
        std::make_pair(512*32, 1024*32),
        std::make_pair(512*64, 1024*64),
    };

    std::cout << "\nLocal alignment benchmark:\n";

    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        fillDummyRequest(request, numRows, numCols);
        request.alignmentType = SequenceAlignment::programArgs::LOCAL;
        char *M = new char[numRows * numCols];

        int totalTimeCPU = 0;
        if (cpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                fillMatrixSW(M, numRows, numCols, request);
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                auto begin = omp_get_wtime();
                fillMatrixSW(M, numRows, numCols, request);
                auto end = omp_get_wtime();
                totalTimeCPU += (end - begin) * 1000.0;
            }
            std::cout << "CPU = " << (totalTimeCPU/NUM_REPEATS) << " ms\n"
                      << "MCUPS: " << int((numRows * numCols) / ((totalTimeCPU/NUM_REPEATS) * 1000.0)) << "\n\n";
        }

        int totalTimeGPU = 0;
        if (gpu)
        {
            for (int i=0; i<NUM_WARMUPS; ++i)
                SequenceAlignment::alignSequenceGPU(request, &response);
            for (int i=0; i<NUM_REPEATS; ++i)
                totalTimeGPU += SequenceAlignment::alignSequenceGPU(request, &response);
            std::cout << "GPU = " << (totalTimeGPU/NUM_REPEATS) << " ms\n"
                      << "MCUPS: " << int((numRows * numCols) / ((totalTimeGPU/NUM_REPEATS) * 1000.0)) << "\n\n";
        }

        std::cout << "GPU Speedup = " << (double(totalTimeCPU)/double(totalTimeGPU)) << "\n";

        delete [] M;
    }
}


void benchmarkBatch (bool cpu, bool gpu, bool isGlobal, int nBatches)
{
    std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(512, 512),
        std::make_pair(1024, 1024),
        std::make_pair(2048, 2048),
    };

    auto alignType = isGlobal ? "Global" : "Local";
    std::cout << "\n" << alignType << " alignment batch (" << nBatches << "x) benchmark:\n";
    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request requests[nBatches];
        SequenceAlignment::Response responsesCPU[nBatches];
        SequenceAlignment::Response responsesGPU[nBatches];
        for (int i=0; i < nBatches; ++i)
        {
            fillDummyRequest(requests[i], numRows, numCols);
            requests[i].alignmentType = isGlobal ? SequenceAlignment::programArgs::GLOBAL
                                                 : SequenceAlignment::programArgs::LOCAL;
        }

        int totalTimeCPU = 0;
        if (cpu)
        {
            auto begin = omp_get_wtime();
            for (int i=0; i<nBatches; ++i)
                alignSequenceCPU(requests[i], &responsesCPU[i]);

            auto end = omp_get_wtime();
            totalTimeCPU += (end - begin) * 1000.0;
            std::cout << "CPU = " << totalTimeCPU << " ms\n"
                      << "MCUPS: " << int((numRows * numCols * nBatches) / (totalTimeCPU * 1000.0)) << "\n\n"; // /1000 for seconds * 1000000 for MCUPS
        }

        int totalTimeGPU = 0;
        if (gpu)
        {
            auto begin = omp_get_wtime();
            for (int i=0; i<nBatches; ++i)
                SequenceAlignment::alignSequenceGPU(requests[i], &responsesGPU[i]);

            auto end = omp_get_wtime();
            totalTimeGPU += (end - begin) * 1000.0;
            std::cout << "GPU = " << totalTimeGPU << " ms\n"
                      << "MCUPS: " << int((numRows * numCols * nBatches) / (totalTimeGPU * 1000.0)) << "\n\n"; // /1000 for seconds * 1000000 for MCUPS
        }

        std::cout << "GPU Speedup = " << (double(totalTimeCPU)/double(totalTimeGPU)) << "\n";
    }
}


int main(int argc, const char *argv[])
{
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    std::cout << "Benchmark on GPU: " << deviceProp.name << "\n";

    // Don't time CPU on barkla.

    benchmarkFillMatrixNW(false, true, false);

    benchmarkFillMatrixSW(false, true);

    // NOTE:
    // Undefine BENCHMARK macro at top of file to use benchmarkBatch.
    // In the batch benchmark, we measure the execution time of the whole system
    // when multiple sequences are aligned in quick succession. The BENCHMARK macro declares that
    // the time of just one run should be mesured.
    // benchmarkBatch (false, true, true, 1000);
    // benchmarkBatch (false, true, false, 1000);

    return 0;
}

