#define BENCHMARK

/// Best-out-of-N time is picked
#define NUM_REPEATS 10

#include <iostream>
#include <utility>

#include <sys/time.h>
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
        std::make_pair(512, 512),
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

        struct timeval t1, t2, res;

        uint64_t cpuTime = UINT64_MAX;
        if (cpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                gettimeofday(&t1, 0);
                fillMatrixNW(M, numRows, numCols, request);
                gettimeofday(&t2, 0);

                timersub(&t2, &t1, &res);
                cpuTime = std::min(cpuTime, 1000000 * uint64_t(res.tv_sec) + uint64_t(res.tv_usec));
            }

            // Ensure we don't divide by 0.
            cpuTime = std::max(uint64_t(1), cpuTime);
            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "CPU = " << (cpuTime / 1000) << " ms\n" // display in ms
                      << "MCUPS: " << ((numRows * numCols) / cpuTime) << "\n\n";
        }

        uint64_t gpuTime = UINT64_MAX;
        if (gpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
                gpuTime = std::min(gpuTime, SequenceAlignment::alignSequenceGPU(request, &response));

            gpuTime = std::max(uint64_t(1), gpuTime);
            std::cout << "GPU = " << (gpuTime / 1000) << " ms\n"
                      << "MCUPS: " << ((numRows * numCols) / gpuTime) << "\n\n";
        }

        uint64_t gpuTime_horizontal = UINT64_MAX;
        if (old_gpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
                gpuTime_horizontal = std::min(gpuTime_horizontal, old_alignSequenceGPU(request, &response));

            gpuTime_horizontal = std::max(uint64_t(1), gpuTime_horizontal);
            std::cout << "GPU (horizontal) = " << (gpuTime_horizontal / 1000) << " ms\n"
                      << "MCUPS: " << ((numRows * numCols) / gpuTime_horizontal) << "\n\n";
        }

        std::cout << "GPU Speedup = " << (double(cpuTime) / double(gpuTime)) << "\n";

        delete [] M;
    }
}


void benchmarkFillMatrixSW(bool cpu, bool gpu)
{
    std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(512, 1024),
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

        struct timeval t1, t2, res;

        uint64_t cpuTime = UINT64_MAX;
        if (cpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                gettimeofday(&t1, 0);
                fillMatrixSW(M, numRows, numCols, request);
                gettimeofday(&t2, 0);

                timersub(&t2, &t1, &res);
                cpuTime = std::min(cpuTime, 1000000 * uint64_t(res.tv_sec) + uint64_t(res.tv_usec));
            }

            cpuTime = std::max(uint64_t(1), cpuTime);
            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "CPU = " << (cpuTime / 1000) << " ms\n" // display in ms
                      << "MCUPS: " << ((numRows * numCols) / cpuTime) << "\n\n";
        }

        uint64_t gpuTime = UINT64_MAX;
        if (gpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
                gpuTime = std::min(gpuTime, SequenceAlignment::alignSequenceGPU(request, &response));

            gpuTime = std::max(uint64_t(1), gpuTime);
            std::cout << "GPU = " << (gpuTime / 1000) << " ms\n"
                      << "MCUPS: " << ((numRows * numCols) / gpuTime) << "\n\n";
        }

        std::cout << "GPU Speedup = " << (double(cpuTime)/double(gpuTime)) << "\n";

        delete [] M;
    }
}


void benchmarkBatch (bool cpu, bool gpu, bool isGlobal, uint64_t nBatches)
{
    std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(512, 512),
        std::make_pair(1024, 1024),
        std::make_pair(1024*2, 1024*2),
        std::make_pair(1024*4, 1024*4),
        std::make_pair(1024*8, 1024*8),
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

        struct timeval t1, t2, res;

        uint64_t cpuTime = UINT64_MAX;
        if (cpu)
        {
            gettimeofday(&t1, 0);
            for (int i=0; i<nBatches; ++i)
                alignSequenceCPU(requests[i], &responsesCPU[i]);
            gettimeofday(&t2, 0);

            timersub(&t2, &t1, &res);
            cpuTime = std::min(cpuTime, 1000000 * uint64_t(res.tv_sec) + uint64_t(res.tv_usec));
            cpuTime = std::max(uint64_t(1), cpuTime);
            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "CPU = " << (cpuTime / 1000) << " ms\n" // display in ms
                      << "MCUPS: " << ((numRows * numCols * nBatches) / cpuTime) << "\n\n";
        }

        uint64_t gpuTime = UINT64_MAX;
        if (gpu)
        {
            gettimeofday(&t1, 0);
            for (int i=0; i<nBatches; ++i)
                alignSequenceGPU(requests[i], &responsesGPU[i]);
            gettimeofday(&t2, 0);

            timersub(&t2, &t1, &res);
            gpuTime = std::min(gpuTime, 1000000 * uint64_t(res.tv_sec) + uint64_t(res.tv_usec));
            gpuTime = std::max(uint64_t(1), gpuTime);
            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "GPU = " << (gpuTime / 1000) << " ms\n" // display in ms
                      << "MCUPS: " << ((numRows * numCols * nBatches) / gpuTime) << "\n\n";
        }

        std::cout << "GPU Speedup = " << (double(cpuTime)/double(gpuTime)) << "\n";
    }
}


int main(int argc, const char *argv[])
{
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    std::cout << "Benchmark on GPU: " << deviceProp.name << "\n";

    benchmarkFillMatrixNW(true, true, false);

    benchmarkFillMatrixSW(true, true);

    // NOTE:
    // Undefine BENCHMARK macro at top of file to use benchmarkBatch.
    // In the batch benchmark, we measure the execution time of the whole system
    // when multiple sequences are aligned in quick succession. The BENCHMARK macro declares that
    // the time of just one run should be mesured.
    // benchmarkBatch (true, true, true, 1000);
    // benchmarkBatch (true, true, false, 1000);

    return 0;
}

