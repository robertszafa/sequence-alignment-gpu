// Comment out to use time/latency benchmark, instead of throughput benchmark.
#define BENCHMARK

/// Best-out-of-N time is picked when measuring latency & throughput.
#define NUM_REPEATS 5

#include <iostream>
#include <utility>

#include <sys/time.h>
#include <stdlib.h>

#include "../SequenceAlignment.hpp"
#include "old_alignSequenceGPU.cu"

using SequenceAlignment::programArgs;


/// Create a Request with the specified numRows, numCols. Fill the sequences out with
/// random data from the alphabet.
void fillDummyRequest(SequenceAlignment::Request &request, const uint64_t numRows,
                      const uint64_t numCols, programArgs alignType=programArgs::GLOBAL)
{
    request.sequenceType = SequenceAlignment::programArgs::PROTEIN;
    request.alignmentType = alignType;
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

/// Return time interval between 2 Unix timestamps in microseconds.
uint64_t getTime(timeval &t1, timeval &t2)
{
    timeval res;
    timersub(&t2, &t1, &res);

    // Microseconds is 1e-6th of a second.
    return 1000000 * uint64_t(res.tv_sec) + uint64_t(res.tv_usec);
}

/// Show the speedup from having the GPU wavefront on the diagonal axis of the M matrix,
/// compared to on the horizontal axis.
void benchmarkDiagonalVsHorizontalGPU()
{
    std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(512, 512),
        std::make_pair(1024, 1024),
        std::make_pair(1024*2, 1024*2),
    };

    std::cout << "\nFill Matrix NW diagonal vs horizontal wavefront benchmark:\n";
    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        fillDummyRequest(request, numRows, numCols);
        char *M = new char[numRows * numCols];

        uint64_t gpuTime = UINT64_MAX;
        uint64_t gpuTime_horizontal = UINT64_MAX;

        for (int i=0; i<NUM_REPEATS; ++i)
            gpuTime = std::min(gpuTime, SequenceAlignment::alignSequenceGPU(request, &response));

        gpuTime = std::max(uint64_t(1), gpuTime);
        std::cout << "GPU (diagonal) = " << (gpuTime / 1000) << " ms\n"
                  << "MCUPS: " << ((numRows * numCols) / gpuTime) << "\n\n";

        for (int i=0; i<NUM_REPEATS; ++i)
            gpuTime_horizontal = std::min(gpuTime_horizontal, old_alignSequenceGPU(request, &response));

        gpuTime_horizontal = std::max(uint64_t(1), gpuTime_horizontal);
        std::cout << "GPU (horizontal) = " << (gpuTime_horizontal / 1000) << " ms\n"
                    << "MCUPS: " << ((numRows * numCols) / gpuTime_horizontal) << "\n\n";

        std::cout << "GPU diagonal Speedup = " << (double(gpuTime_horizontal) / double(gpuTime)) << "\n";

        delete [] M;
    }
}


/// Measure how many Millions of Cell Updates Per Second (MCUPS) do the cpu & gpu algorithms achieve.
void benchmarkFillMatrixThroughput(const bool cpu, const bool gpu, const programArgs alignType)
{
    const std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizesNW =
    {
        std::make_pair(256, 256),
        std::make_pair(512, 512),
        std::make_pair(1024, 1024),
        std::make_pair(1024*2, 1024*2),
        std::make_pair(1024*4, 1024*4),
        std::make_pair(1024*8, 1024*8),
        std::make_pair(1024*16, 1024*16),
        std::make_pair(1024*32, 1024*32),
        std::make_pair(1024*64, 1024*64),
    };
    const std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizesSW =
    {
        std::make_pair(256, 1024*32),
        std::make_pair(512, 1024*32),
        std::make_pair(1024, 1024*32),
        std::make_pair(1024*2, 1024*32),
        std::make_pair(1024*4, 1024*32),
        std::make_pair(1024*8, 1024*32),
        std::make_pair(1024*16, 1024*32),
        std::make_pair(1024*32, 1024*32),
    };

    const auto benchmarkSizes = (alignType == programArgs::GLOBAL) ? benchmarkSizesNW : benchmarkSizesSW;
    std::string alignTypeStr = (alignType == programArgs::GLOBAL) ? "Global" : "Local";

    std::cout << "\n" << alignTypeStr << " alignment benchmark:\n";

    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        fillDummyRequest(request, numRows, numCols, alignType);
        char *M = new char[numRows * numCols];

        struct timeval t1, t2;

        uint64_t cpuTime = UINT64_MAX;
        if (cpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                gettimeofday(&t1, 0);

                if (alignType == programArgs::GLOBAL) fillMatrixNW(M, numRows, numCols, request);
                else fillMatrixSW(M, numRows, numCols, request);

                gettimeofday(&t2, 0);

                cpuTime = std::min(cpuTime, getTime(t1, t2));
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
            // For throughput measurements, we only measure the computation on the GPU
            // and data transfer back, not mem allocation or setup. To achieve this, the
            // BENCHMARK macro is used and the mesurements happen inside alignSequenceGPU.
            for (int i=0; i<NUM_REPEATS; ++i)
                gpuTime = std::min(gpuTime, SequenceAlignment::alignSequenceGPU(request, &response));

            gpuTime = std::max(uint64_t(1), gpuTime);
            std::cout << "GPU = " << (gpuTime / 1000) << " ms\n"
                      << "MCUPS: " << ((numRows * numCols) / gpuTime) << "\n\n";
        }

        if (gpu && cpu)
            std::cout << "GPU Speedup = " << (double(cpuTime) / double(gpuTime)) << "\n";

        delete [] M;
    }
}


/// Measure the time to process nBatches requests in sequence (end-to-end).
void benchmarkEndToEndLatency (const bool cpu, const bool gpu, const programArgs alignType, const uint64_t nBatches)
{
    const std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizesNW =
    {
        std::make_pair(256, 256),
        std::make_pair(512, 512),
        std::make_pair(1024, 1024),
        std::make_pair(1024*4, 1024*4),
        std::make_pair(1024*8, 1024*8),
        std::make_pair(1024*16, 1024*16),
        std::make_pair(1024*32, 1024*32),
        std::make_pair(1024*64, 1024*64),
    };
    const std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizesSW =
    {
        std::make_pair(256, 1024*32),
        std::make_pair(512, 1024*32),
        std::make_pair(1024, 1024*32),
        std::make_pair(1024*2, 1024*32),
        std::make_pair(1024*4, 1024*32),
        std::make_pair(1024*8, 1024*32),
        std::make_pair(1024*16, 1024*32),
        std::make_pair(1024*32, 1024*32),
    };

    const auto benchmarkSizes = (alignType == programArgs::GLOBAL) ? benchmarkSizesNW : benchmarkSizesSW;
    std::string alignTypeStr = (alignType == programArgs::GLOBAL) ? "Global" : "Local";

    std::cout << "\n" << alignTypeStr << " alignment latency (end-to-end) benchmark:\n";

    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request request;
        SequenceAlignment::Response response;
        fillDummyRequest(request, numRows, numCols, alignType);

        struct timeval t1, t2;

        uint64_t cpuTime = UINT64_MAX;
        if (cpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                gettimeofday(&t1, 0);
                alignSequenceCPU(request, &response);
                gettimeofday(&t2, 0);

                cpuTime = std::min(cpuTime, getTime(t1, t2));
            }

            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "CPU = " << (cpuTime / 1000) << " ms\n";
        }

        uint64_t gpuTime = UINT64_MAX;
        if (gpu)
        {
            for (int i=0; i<NUM_REPEATS; ++i)
            {
                gettimeofday(&t1, 0);
                alignSequenceGPU(request, &response);
                gettimeofday(&t2, 0);

                gpuTime = std::min(gpuTime, getTime(t1, t2));
            }

            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "GPU = " << (gpuTime / 1000) << " ms\n";
        }

        if (gpu && cpu)
            std::cout << "GPU Speedup = " << (double(cpuTime)/double(gpuTime)) << "\n";
    }
}

/// Measure the time to process nBatches requests in sequence (end-to-end).
void benchmarkEndToEndBatch (const bool cpu, const bool gpu, const programArgs alignType, const uint64_t nBatches)
{
    // In this benchmark, we fix the sequence size and measure how the GPU speedup
    // behaves as the batch size grows.
    const std::vector<std::pair<uint64_t, uint64_t>> benchmarkSizes =
    {
        std::make_pair(1024*8, 1024*8),
    };

    std::string alignTypeStr = (alignType == programArgs::GLOBAL) ? "Global" : "Local";

    std::cout << "\n" << alignTypeStr << " alignment batch (" << nBatches << "x) benchmark:\n";

    for (const auto &sizePair : benchmarkSizes)
    {
        uint64_t numRows = sizePair.first;
        uint64_t numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

        SequenceAlignment::Request requests[nBatches];
        SequenceAlignment::Response responsesCPU[nBatches];
        SequenceAlignment::Response responsesGPU[nBatches];
        for (int i=0; i < nBatches; ++i)
            fillDummyRequest(requests[i], numRows, numCols, alignType);

        struct timeval t1, t2;

        uint64_t cpuTime = 0;
        if (cpu)
        {
            // warmup
            alignSequenceCPU(requests[0], &responsesCPU[0]);

            gettimeofday(&t1, 0);
            for (int i=0; i<nBatches; ++i)
                alignSequenceCPU(requests[i], &responsesCPU[i]);
            gettimeofday(&t2, 0);

            cpuTime = getTime(t1, t2);
            // we measure in microseconds (1e6) - no more conversion needed for MCUPS
            std::cout << "CPU = " << (cpuTime / 1000) << " ms\n";
        }

        uint64_t gpuTime = 0;
        if (gpu)
        {
            // warmup
            alignSequenceGPU(requests[0], &responsesGPU[0]);

            gettimeofday(&t1, 0);
            for (int i=0; i<nBatches; ++i)
                alignSequenceGPU(requests[i], &responsesGPU[i]);
            gettimeofday(&t2, 0);

            gpuTime = getTime(t1, t2);
            std::cout << "GPU = " << (gpuTime / 1000) << " ms\n";
        }

        if (gpu && cpu)
            std::cout << "GPU Speedup = " << (double(cpuTime)/double(gpuTime)) << "\n";
    }
}



int main(int argc, const char *argv[])
{
    // benchmarkDiagonalVsHorizontalGPU();

    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    std::cout << "Benchmark on GPU: " << deviceProp.name << "\n";

    // For these, make sure BENCHMARK is defined.
    #ifdef BENCHMARK
        benchmarkFillMatrixThroughput(true, true, programArgs::GLOBAL);
        benchmarkFillMatrixThroughput(true, true, programArgs::LOCAL);
    #endif

    // For the following benchmarks, make sure the BENCHMARK macro (top of file) is not defined.
    #ifndef BENCHMARK
        // Measure latency of 1 alignment.
        benchmarkEndToEndLatency(true, true, programArgs::GLOBAL, 1);
        benchmarkEndToEndLatency(true, true, programArgs::LOCAL, 1);

        // Measure time to do N alignments.
        // In the batch we measure the execution time of the whole system
        // when multiple sequences are aligned in quick succession.
        benchmarkEndToEndBatch(true, true, programArgs::GLOBAL, 1);
        benchmarkEndToEndBatch(true, true, programArgs::GLOBAL, 2);
        benchmarkEndToEndBatch(true, true, programArgs::GLOBAL, 4);
        benchmarkEndToEndBatch(true, true, programArgs::GLOBAL, 8);
        benchmarkEndToEndBatch(true, true, programArgs::GLOBAL, 16);
        benchmarkEndToEndBatch(true, true, programArgs::GLOBAL, 32);
    #endif

    return 0;
}

