#include "../SequenceAlignment.hpp"

#include <iostream>
#include <utility>
#include <chrono>

/// Number of times an experiment is ran. Then the average is taken.
#define NUM_REPEATS 20


/// Initial version of the Needleman-Wunsch DP matrix filling algorithm.
/// This version exploits parallelism across a single row of of the matrix:
///     1. Check fromDiagScore in parallel.
///     2. Check fromTopScore in parallel.
///     3. Check fromLeftScore in sequence from left to right.
///
/// This version is slower than the version which exploits parallism across the diagonal
/// of the matrix, however, it is included here in the benchmarks to show that it was explored.
__global__ void cuda_fillMatrixNW_horizontal(const char *textBytes, const char *patternBytes,
                                             const int *scoreMatrix, const int alphabetSize,
                                             const int gapPenalty, const int numRows, const int numCols,
                                             char *M, int *finalScore)
{
    using SequenceAlignment::DIR;

    extern __shared__ int _shared[];
    int *_thisScores = _shared;
    int *_prevScores = _shared + numCols;

    const int tid = threadIdx.x;

    // Each thread copies one text letter.
    const char textByte = (tid > 0) ? textBytes[tid - 1] : alphabetSize;
    // Init first row.
    _thisScores[tid] = tid * -gapPenalty;
    M[tid] = DIR::LEFT;

    __syncthreads();

    // Dynamic programming loop.
    auto thisRowM = M + numCols;
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Advance one row.
         int *tmp = _thisScores;
        _thisScores = _prevScores;
        _prevScores = tmp;

        if (tid == 0)
        {
            _thisScores[tid] = -(i_pattern * gapPenalty);
            thisRowM[0] = DIR::TOP;
            continue;
        }

        const char patternByte = patternBytes[i_pattern - 1];
        const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

        // We are accessing the previous row - wait for all columns to finish.
        const int fromTopScore = _prevScores[tid] - gapPenalty;
        const int fromDiagScore = _prevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

        const bool isDiagGreaterThanTop = (fromDiagScore > fromTopScore);
        const int maxFromPrev = isDiagGreaterThanTop ? fromDiagScore : fromTopScore;
        const auto tmpDir = isDiagGreaterThanTop ? DIR::DIAG : DIR::TOP;

        for (int i_text = 1; i_text < numCols; ++i_text)
        {
            // We are accessing the previous column within a row.
            if (tid == i_text)
            {
                const int fromLeftScore = _thisScores[tid - 1] - gapPenalty;
                const bool isPrevGreater = (maxFromPrev > fromLeftScore);

                _thisScores[tid] = isPrevGreater ? maxFromPrev : fromLeftScore;
                thisRowM[i_text] = isPrevGreater ? tmpDir : DIR::LEFT;
            }
            __syncthreads();
        }

        thisRowM += numCols;
    }

    if (tid == (numCols - 1))
        finalScore[0] = _thisScores[tid];
}

unsigned int wrapperCuda_fillMatrixNW(char *M, const unsigned int numRows, const unsigned int numCols,
                                     const SequenceAlignment::Request &request, bool diagonal = true)
{
    int *d_finalScore, *d_scoreMatrix;
    char *d_textBytes, *d_patternBytes, *d_M, *d_columnState;

    auto freeMemory = [&]()
    {
        cudaFree(d_finalScore);
        cudaFree(d_scoreMatrix);
        cudaFree(d_M);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        cudaFree(d_columnState);
    };

    /** Allocate and transfer memory to device */
    if (cudaMalloc(&d_finalScore, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols) != cudaSuccess)
    {
        freeMemory();
        return 0;
    }
    if (cudaMemcpy(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_scoreMatrix, request.scoreMatrix, sizeof(int) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice) != cudaSuccess)
    {
        freeMemory();
        return 0;
    }
    /** End Allocate and transfer memory */

    auto begin = std::chrono::steady_clock::now();

    if (diagonal)
    {
        const unsigned int sharedMemSize = 3 * numRows * sizeof(int) +
                                           request.alphabetSize * request.alphabetSize * sizeof(int);
        const auto workerId = 0, startRow = 0;
        cuda_fillMatrixNW<<<1, numRows, sharedMemSize>>>(d_textBytes, d_patternBytes,
                                                         d_scoreMatrix, request.alphabetSize,
                                                         request.gapPenalty, startRow, numRows,
                                                         numCols, workerId, d_columnState,
                                                         d_M, d_finalScore);
    }
    else
    {
        const unsigned int sharedMemSize = 2 * numCols * sizeof(int);
        cuda_fillMatrixNW_horizontal<<<1, numCols, sharedMemSize>>>(d_textBytes, d_patternBytes,
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

    freeMemory();

    return (unsigned int) duration;
}




void benchmarkCudaFillMatrixNW_diagonal_vs_horizontal()
{
    std::vector<std::pair<int, int>> benchmarkSizes =
    {
        std::make_pair(1024, 1024),
        std::make_pair(1024, 1024*2),
        std::make_pair(1024, 1024*4),
        std::make_pair(1024, 1024*8),
        // std::make_pair(1024, 1024*16),
        // std::make_pair(1024, 1024*32),
    };

    for (const auto &sizePair : benchmarkSizes)
    {
        auto numRows = sizePair.first;
        auto numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

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
        char *M = &(std::vector<char>(numRows * numCols))[0]; // Automatically clean-up M.

        int totalTimeGPU_diagonal = 0;
        for (int i=0; i<NUM_REPEATS; ++i)
            totalTimeGPU_diagonal += wrapperCuda_fillMatrixNW(M, numRows, numCols, request);

        int totalTimeGPU_horizontal = 0;
        for (int i=0; i<NUM_REPEATS; ++i)
            totalTimeGPU_horizontal += wrapperCuda_fillMatrixNW(M, numCols, numRows, request, false);

        std::cout << "GPU (diagonal) = " << (totalTimeGPU_diagonal/NUM_REPEATS) << " ms\n";
        std::cout << "GPU (horizontal) = " << (totalTimeGPU_horizontal/NUM_REPEATS) << " ms\n";
    }

}

void benchmarkFillMatrixNW_GPU_vs_CPU()
{
    std::vector<std::pair<int, int>> benchmarkSizes =
    {
        std::make_pair(1024*2, 1024*2),
    };

    for (const auto &sizePair : benchmarkSizes)
    {
        auto numRows = sizePair.first;
        auto numCols = sizePair.second;
        std::cout << "-----  " << numRows << " x " << numCols << "  -----\n";

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
        char *M = &(std::vector<char>(numRows * numCols))[0]; // Automatically clean-up M.

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
            totalTimeGPU += wrapperCuda_fillMatrixNW(M, numRows, numCols, request);

        std::cout << "CPU = " << (totalTimeCPU/NUM_REPEATS) << " ms\n";
        std::cout << "GPU = " << (totalTimeGPU/NUM_REPEATS) << " ms\n";
    }

}



int main(int argc, const char *argv[])
{
    benchmarkCudaFillMatrixNW_diagonal_vs_horizontal();

    // benchmarkFillMatrixNW_GPU_vs_CPU();


    return 0;
}

