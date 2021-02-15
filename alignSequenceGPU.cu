#include <cuda.h>

#include "SequenceAlignment.hpp"

constexpr unsigned int MAX_THREADS_PER_BLOCK = 1024;

using SequenceAlignment::DIR;

struct columnState { int score; int kernelId; };

__device__ __forceinline__ void set_done(columnState* volatile colState, const int col,
                                         const int score, const int kernelId)
{
    colState[col] = {score, kernelId+1};
}

__device__ __forceinline__ void busy_wait(columnState* volatile colState, const int col,
                                          const int kernelId)
{
    while (colState[col].kernelId != kernelId) continue;
}

__device__ __forceinline__ void ping_pong_buffers(int *&first, int *&second, int *&third)
{
    auto tmpSecond = second;
    second = first;
    first = third;
    third = tmpSecond;
}

__device__ __forceinline__ int choose_direction(const int leftScore,
                                                const int topScore,
                                                const int diagScore,
                                                const int gapPenalty,
                                                const int letterScore,
                                                char *directionStore)
{
    const int fromLeftScore = leftScore - gapPenalty;
    const int fromTopScore = topScore - gapPenalty;
    const int fromDiagScore = diagScore + letterScore;

    const int maxWithGap = max(fromLeftScore, fromTopScore);
    const int maxOverall = max(maxWithGap, fromDiagScore);

    const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
    directionStore[0] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;

    return maxOverall;
}

__global__ void cuda_fillMatrixNW(const char *textBytes, const char *patternBytes,
                                  const int *scoreMatrix, const int alphabetSize, const int gapPenalty,
                                  const int startRow, const int endRow, const int numCols,
                                  const int kernelId, columnState* colState, char *M)
{
    const int numRows = blockDim.x;

    extern __shared__ int _shared[];
    int *_scoreMatrix = _shared;
    int *_thisScores = _shared + alphabetSize*alphabetSize;
    int *_prevScores = _thisScores + numRows;
    int *_prevPrevScores = _prevScores + numRows;

    const int tid = threadIdx.x;

    // Transfer score matrix to shared memory.
    for (int offset=0; offset < alphabetSize*alphabetSize; offset += numRows)
    {
        if ((offset + tid) < alphabetSize*alphabetSize)
            _scoreMatrix[offset + tid] = scoreMatrix[offset + tid];
    }

    const char patternByte = patternBytes[max(0, tid - 1 + startRow)];

    __syncthreads();

    // First half of matrix filling
    int diag_size = 0;
    for (int i_text = 0; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        diag_size = min(diag_size+1, numRows);
        const int threadInRowIdx = i_text - tid;

        if (startRow > 0)
            busy_wait(colState, i_text, kernelId);

        if (tid == 0 && startRow == 0)
        {
            // First row.
            _thisScores[tid] = -(i_text * gapPenalty);
            M[i_text] = DIR::LEFT;
        }
        else if (tid == (diag_size-1) && i_text < numRows)
        {
            // First column.
            _thisScores[tid] = -((tid + startRow) * gapPenalty);
            M[tid * numCols] = DIR::TOP;
        }
        else if (tid < diag_size)
        {
            const char textByte = textBytes[threadInRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            _thisScores[tid] = choose_direction(_prevScores[tid],
                                                _prevScores[tid - 1],
                                                _prevPrevScores[tid - 1],
                                                gapPenalty, _scoreMatrix[scoreMatrixIdx],
                                                (M + tid*numCols + threadInRowIdx));
        }
        __syncthreads();

        if (tid == endRow && tid < diag_size)
            set_done(colState, threadInRowIdx, _thisScores[tid], kernelId);
    }

    // Second half of matrix filling.
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        const int threadInRowIdx = numCols-1 - tid + i_pattern;

        if (tid >= i_pattern)
        {
            const char textByte = textBytes[threadInRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            _thisScores[tid] = choose_direction(_prevScores[tid],
                                                _prevScores[tid - 1],
                                                _prevPrevScores[tid - 1],
                                                gapPenalty, _scoreMatrix[scoreMatrixIdx],
                                                (M + tid*numCols + threadInRowIdx));

            if (tid == endRow)
                set_done(colState, threadInRowIdx, _thisScores[tid], kernelId);
        }

        __syncthreads();
    }

}


void SequenceAlignment::alignSequenceGlobalGPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{
    const unsigned int numCols = request.textNumBytes + 1;
    const unsigned int numRows = request.patternNumBytes + 1;
    char *M;

    /** Allocate host memory */
    try
    {
        M = new char[numRows * numCols];
        response->alignedTextBytes = new char[2 * request.textNumBytes];
        response->alignedPatternBytes = new char[2 * request.textNumBytes];
    }
    catch(const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        return;
    }
    /** End Allocate host memory */

    char *d_textBytes, *d_patternBytes, *d_M;
    int *d_scoreMatrix;
    columnState *d_columnState;

    auto freeMemory = [&]()
    {
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        cudaFree(d_M);
        cudaFree(d_scoreMatrix);
        cudaFree(d_columnState);
        delete [] M;
    };

    /** Allocate and transfer memory to device */
    if (cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols * sizeof(columnState)) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        freeMemory();
        return;
    }
    if (cudaMemcpy(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_scoreMatrix, request.scoreMatrix, sizeof(int) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        freeMemory();
        return;
    }
    /** End Allocate and transfer memory */

    const unsigned int sharedMemSize = 3 * numRows * sizeof(int) +
                                       request.alphabetSize * request.alphabetSize * sizeof(int);

    int startRow = 0;
    for (int i_kernel=0; i_kernel < (numRows/MAX_THREADS_PER_BLOCK + 1); ++i_kernel)
    {
        const int numThreads = std::min(MAX_THREADS_PER_BLOCK, numRows - startRow);
        const int endRow = startRow + numThreads - 1;

        cuda_fillMatrixNW<<<1, numThreads, sharedMemSize>>>(d_textBytes, d_patternBytes,
                                                            d_scoreMatrix, request.alphabetSize,
                                                            request.gapPenalty, startRow, endRow,
                                                            numCols, i_kernel, d_columnState, d_M);

        startRow = endRow + 1;
    }

    // std::cout << "Num rows: " << numRows << "\n";
    // std::cout << "Num col: " << numCols << "\n";
    // std::cout << "Num bytes in shared: " << sharedMemSize << "\n";

    cudaDeviceSynchronize();
    if (cudaMemcpy(&(response->score), &(d_columnState[numCols - 1].score), sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess ||
        cudaMemcpy(M, d_M, numRows*numCols, cudaMemcpyDeviceToHost) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        freeMemory();
        return;
    }

    traceBack(M, numRows, numCols, request, response);

    freeMemory();

}
