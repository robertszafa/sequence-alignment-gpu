#include <cuda.h>

#include "SequenceAlignment.hpp"

constexpr unsigned int MAX_THREADS_PER_BLOCK = 1024;


enum COLUMN_STATE { UNDEFINED, READY_FOR_0, READY_FOR_1, };

__device__ void yield(char* columnState, const int worker_id, const int col)
{
    columnState[col] = (worker_id == 0) ? COLUMN_STATE::READY_FOR_1 : COLUMN_STATE::READY_FOR_0;
}
__device__ void busy_wait(const char* columnState, const int worker_id, const int col)
{
    auto waitingFor = (worker_id == 0) ? COLUMN_STATE::READY_FOR_0 : COLUMN_STATE::READY_FOR_1;

    while (columnState[col] != waitingFor) continue;
}


__global__ void cuda_fillMatrixNW(const char *textBytes, const char *patternBytes,
                                  const int *scoreMatrix, const int alphabetSize,
                                  const int gapPenalty, const int startRow, const int endRow,
                                  const int numCols, const int workerId, char *columnState,
                                  char *M, int *finalScore)
{
    using SequenceAlignment::DIR;

    const int numRows = endRow - startRow;

    extern __shared__ int _shared[];
    int *_scoreMatrix = _shared;
    int *_thisScores = _shared + alphabetSize*alphabetSize;
    int *_prevScores = _thisScores + numRows;
    int *_prevPrevScores = _prevScores + numRows;

    const int tid = threadIdx.x;

    // Transfer score matrix to shared memory.
    for (int offset=0; offset < alphabetSize*alphabetSize; offset += blockDim.x)
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
        // Advance one diag.
        diag_size = min(diag_size+1, numRows);
        auto tmp = _prevScores;
        _prevScores = _thisScores;
        _thisScores = _prevPrevScores;
        _prevPrevScores = tmp;

        const int threadInRowIdx = i_text - tid;

        if (startRow > 0)
            busy_wait(columnState, workerId, i_text);

        if (tid == 0)
        {
            // First row.
            _thisScores[tid] = -(i_text * gapPenalty);
            M[i_text] = DIR::LEFT;
        }
        else if (tid == (diag_size-1) && i_text < numRows)
        {
            // First column.
            _thisScores[tid] = -(tid * gapPenalty);
            M[tid * numCols] = DIR::TOP;
        }
        else if (tid < diag_size)
        {
            const char textByte = textBytes[threadInRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            const int fromLeftScore = _prevScores[tid] - gapPenalty;
            const int fromTopScore = _prevScores[tid - 1] - gapPenalty;
            const int fromDiagScore = _prevPrevScores[tid - 1] + _scoreMatrix[scoreMatrixIdx];

            const int maxWithGap = max(fromLeftScore, fromTopScore);
            _thisScores[tid] = max(maxWithGap, fromDiagScore);

            const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
            M[tid*numCols + threadInRowIdx] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
        }

        if (tid == (diag_size-1) && diag_size == numRows)
            yield(columnState, workerId, threadInRowIdx);

        __syncthreads();
    }

    // Second half of matrix filling.
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Advance one diag.
        auto tmp = _prevScores;
        _prevScores = _thisScores;
        _thisScores = _prevPrevScores;
        _prevPrevScores = tmp;

        const int threadInRowIdx = numCols-1 - tid + i_pattern;

        if (tid >= i_pattern)
        {
            const char textByte = textBytes[threadInRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            const int fromLeftScore = _prevScores[tid] - gapPenalty;
            const int fromTopScore = _prevScores[tid - 1] - gapPenalty;
            const int fromDiagScore = _prevPrevScores[tid - 1] + _scoreMatrix[scoreMatrixIdx];

            const int maxWithGap = max(fromLeftScore, fromTopScore);
            _thisScores[tid] = max(maxWithGap, fromDiagScore);

            const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
            M[tid*numCols + threadInRowIdx] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
        }

        if (tid == (diag_size-1))
            yield(columnState, workerId, threadInRowIdx);

        __syncthreads();
    }

    if (tid == (numRows - 1))
        finalScore[0] = _thisScores[tid];
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
        delete [] M;
    };

    /** Allocate and transfer memory to device */
    if (cudaMalloc(&d_finalScore, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols) != cudaSuccess)
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
    for (int i=0; i < (numRows/MAX_THREADS_PER_BLOCK + 1); ++i)
    {
        const int numThreads = std::min(MAX_THREADS_PER_BLOCK, numRows - startRow);
        const int endRow = startRow + numThreads;

        const int workerId = (i%2 == 0) ? 0 : 1;
        cuda_fillMatrixNW<<<1, numThreads, sharedMemSize>>>(d_textBytes, d_patternBytes,
                                                            d_scoreMatrix, request.alphabetSize,
                                                            request.gapPenalty, startRow, endRow,
                                                            numCols, workerId, d_columnState,
                                                            d_M, d_finalScore);

        startRow = std::min(startRow + MAX_THREADS_PER_BLOCK, numRows-1);
    }

    // std::cout << "Num rows: " << numRows << "\n";
    // std::cout << "Num col: " << numCols << "\n";
    // std::cout << "Num bytes in shared: " << sharedMemSize << "\n";

    if (cudaMemcpy(&(response->score), d_finalScore, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess ||
        cudaMemcpy(M, d_M, numRows*numCols, cudaMemcpyDeviceToHost) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        return;
    }

    traceBack(M, numRows, numCols, request, response);

    freeMemory();

}
