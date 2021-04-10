#include <cuda.h>

#include "SequenceAlignment.hpp"

#define MAX_THREADS_PER_BLOCK 1024
#define MAX_CONCURRENT_KERNELS 32

using SequenceAlignment::DIRECTION;

struct columnState { int score; int kernelId; };


/// Given a column and a kernelId, set the kernelId of that column to the supplied kernelId.
__device__ __forceinline__ void set_done(columnState* volatile colState, const unsigned int col,
                                         const int score, const int kernelId)
{
    colState[col] = {score, kernelId+1};
}

/// Given a column and a kernelId, spin while the kernelId of that column
/// is not equal to the supplied kernelId.
/// Once finished waiting, return the score of the columnState.
__device__ __forceinline__ int get_prev_score(columnState* volatile colState, const unsigned int col,
                                              const int kernelId)
{
    volatile int currKernelId = colState[col].kernelId;
    while (currKernelId != kernelId)
        currKernelId = colState[col].kernelId;

    return colState[col].score;
}

/// Swap 3 pointers: third will point to old_second, second to old_first, first to old_third.
__device__ __forceinline__ void ping_pong_buffers(int *&first, int *&second, int *&third)
{
    auto tmpSecond = second;
    second = first;
    first = third;
    third = tmpSecond;
}

/// Given scores of North, West and North-West neigbours,
/// return the best score and direction.
__device__ __forceinline__ std::pair<int, DIRECTION> choose_direction_NW(const int leftScore,
                                                                   const int topScore,
                                                                   const int diagScore,
                                                                   const int gapPenalty,
                                                                   const int letterScore)
{
    const int fromLeftScore = leftScore - gapPenalty;
    const int fromTopScore = topScore - gapPenalty;
    const int fromDiagScore = diagScore + letterScore;

    const int maxWithGap = max(fromLeftScore, fromTopScore);
    const int bestScore = max(maxWithGap, fromDiagScore);

    const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIRECTION::TOP : DIRECTION::LEFT;
    const auto bestDir = (fromDiagScore > maxWithGap) ? DIRECTION::DIAG : dirWithGap;

    return {bestScore, bestDir};
}

__global__ void cuda_fillMatrixNW(const char* __restrict__ textBytes,
                                  const char* __restrict__ patternBytes,
                                  const int* __restrict__ scoreMatrix,
                                  const int alphabetSize, const int gapPenalty,
                                  const int startRow, const int endRow, const int numCols,
                                  const int kernelId,
                                  columnState* __restrict__ colState,
                                  char* __restrict__ M)
{
    const int numRows = blockDim.x;
    const int tid = threadIdx.x;

    // Set up shared memory pointers.
    extern __shared__ int _shared[];
    int *_scoreMatrix = _shared;
    int *_thisScores = _shared + alphabetSize*alphabetSize;
    int *_prevScores = _thisScores + numRows;
    int *_prevPrevScores = _prevScores + numRows;

    // Transfer score matrix to shared memory.
    for (int offset=0; offset < alphabetSize * alphabetSize; offset += numRows)
    {
        if (offset + tid < alphabetSize * alphabetSize)
            _scoreMatrix[offset + tid] = scoreMatrix[offset + tid];
    }

    // Each thread gets one row - one pattern letter (note that startRow >= 1).
    const char patternByte = patternBytes[tid + startRow - 1];
    M[tid * numCols] = DIRECTION::TOP;

    /* First half of matrix filling */
    int diagonalSize = 0;
    int fromLeft = 0, fromDiag = 0, fromTop = 0, fromPrevKernel = 0;
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        diagonalSize = min(diagonalSize + 1, numRows);
        const int idxInRow = i_text - tid;

        // Need to wait for previous kernel to finish with this column.
        if (tid == 0 && kernelId > 0)
            fromPrevKernel = get_prev_score(colState, i_text, kernelId);
        else if (tid == 0)
            fromPrevKernel = i_text * -gapPenalty;

        // Join previous iteration's threads and thread 0.
        __syncthreads();

        if (tid < diagonalSize)
        {
            fromLeft = (idxInRow == 1)                  // Am I in the first column to the left?
                       ? (tid + startRow) * -gapPenalty // If yes, then column to left is the empty char.
                       : _prevScores[tid];              // Otherwise, I'm somewhere in the middle.
            fromDiag = fromTop;
            fromTop = (tid == 0)                        // Am I the first row?
                      ? fromPrevKernel                  // If yes, then need to get score from previous kernel.
                      : _prevScores[tid - 1];           // Otherwise, I'm somewhere in the middle.

            const char textByte = textBytes[idxInRow - 1];
            const int scoreMatrixIdx = textByte * alphabetSize + patternByte;

            // Choose best direction, save score to shared and direction to global memory.
            auto scoreDirPair = choose_direction_NW(fromLeft, fromTop, fromDiag, gapPenalty,
                                                    _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            // If I'm I have finishes with a cell in the last row of this kernel,
            // then this column is done. Set done state and save score in global memory.
            if ((tid + startRow) == endRow)
                set_done(colState, idxInRow, scoreDirPair.first, kernelId);
        }
    }

    /* Second half of matrix filling */
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        const int idxInRow = (numCols - 1) + (i_pattern - tid);

        __syncthreads();

        if (tid >= i_pattern)
        {
            const char textByte = textBytes[idxInRow - 1];
            const int scoreMatrixIdx = textByte * alphabetSize + patternByte;

            // Choose best direction, save score to shared and direction to global memory.
            auto scoreDirPair = choose_direction_NW(_prevScores[tid], _prevScores[tid - 1], _prevPrevScores[tid - 1],
                                                    gapPenalty, _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            // If I'm I have finishes with a cell in the last row of this kernel,
            // then this column is done. Set done state and save score in global memory.
            if ((tid + startRow) == endRow)
                set_done(colState, idxInRow, scoreDirPair.first, kernelId);
        }
    }

}

/// Given scores of North, West and North-West neigbours,
/// return the best score and direction. If the best score is 0, then return DIRECTION::STOP.
__device__ __forceinline__ std::pair<int, DIRECTION> choose_direction_SW(const int leftScore,
                                                                   const int topScore,
                                                                   const int diagScore,
                                                                   const int gapPenalty,
                                                                   const int letterScore)
{
    const int fromLeftScore = leftScore - gapPenalty;
    const int fromTopScore = topScore - gapPenalty;
    const int fromDiagScore = diagScore + letterScore;

    const int maxWithGap = max(fromLeftScore, fromTopScore);
    const int bestScore = max(max(maxWithGap, fromDiagScore), 0);

    const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIRECTION::TOP : DIRECTION::LEFT;
    const auto bestDirNonZero = (fromDiagScore > maxWithGap) ? DIRECTION::DIAG : dirWithGap;
    const auto bestDir = bestScore > 0 ? bestDirNonZero : DIRECTION::STOP;

    return {bestScore, bestDir};
}

/// Given an array values, compute the maximum of of the N first items and store it in values[0].
__device__ __forceinline__ void max_reduce(int *values, const int N)
{
    const int tid = threadIdx.x;

    // Tree reduction with log2 levels.
    for (int pow2 = 1; pow2 < N; pow2 <<= 1)
    {
        // Only threads at pow2 indexes do work.
        if ((tid & pow2) == 0 && (tid + pow2) < N)
            values[tid] = max(values[tid], values[tid + pow2]);

        __syncthreads();
    }
}

__global__ void cuda_fillMatrixSW(const char* __restrict__ textBytes,
                                  const char* __restrict__ patternBytes,
                                  const int* __restrict__ scoreMatrix,
                                  const int alphabetSize, const int gapPenalty,
                                  const int startRow, const int endRow, const int numCols,
                                  const int kernelId,
                                  columnState* __restrict__ colState,
                                  int* __restrict__ p_bestScore, int* __restrict__ p_bestScoreIdx,
                                  char* __restrict__ M)
{
    const int numRows = blockDim.x;
    const int tid = threadIdx.x;

    // Set up shared memory pointers.
    extern __shared__ int _shared[];
    int *_scoreMatrix = _shared;
    int *_thisScores = _shared + alphabetSize*alphabetSize;
    int *_prevScores = _thisScores + numRows;
    int *_prevPrevScores = _prevScores + numRows;

    // Transfer score matrix to shared memory.
    for (int offset=0; offset < alphabetSize * alphabetSize; offset += numRows)
    {
        if (offset + tid < alphabetSize * alphabetSize)
            _scoreMatrix[offset + tid] = scoreMatrix[offset + tid];
    }

    // Each thread gets one row - one pattern letter (note that startRow >= 1).
    const char patternByte = patternBytes[tid + startRow - 1];
    M[tid * numCols] = DIRECTION::STOP;

    /* First half of matrix filling */
    int diagonalSize = 0;
    int fromLeft = 0, fromDiag = 0, fromTop = 0, fromPrevKernel = 0;
    int thisBestScore = 0, thisBestScoreIdx = 0;
    // The fromDiagonal score at iteration_i, becomes the fromTop score at itreation_i+1.
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        diagonalSize = min(diagonalSize + 1, numRows);
        const int idxInRow = i_text - tid;

        // Need to wait for previous kernel to finish with this column.
        if (tid == 0 && kernelId > 0)
            fromPrevKernel = get_prev_score(colState, i_text, kernelId);
        else if (tid == 0)
            fromPrevKernel = 0;

        // Join previous iteration's threads and thread 0.
        __syncthreads();

        if (tid < diagonalSize)
        {
            fromLeft = (idxInRow == 1)                  // Am I in the first column to the left?
                       ? 0                              // If yes, then column to left has score 0.
                       : _prevScores[tid];              // Otherwise, I'm somewhere in the middle.
            fromDiag = fromTop;
            fromTop = (tid == 0)                        // Am I the first row?
                      ? fromPrevKernel                  // If yes, then need to get score from previous chunk.
                      : _prevScores[tid - 1];           // Otherwise, I'm somewhere in the middle.

            const char textByte = textBytes[idxInRow - 1];
            const int scoreMatrixIdx = textByte * alphabetSize + patternByte;

            // Choose best direction, save score to shared and direction to global memory.
            auto scoreDirPair = choose_direction_SW(fromLeft, fromTop, fromDiag, gapPenalty,
                                                    _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            thisBestScoreIdx = scoreDirPair.first > thisBestScore
                               ? ((tid + startRow) * numCols + idxInRow) // Idx in the whole M matrix.
                               : thisBestScoreIdx;
            thisBestScore = max(thisBestScore, scoreDirPair.first);

            // If I'm I have finishes with a cell in the last row of this kernel,
            // then this column is done. Set done state and save score in global memory.
            if ((tid + startRow) == endRow)
                set_done(colState, idxInRow, scoreDirPair.first, kernelId);
        }
    }

    /* Second half of matrix filling */
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        const int idxInRow = (numCols - 1) + (i_pattern - tid);

        __syncthreads();

        if (tid >= i_pattern)
        {
            const char textByte = textBytes[idxInRow - 1];
            const int scoreMatrixIdx = textByte * alphabetSize + patternByte;

            // Choose best direction, save score to shared and direction to global memory.
            auto scoreDirPair = choose_direction_SW(_prevScores[tid], _prevScores[tid - 1], _prevPrevScores[tid - 1],
                                                    gapPenalty, _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            thisBestScoreIdx = scoreDirPair.first > thisBestScore
                               ? ((tid + startRow) * numCols + idxInRow) // Idx in the whole M matrix.
                               : thisBestScoreIdx;
            thisBestScore = max(thisBestScore, scoreDirPair.first);

            // If I'm I have finishes with a cell in the last row of this kernel,
            // then this column is done. Set done state and save score in global memory.
            if ((tid + startRow) == endRow)
                set_done(colState, idxInRow, scoreDirPair.first, kernelId);
        }
    }

    // Find the maximum score (and its index) accross all threads in the kernel.
    // Reuse the thisScores shared memory buffer to put all threads bestScores into consecutive
    // memory.
    _thisScores[tid] = thisBestScore;
    __syncthreads();
    max_reduce(_thisScores, blockDim.x);

    // Store best score and corresponding index in global memory. There could be multiple threads
    // with the highest score, but it doesn't matter which one is selected.
    if (_thisScores[0] == thisBestScore)
    {
        *p_bestScore = thisBestScore;
        *p_bestScoreIdx = thisBestScoreIdx;
        return;
    }

}


uint64_t initMemory(const SequenceAlignment::Request &request, SequenceAlignment::Response *response,
                    char *&d_textBytes, char *&d_patternBytes, int *&d_scoreMatrix,
                    columnState *&d_columnState, char *&os_M,
                    std::vector<char*> &d_M, std::vector<int*> &d_maxScore, std::vector<int*> &d_maxScoreIdx,
                    std::vector<char*> &h_M, std::vector<int*> &h_maxScore, std::vector<int*> &h_maxScoreIdx,
                    cudaStream_t &cuStream, const int numCuStreams)
{
    const uint64_t numCols = request.textNumBytes + 1;
    const uint64_t numRows = request.patternNumBytes + 1;

    try
    {
        os_M = new char[numRows * numCols];
        response->alignedTextBytes = new char[2 * request.textNumBytes];
        response->alignedPatternBytes = new char[2 * request.textNumBytes];
    }
    catch(const std::bad_alloc& e)
    {
        return 0;
    }

    // Select a number of threads per block such that we fit into global memory.
    auto calcGlobalMem = [&] (int numThreads)
    {
        return sizeof(int) * request.alphabetSize * request.alphabetSize +  // scoreMatrix
               request.textNumBytes + request.patternNumBytes +             // sequences
               sizeof(columnState) * numCols +                              // columState
               numCuStreams * numThreads * numCols +                        // device M per stream
               sizeof(int) * numCuStreams * 2;                              // bestScore + bestScoreIdx per stream
    };
    uint64_t numThreads = MAX_THREADS_PER_BLOCK;
    uint64_t freeGlobalMem = 0;
    cudaMemGetInfo((size_t*) &freeGlobalMem, 0);
    while (freeGlobalMem < calcGlobalMem(numThreads))
    {
        numThreads -= 32;
        if (numThreads < 32)
            return 0;
    }

    // Shared by all streams.
    if (cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols * sizeof(columnState)) != cudaSuccess)
    {
        return 0;
    }

    // One per stream.
    for (int i=0; i<numCuStreams; ++i)
    {
        if (cudaMalloc(&(d_M[i]), numThreads * numCols) != cudaSuccess ||
            cudaMalloc(&(d_maxScore[i]), sizeof(int)) != cudaSuccess ||
            cudaMalloc(&(d_maxScoreIdx[i]), sizeof(int)) != cudaSuccess ||
            cudaMallocHost(&(h_M[i]), numThreads * numCols) != cudaSuccess ||
            cudaMallocHost(&(h_maxScore[i]), sizeof(int)) != cudaSuccess ||
            cudaMallocHost(&(h_maxScoreIdx[i]), sizeof(int)) != cudaSuccess)
        {
            return 0;
        }
    }

    // Initialize the first row directions of M.
    // The init value will be different for global and local alignments.
    if (request.alignmentType == SequenceAlignment::programArgs::GLOBAL)
        std::fill_n(os_M, numCols, DIRECTION::LEFT);
    else if (request.alignmentType == SequenceAlignment::programArgs::LOCAL)
        std::fill_n(os_M, numCols, DIRECTION::STOP);

    // Zero columnState.
    cudaMemsetAsync(d_columnState, 0, sizeof(columnState) * numCols, cuStream);

    if (cudaMemcpyAsync(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice, cuStream) != cudaSuccess ||
        cudaMemcpyAsync(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice, cuStream) != cudaSuccess ||
        cudaMemcpyAsync(d_scoreMatrix, request.scoreMatrix, sizeof(int) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice, cuStream) != cudaSuccess)
    {
        return 0;
    }

    cudaStreamSynchronize(cuStream);

    return numThreads;
}

int SequenceAlignment::alignSequenceGPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{
    const uint64_t numCols = request.textNumBytes + 1;
    const uint64_t numRows = request.patternNumBytes + 1;

    // One CUDA thread block is limited to 1024 rows (thread limit).
    // Larger sequences are broken up into 1024 rows sized chunks.
    // Seperate CUDA streams for each SMs are started, and each stream gets assigned some rows, e.g.
    // stream0 gets rows 0-1023, and stream1 rows 1024-2047, ...
    // There is some degree of pipelining between streams - once stream{i} fills out all values in
    // a column for all its rows, then stream{i+1} can start filling out that column for its rows.
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    const int numCuStreams = std::min(MAX_CONCURRENT_KERNELS, deviceProp.multiProcessorCount);

    std::vector<cudaStream_t> cuStreams(numCuStreams);
    for (int i=0; i<numCuStreams; ++i)
        cudaStreamCreate(&(cuStreams[i]));
    // Used to identify current stream for CUDA operations. Start with stream0.
    cudaStream_t currStream = cuStreams[0];

    /** Memory allocation and transfer. */
    // OS managed memory (can be swapped to disk).
    char *os_M = nullptr;
    // CUDA managed host memory, pinned to physical mem address and not swappable.
    std::vector<char*> h_M(numCuStreams);
    std::vector<int*> h_maxScore(numCuStreams);
    std::vector<int*> h_maxScoreIdx(numCuStreams);
    // Device memory.
    char *d_textBytes = nullptr, *d_patternBytes = nullptr;
    int *d_scoreMatrix = nullptr;
    columnState *d_columnState = nullptr;
    std::vector<char*> d_M(numCuStreams);
    std::vector<int*> d_maxScore(numCuStreams);
    std::vector<int*> d_maxScoreIdx(numCuStreams);

    auto cleanUp = [&]()
    {
        for (int i=0; i<numCuStreams; ++i)
        {
            cudaStreamDestroy(cuStreams[i]);
            if (d_M[i]) cudaFree(d_M[i]);
            if (d_maxScore[i]) cudaFreeHost(d_maxScore[i]);
            if (d_maxScoreIdx[i]) cudaFreeHost(d_maxScoreIdx[i]);
            if (h_M[i]) cudaFreeHost(h_M[i]);
            if (h_maxScore[i]) cudaFreeHost(h_maxScore[i]);
            if (h_maxScoreIdx[i]) cudaFreeHost(h_maxScoreIdx[i]);
            d_M[i] = nullptr;
            d_maxScore[i] = nullptr;
            d_maxScoreIdx[i] = nullptr;
            h_M[i] = nullptr;
            h_maxScore[i] = nullptr;
            h_maxScoreIdx[i] = nullptr;
        }

        if (d_textBytes) cudaFree(d_textBytes);
        if (d_patternBytes) cudaFree(d_patternBytes);
        if (d_scoreMatrix) cudaFree(d_scoreMatrix);
        if (d_columnState) cudaFree(d_columnState);
        if (os_M) delete [] os_M;

        d_textBytes = nullptr;
        d_patternBytes = nullptr;
        d_scoreMatrix = nullptr;
        d_columnState = nullptr;
        os_M = nullptr;
    };

    const uint64_t NUM_THREADS_PER_BLOCK = initMemory(request, response, d_textBytes, d_patternBytes,
                                                      d_scoreMatrix, d_columnState, os_M,
                                                      d_M, d_maxScore, d_maxScoreIdx,
                                                      h_M, h_maxScore, h_maxScoreIdx,
                                                      currStream, numCuStreams);
    if (NUM_THREADS_PER_BLOCK == 0)
    {
        std::cout << MEM_ERROR;
        cleanUp();
        return -1;
    }
    /** End Allocate and transfer memory */

    // Three score buffers and the score matrix are kept in shared memory.
    const unsigned int sharedMemSize = 3 * std::min(NUM_THREADS_PER_BLOCK, numRows) * sizeof(int) +
                                       request.alphabetSize * request.alphabetSize * sizeof(int);

    #ifdef BENCHMARK
        auto begin = std::chrono::steady_clock::now();
    #endif

    // First row follows from inititalisation, start from second.
    int startRow = 1;
    auto curr_os_M = os_M + numCols;
    int maxScore = 0, maxScoreIdx = 0;
    for (int i_kernel=0; i_kernel < (numRows/NUM_THREADS_PER_BLOCK + 1); ++i_kernel)
    {
        const int numThreads = std::min(NUM_THREADS_PER_BLOCK, numRows - startRow);
        const int endRow = startRow + numThreads - 1;

        auto i_stream = i_kernel % numCuStreams;
        currStream = cuStreams[i_stream];

        if (request.alignmentType == programArgs::GLOBAL)
        {
            cuda_fillMatrixNW<<<1, numThreads, sharedMemSize, currStream>>>
                (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
                startRow, endRow, numCols, i_kernel, d_columnState, d_M[i_stream]);
        }
        else if (request.alignmentType == programArgs::LOCAL)
        {
            cuda_fillMatrixSW<<<1, numThreads, sharedMemSize, currStream>>>
                (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
                startRow, endRow, numCols, i_kernel, d_columnState, d_maxScore[i_stream],
                d_maxScoreIdx[i_stream], d_M[i_stream]);

            if (cudaMemcpyAsync(h_maxScore[i_stream], d_maxScore[i_stream], sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess ||
                cudaMemcpyAsync(h_maxScoreIdx[i_stream], d_maxScoreIdx[i_stream], sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
            {
                std::cout << "error: could not copy from device memory\n";
                cudaDeviceSynchronize();
                cleanUp();
                return -1;
            }
        }

        // Get the local M matrix of currStream into the large os_M matrix.
        // From device memory -> CUDA managed pinned memory -> OS managed swappable memory.
        if (cudaMemcpyAsync(h_M[i_stream], d_M[i_stream], numThreads*numCols, cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory\n";
            cudaDeviceSynchronize();
            cleanUp();
            return -1;
        }
        cudaStreamSynchronize(currStream);
        std::copy_n(h_M[i_stream], numThreads*numCols, curr_os_M);
        curr_os_M += numThreads*numCols;

        // Kepp track of maximum score and its index.
        if (request.alignmentType == programArgs::LOCAL)
        {
            maxScoreIdx = *(h_maxScore[i_stream]) > maxScore ? *(h_maxScoreIdx[i_stream]) : maxScoreIdx;
            maxScore = std::max(*(h_maxScore[i_stream]), maxScore);
        }

        startRow = endRow + 1;
    }

    #ifdef BENCHMARK
        // If benchmraking, return the time taken instead of error code.
        auto end = std::chrono::steady_clock::now();
        cleanUp();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    #endif

    if (request.alignmentType == programArgs::GLOBAL)
    {
        // At the end, pull out the score from the columnState.
        if (cudaMemcpyAsync(h_maxScore[0], &(d_columnState[numCols - 1].score), sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory\n";
            cleanUp();
            return -1;
        }

        cudaStreamSynchronize(currStream);
        response->score = *(h_maxScore[0]);
        traceBackNW(os_M, numRows, numCols, request, response);
    }
    else if (request.alignmentType == programArgs::LOCAL)
    {
        response->score = maxScore;
        traceBackSW(os_M, maxScoreIdx, numRows, numCols, request, response);
    }

    cleanUp();

    return 0;
}
