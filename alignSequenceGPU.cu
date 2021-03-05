#include <cuda.h>

#include "SequenceAlignment.hpp"

constexpr uint64_t MAX_THREADS_PER_BLOCK = 1024;

using SequenceAlignment::DIR;

struct columnState { int score; int kernelId; };


/// Given a column and a kernelId, set the kernelId of that column to the supplied kernelId.
__device__ __forceinline__ void set_done(columnState* volatile colState, const unsigned int col,
                                         const int score, const int kernelId)
{
    colState[col] = {score, kernelId+1};
}

/// Given a column and a kernelId, spin while the kernelId of that column
/// is not equal to the supplied kernelId.
__device__ __forceinline__ void busy_wait(columnState* volatile colState, const unsigned int col,
                                          const int kernelId)
{
    volatile int currKernelId = colState[col].kernelId;
    while (currKernelId != kernelId)
        currKernelId = colState[col].kernelId;
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
__device__ __forceinline__ std::pair<int, DIR> choose_direction_NW(const int leftScore,
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

    const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
    const auto bestDir = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;

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

    // Each thread gets one row (one pattern letter).
    const char patternByte = patternBytes[max(0, (tid + startRow) - 1)];
    M[tid * numCols] = DIR::TOP;

    /* First half of matrix filling */
    int diagonalSize = 0;
    int fromLeft, fromDiag;
    // The fromDiagonal score at iteration_i, becomes the fromTop score at itreation_i+1.
    int fromTop = (tid + startRow-1) * gapPenalty;
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        diagonalSize = min(diagonalSize + 1, numRows);
        const int idxInRow = i_text - tid;

        // Need to wait for previous chunk to finish with this column.
        if (tid == 0)
            busy_wait(colState, i_text, kernelId);

        // Join previous iteration's threads and thread 0.
        __syncthreads();

        if (tid < diagonalSize)
        {
            fromLeft = (idxInRow == 1)                  // Am I in the first column to the left?
                       ? (tid + startRow) * gapPenalty  // If yes, then column to left is the empty char.
                       : _prevScores[tid];              // Otherwise, I'm somewhere in the middle.
            fromDiag = fromTop;
            fromTop = (tid == 0)                        // Am I the first row?
                      ? colState[i_text].score          // If yes, then need to get score from previous chunk.
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
/// return the best score and direction. If the best score is 0, then return DIR::STOP.
__device__ __forceinline__ std::pair<int, DIR> choose_direction_SW(const int leftScore,
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

    const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
    const auto bestDirNonZero = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
    const auto bestDir = bestScore > 0 ? bestDirNonZero : DIR::STOP;

    return {bestScore, bestDir};
}

__device__ __forceinline__ std::pair<int, int> max_reduce(const int threadBestScore,
                                                          const int threadBestScoreIdx)
{
    int tid = threadIdx.x;

    __shared__ volatile int s_bestScore;
    __shared__ volatile int s_bestScoreIdx;

    if (tid == 0)
    {
        s_bestScore = 0;
        s_bestScoreIdx = 0;
    }

    // Do all warps.
    for (int i_warp = 0; i_warp < MAX_THREADS_PER_BLOCK; i_warp += 32)
    {
        int i = i_warp + tid;
        if (i >= blockDim.x) break;

        // Reduce max inside warp. This while loop is guaranteed to terminate.
        while (s_bestScore < threadBestScore)
        {
            s_bestScore = threadBestScore;
            s_bestScoreIdx = threadBestScoreIdx;
        }

        __syncthreads();
    }

    return {s_bestScore, s_bestScoreIdx};
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

    // Each thread gets one row (one pattern letter).
    const char patternByte = patternBytes[max(0, (tid + startRow) - 1)];
    M[tid * numCols] = DIR::STOP;

    /* First half of matrix filling */
    int diagonalSize = 0;
    int fromLeft = 0, fromDiag = 0, fromTop = 0;
    int thisBestScore = 0;
    int thisBestScoreIdx = startRow * numCols;
    // The fromDiagonal score at iteration_i, becomes the fromTop score at itreation_i+1.
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        diagonalSize = min(diagonalSize + 1, numRows);
        const int idxInRow = i_text - tid;

        // Need to wait for previous chunk to finish with this column.
        if (tid == 0)
            busy_wait(colState, i_text, kernelId);

        // Join previous iteration's threads and thread 0.
        __syncthreads();

        if (tid < diagonalSize)
        {
            fromLeft = (idxInRow == 1)                  // Am I in the first column to the left?
                       ? 0                              // If yes, then column to left has score 0.
                       : _prevScores[tid];              // Otherwise, I'm somewhere in the middle.
            fromDiag = fromTop;
            fromTop = (tid == 0)                        // Am I the first row?
                      ? colState[i_text].score          // If yes, then need to get score from previous chunk.
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

    // Find maximum score (and its index) accross all threads in block.
    auto bestScoreAndIdxPair = max_reduce(thisBestScore, thisBestScoreIdx);
    if (tid == 0)
    {
        *p_bestScore = bestScoreAndIdxPair.first;
        *p_bestScoreIdx = bestScoreAndIdxPair.second;
    }
}


uint64_t initMemory(const SequenceAlignment::Request &request, SequenceAlignment::Response *response,
                    char *&d_M0, char *&d_M1, char *&d_textBytes, char *&d_patternBytes,
                    int *&d_scoreMatrix, columnState *&d_columnState, int *&d_maxScore0, int *&d_maxScore1,
                    int *&d_maxScoreIdx0, int *&d_maxScoreIdx1, char *&h_M0, char *&h_M1, int *&h_maxScore0,
                    int *&h_maxScore1, int *&h_maxScoreIdx0, int *&h_maxScoreIdx1, char *&os_M, cudaStream_t &cuStream)
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
               2 * numThreads * numCols +                                   // M0, M1
               request.textNumBytes + request.patternNumBytes +             // sequences
               sizeof(columnState) * numCols +                              // columState
               sizeof(int) * 4;                                             // 2*bestScore + 2*bestScoreIdx
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

    if (cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M0, numThreads * numCols) != cudaSuccess ||
        cudaMalloc(&d_M1, numThreads * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols * sizeof(columnState)) != cudaSuccess ||
        cudaMalloc(&d_maxScore0, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_maxScore1, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_maxScoreIdx0, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_maxScoreIdx1, sizeof(int)) != cudaSuccess ||
        cudaMallocHost(&h_M0, numThreads * numCols) != cudaSuccess ||
        cudaMallocHost(&h_M1, numThreads * numCols) != cudaSuccess ||
        cudaMallocHost(&h_maxScore0, sizeof(int)) != cudaSuccess ||
        cudaMallocHost(&h_maxScore1, sizeof(int)) != cudaSuccess ||
        cudaMallocHost(&h_maxScoreIdx0, sizeof(int)) != cudaSuccess ||
        cudaMallocHost(&h_maxScoreIdx1, sizeof(int)) != cudaSuccess)
    {
        return 0;
    }

    // Initialize the very first row scores and directions of M. Also, initialize columnState.
    // The init value will be different for global and local alignments.
    if (request.alignmentType == SequenceAlignment::programArgs::GLOBAL)
    {
        std::fill_n(os_M, numCols, DIR::LEFT);
        std::vector<columnState> initState(numCols);
        for (int i=0; i<numCols; ++i)
        {
            initState[i].score = i * -request.gapPenalty;
            initState[i].kernelId = 0;
        }
        if (cudaMemcpyAsync(d_columnState, &(initState[0]), sizeof(columnState) * numCols, cudaMemcpyHostToDevice, cuStream) != cudaSuccess)
            return 0;
    }
    else if (request.alignmentType == SequenceAlignment::programArgs::LOCAL)
    {
        std::fill_n(os_M, numCols, DIR::STOP);
        cudaMemsetAsync(d_columnState, 0, sizeof(columnState) * numCols, cuStream);
    }

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
    // Two seperate CUDA streams execute chunk{i} and chunk{i+1}, using the fact that once
    // stream{j} fills out all values in a column, then stream{j+1} can start filling out that column
    // for later rows - achieving some amount of pipelining.
    cudaStream_t stream0, stream1;
    cudaStreamCreate(&stream0);
    cudaStreamCreate(&stream1);
    // Used to identify current stream for CUDA operations. Start with stream0.
    cudaStream_t currStream = stream0;

    /** Memory allocation and transfer. */
    // OS managed memory (can be swapped to disk).
    char *os_M = nullptr;
    // CUDA managed host memory, pinned to physical mem address and not swappable.
    char *h_M0 = nullptr, *h_M1 = nullptr;
    int *h_maxScore0 = nullptr, *h_maxScore1 = nullptr;
    int *h_maxScoreIdx0 = nullptr, *h_maxScoreIdx1 = nullptr;
    // Device memory.
    char *d_textBytes = nullptr, *d_patternBytes = nullptr, *d_M0 = nullptr, *d_M1 = nullptr;
    int *d_scoreMatrix = nullptr;
    columnState *d_columnState = nullptr;
    int *d_maxScore0 = nullptr, *d_maxScore1 = nullptr;
    int *d_maxScoreIdx0 = nullptr, *d_maxScoreIdx1 = nullptr;

    auto cleanUp = [&]()
    {
        cudaStreamDestroy(stream0);
        cudaStreamDestroy(stream1);

        if (d_textBytes) cudaFree(d_textBytes);
        if (d_patternBytes) cudaFree(d_patternBytes);
        if (d_M0) cudaFree(d_M0);
        if (d_M1) cudaFree(d_M1);
        if (d_scoreMatrix) cudaFree(d_scoreMatrix);
        if (d_maxScore0) cudaFreeHost(d_maxScore0);
        if (d_maxScore1) cudaFreeHost(d_maxScore1);
        if (d_maxScoreIdx0) cudaFreeHost(d_maxScoreIdx0);
        if (d_maxScoreIdx1) cudaFreeHost(d_maxScoreIdx1);
        if (d_columnState) cudaFree(d_columnState);
        if (h_M0) cudaFreeHost(h_M0);
        if (h_M1) cudaFreeHost(h_M1);
        if (h_maxScore0) cudaFreeHost(h_maxScore0);
        if (h_maxScore1) cudaFreeHost(h_maxScore1);
        if (h_maxScoreIdx0) cudaFreeHost(h_maxScoreIdx0);
        if (h_maxScoreIdx1) cudaFreeHost(h_maxScoreIdx1);
        if (os_M) delete [] os_M;

        d_textBytes = nullptr;
        d_patternBytes = nullptr;
        d_M0 = nullptr;
        d_M1 = nullptr;
        d_scoreMatrix = nullptr;
        d_columnState = nullptr;
        d_maxScore0 = nullptr;
        d_maxScore1 = nullptr;
        d_maxScoreIdx0 = nullptr;
        d_maxScoreIdx1 = nullptr;
        h_M0 = nullptr;
        h_M1 = nullptr;
        h_maxScore0 = nullptr;
        h_maxScore1 = nullptr;
        h_maxScoreIdx0 = nullptr;
        h_maxScoreIdx1 = nullptr;
        os_M = nullptr;
    };

    const uint64_t NUM_THREADS_PER_BLOCK = initMemory(request, response, d_M0, d_M1, d_textBytes,
                                                      d_patternBytes, d_scoreMatrix, d_columnState,
                                                      d_maxScore0, d_maxScore1, d_maxScoreIdx0, d_maxScoreIdx1,
                                                      h_M0, h_M1, h_maxScore0, h_maxScore1, h_maxScoreIdx0,
                                                      h_maxScoreIdx1, os_M, currStream);
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

    // First row is already inititalized, start from second.
    int startRow = 1;
    auto curr_os_M = os_M + numCols;
    int maxScore = 0, maxScoreIdx = 0;
    for (int i_kernel=0; i_kernel < (numRows/NUM_THREADS_PER_BLOCK + 1); ++i_kernel)
    {
        const int numThreads = std::min(NUM_THREADS_PER_BLOCK, numRows - startRow);
        const int endRow = startRow + numThreads - 1;

        currStream = (i_kernel % 2 == 0) ? stream0 : stream1;
        // Each stream has its local M matrix in pinned memory.
        // This allows concurrent copy operations.
        auto curr_d_M = (i_kernel % 2 == 0) ? d_M0 : d_M1;
        auto curr_h_M = (i_kernel % 2 == 0) ? h_M0 : h_M1;
        // In local alignment, we need to keep track of the maximum score (and idx)
        // inside the whole M matrix. Thus, the below are used for local alignment only.
        auto curr_d_score = (i_kernel % 2 == 0) ? d_maxScore0 : d_maxScore1;
        auto curr_d_idx = (i_kernel % 2 == 0) ? d_maxScoreIdx0 : d_maxScoreIdx1;
        auto curr_h_score = (i_kernel % 2 == 0) ? h_maxScore0 : h_maxScore1;
        auto curr_h_idx = (i_kernel % 2 == 0) ? h_maxScoreIdx0 : h_maxScoreIdx1;

        if (request.alignmentType == programArgs::GLOBAL)
        {
            cuda_fillMatrixNW<<<1, numThreads, sharedMemSize, currStream>>>
                (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
                startRow, endRow, numCols, i_kernel, d_columnState, curr_d_M);
        }
        else if (request.alignmentType == programArgs::LOCAL)
        {
            cuda_fillMatrixSW<<<1, numThreads, sharedMemSize, currStream>>>
                (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
                startRow, endRow, numCols, i_kernel, d_columnState, curr_d_score, curr_d_idx, curr_d_M);

            if (cudaMemcpyAsync(curr_h_score, curr_d_score, sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess ||
                cudaMemcpyAsync(curr_h_idx, curr_d_idx, sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
            {
                std::cout << "error: could not copy from device memory\n";
                cudaDeviceSynchronize();
                cleanUp();
                return -1;
            }
        }

        // Get the local M matrix of currStream into the large os_M matrix.
        // From device memory -> CUDA managed pinned memory -> OS managed swappable memory.
        if (cudaMemcpyAsync(curr_h_M, curr_d_M, numThreads*numCols, cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory\n";
            cudaDeviceSynchronize();
            cleanUp();
            return -1;
        }
        cudaStreamSynchronize(currStream);
        std::copy_n(curr_h_M, numThreads*numCols, curr_os_M);
        curr_os_M += numThreads*numCols;

        // Kepp track of maximum score and its index.
        if (request.alignmentType == programArgs::LOCAL)
        {
            maxScoreIdx = *curr_h_score > maxScore ? *curr_h_idx : maxScoreIdx;
            maxScore = std::max(*curr_h_score, maxScore);
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
        if (cudaMemcpyAsync(h_maxScore0, &(d_columnState[numCols - 1].score), sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory\n";
            cleanUp();
            return -1;
        }

        cudaStreamSynchronize(currStream);
        response->score = *h_maxScore0;
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
