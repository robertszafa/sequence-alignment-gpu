#include <cuda.h>

#include "SequenceAlignment.hpp"

#define MAX_THREADS_PER_BLOCK 1024
#define MAX_CONCURRENT_KERNELS 32

using SequenceAlignment::DIRECTION;

struct __align__(8) columnState { int score; int kernelId; };


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

/// Global alignment using the Needleman-Wunsch algorithm.
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
    // In general, we look at 3 neighbours at each cell.
    // Additionally, we might have to look at the previous kernel's cell scores.
    int fromLeft = 0, fromDiag = 0, fromTop = 0, fromPrevKernel = 0;
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        // The diagonal of parallel executing cells grows by one, until reaching maximum.
        diagonalSize = min(diagonalSize + 1, numRows);
        const int idxInRow = i_text - tid;

        // Check if need to wait for previous kernel to start with row_0 of this kernel.
        if (tid == 0 && kernelId > 0)
            fromPrevKernel = get_prev_score(colState, i_text, kernelId);
        // Otherwise, can initialise.
        else if (tid == 0)
            fromPrevKernel = i_text * -gapPenalty;

        // Wait for scoreMatrix, thread 0 to return from get_prev_score, and/or previous iteration.
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

            // Choose best direction, save score to shared, and direction to global memory.
            auto scoreDirPair = choose_direction_NW(fromLeft, fromTop, fromDiag, gapPenalty,
                                                    _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            // Check if the 'idxInRow' column is done.
            // If yes, set done state and save score in global memory for next kernel.
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

            // Choose best direction, save score to shared, and direction to global memory.
            auto scoreDirPair = choose_direction_NW(_prevScores[tid], _prevScores[tid - 1], _prevPrevScores[tid - 1],
                                                    gapPenalty, _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            // Check if the 'idxInRow' column is done.
            // If yes, set done state and save score in global memory for next kernel.
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

/// Local alignment using the Smith-Waterman algorithm.
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
    // In general, we look at 3 neighbours at each cell.
    // Additionally, we might have to look at the previous kernel's cell scores.
    int fromLeft = 0, fromDiag = 0, fromTop = 0, fromPrevKernel = 0;
    // Each thread keeps track of the best score, and its index, in the whole kernel.
    int thisBestScore = 0, thisBestScoreIdx = 0;
    // The fromDiagonal score at iteration_i, becomes the fromTop score at itreation_i+1.
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        // The diagonal of parallel executing cells grows by one, until reaching maximum.
        diagonalSize = min(diagonalSize + 1, numRows);
        const int idxInRow = i_text - tid;

        // Check if need to wait for previous kernel to start with row_0 of this kernel.
        if (tid == 0 && kernelId > 0)
            fromPrevKernel = get_prev_score(colState, i_text, kernelId);
        // Otherwise, can initialise.
        else if (tid == 0)
            fromPrevKernel = 0;

        // Wait for scoreMatrix, thread 0 to return from get_prev_score, and/or previous iteration.
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

            // Choose best direction, save score to shared, and direction to global memory.
            auto scoreDirPair = choose_direction_SW(fromLeft, fromTop, fromDiag, gapPenalty,
                                                    _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            thisBestScoreIdx = scoreDirPair.first > thisBestScore
                               ? ((tid + startRow) * numCols + idxInRow) // Idx in the whole M matrix.
                               : thisBestScoreIdx;
            thisBestScore = max(thisBestScore, scoreDirPair.first);

            // Check if the 'idxInRow' column is done.
            // If yes, set done state and save score in global memory for next kernel.
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

            // Choose best direction, save score to shared, and direction to global memory.
            auto scoreDirPair = choose_direction_SW(_prevScores[tid], _prevScores[tid - 1], _prevPrevScores[tid - 1],
                                                    gapPenalty, _scoreMatrix[scoreMatrixIdx]);
            _thisScores[tid] = scoreDirPair.first;
            M[tid*numCols + idxInRow] = scoreDirPair.second;

            thisBestScoreIdx = scoreDirPair.first > thisBestScore
                               ? ((tid + startRow) * numCols + idxInRow) // Idx in the whole M matrix.
                               : thisBestScoreIdx;
            thisBestScore = max(thisBestScore, scoreDirPair.first);

            // Check if the 'idxInRow' column is done.
            // If yes, set done state and save score in global memory for next kernel.
            if ((tid + startRow) == endRow)
                set_done(colState, idxInRow, scoreDirPair.first, kernelId);
        }
    }

    // Find the maximum score (and its index) accross all threads in the kernel.
    // Reuse the thisScores buffer to put all threads' bestScores into addressable mem.
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


/// Allocate memory and initialise all data structures used in the NW, SW GPU algorithms.
/// Some buffers are replicated per CUDA stream - those are passed in as vectors of pointers.
///
/// Return the number of possible threads per block, given the size of text and pattern sequences.
/// Return 0 if not enough memory to create the minimum sized kernel of 32 x numCols.
uint64_t initMemory(const SequenceAlignment::Request &request, SequenceAlignment::Response *response,
                    char *&d_textBytes, char *&d_patternBytes, int *&d_scoreMatrix,
                    columnState *&d_columnState, std::vector<char*> &d_M,
                    std::vector<int*> &d_maxScore, std::vector<int*> &d_maxScoreIdx,
                    char *&h_M, int *&h_maxScores, int *&h_maxScoreIdxs,
                    cudaEvent_t *&finishedKernelEvents, cudaStream_t &cuStream, const int numComputeStreams)
{
    const uint64_t numCols = request.textNumBytes + 1;
    const uint64_t numRows = request.patternNumBytes + 1;

    // Calculate number of threads that can be supported with the input size. If <32, then abort.
    auto calcGlobalMem = [&] (int numThreads)
    {
        return sizeof(int) * request.alphabetSize * request.alphabetSize +  // scoreMatrix
               request.textNumBytes + request.patternNumBytes +             // sequences
               sizeof(columnState) * numCols +                              // columState
               numComputeStreams * numThreads * numCols +                   // one M per stream
               sizeof(int) * numComputeStreams * 2;                         // bestScore + bestScoreIdx per stream
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
    // How many GPU kernels will be started?
    const auto numKernels = numRows/numThreads + 1;

    // Start allocating memory

    /** Host Memory */
    try
    {
        response->alignedTextBytes = new char[2 * request.textNumBytes];
        response->alignedPatternBytes = new char[2 * request.textNumBytes];
    }
    catch(const std::bad_alloc& e)
    {
        return 0;
    }

    // Streams transfer their filled out M matrix portions to pinned memory. This allows the DMA
    // engine on the GPU to be used, and enables interleaving data transfer and computation.
    //
    // The downside is that the OS cannot swap this memory - the size limit of h_M is
    // the physical memory size, rather the virtual memory size.
    if (cudaMallocHost(&h_M, numRows * numCols) != cudaSuccess)
        return 0;
    // One event per kernel to start DMA transfer.
    finishedKernelEvents = new cudaEvent_t[numKernels];
    // For local alignment:
    if (cudaMallocHost(&(h_maxScores), numKernels * sizeof(int)) != cudaSuccess ||
        cudaMallocHost(&(h_maxScoreIdxs), numKernels * sizeof(int)) != cudaSuccess)
    {
        return 0;
    }

    // Initialize the first row directions of M (different for global and local alignments).
    if (request.alignmentType == SequenceAlignment::programArgs::GLOBAL)
        std::fill_n(h_M, numCols, DIRECTION::LEFT);
    else if (request.alignmentType == SequenceAlignment::programArgs::LOCAL)
        std::fill_n(h_M, numCols, DIRECTION::STOP);
    /** End Host Memory */


    /** Device Memory */
    // Shared by all streams.
    if (cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols * sizeof(columnState)) != cudaSuccess)
    {
        return 0;
    }

    // One per stream.
    for (int i=0; i<numComputeStreams; ++i)
    {
        if (cudaMalloc(&(d_M[i]), numThreads * numCols) != cudaSuccess ||
            cudaMalloc(&(d_maxScore[i]), sizeof(int)) != cudaSuccess ||
            cudaMalloc(&(d_maxScoreIdx[i]), sizeof(int)) != cudaSuccess)
        {
            return 0;
        }
    }

    // Zero columnState.
    cudaMemsetAsync(d_columnState, 0, sizeof(columnState) * numCols, cuStream);

    if (cudaMemcpyAsync(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice, cuStream) != cudaSuccess ||
        cudaMemcpyAsync(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice, cuStream) != cudaSuccess ||
        cudaMemcpyAsync(d_scoreMatrix, request.scoreMatrix, sizeof(int) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice, cuStream) != cudaSuccess)
    {
        return 0;
    }
    // Wait for all mem ops to finish.
    cudaStreamSynchronize(cuStream);
    /** End Device Memory */

    return numThreads;
}

int SequenceAlignment::alignSequenceGPU(const SequenceAlignment::Request &request,
                                        SequenceAlignment::Response *response)
{
    const uint64_t numCols = request.textNumBytes + 1;
    const uint64_t numRows = request.patternNumBytes + 1;

    // One CUDA thread block is limited to 1024 rows (thread limit).
    // Larger sequences are broken up into 1024-sized chunks.
    // Seperate CUDA streams for each SMs are started, and each stream gets assigned some rows, e.g.
    // stream0 gets rows 0-1023, and stream1 rows 1024-2047, ...
    // There is some degree of pipelining between streams - once stream{i} fills out all values in
    // a column for all its rows, then stream{i+1} can start filling out that column for its rows.
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    // Ensure range 1 >= numComputeStreams <= MAX_CONCURRENT_KERNELS, with the conditions:
    //      - one stream per 1024 rows,
    //      - no more streams than the number of available SMs.
    const int numComputeStreams = std::min(std::max((int) (numRows + 1) / MAX_THREADS_PER_BLOCK, 1),
                                           std::min(MAX_CONCURRENT_KERNELS, deviceProp.multiProcessorCount));
    // Will hold the number of kernel invocations required for the input sequences.
    int numKernels = 0;
    std::vector<cudaStream_t> computeStreams(numComputeStreams);
    cudaStream_t memoryStream;

    cudaStreamCreate(&(memoryStream));
    for (int i=0; i<numComputeStreams; ++i)
        cudaStreamCreate(&(computeStreams[i]));

    // Used to identify current stream for CUDA operations. Start with stream0.
    cudaStream_t currComputeStream = computeStreams[0];
    // Each compute stream will need to notify the memoryStream that it finished with a kernel
    // and its M matrix is ready to be transfered to host memory. This is done through CUDA events.
    cudaEvent_t *finishedKernelEvents = nullptr;

    /** Memory allocation and transfer. */
    // CUDA managed host memory, pinned to physical mem address and not swappable by the OS:
    char *h_M = nullptr;
    // Local alignment: each kernel will result in a bestScore and bestIdx instance.
    int *h_maxScores = nullptr, *h_maxScoreIdxs = nullptr;
    // Device memory:
    // Shared by all streams.
    char *d_textBytes = nullptr, *d_patternBytes = nullptr;
    int *d_scoreMatrix = nullptr;
    columnState *d_columnState = nullptr;
    // Private to a stream.
    std::vector<char*> d_M(numComputeStreams);
    std::vector<int*> d_maxScore(numComputeStreams);
    std::vector<int*> d_maxScoreIdx(numComputeStreams);

    // Define how the allocated resources are destroyed.
    auto cleanUp = [&]()
    {
        cudaStreamDestroy(memoryStream);
        for (int i=0; i<numComputeStreams; ++i)
        {
            cudaStreamDestroy(computeStreams[i]);
            if (d_M[i]) cudaFree(d_M[i]);
            if (d_maxScore[i]) cudaFreeHost(d_maxScore[i]);
            if (d_maxScoreIdx[i]) cudaFreeHost(d_maxScoreIdx[i]);
            d_M[i] = nullptr;
            d_maxScore[i] = nullptr;
            d_maxScoreIdx[i] = nullptr;
        }

        for (int i=0; i<numKernels; ++i)
            cudaEventDestroy(finishedKernelEvents[i]);

        if (finishedKernelEvents) delete [] finishedKernelEvents;
        if (h_M) cudaFreeHost(h_M);
        if (h_maxScores) cudaFreeHost(h_maxScores);
        if (h_maxScoreIdxs) cudaFreeHost(h_maxScoreIdxs);
        if (d_textBytes) cudaFree(d_textBytes);
        if (d_patternBytes) cudaFree(d_patternBytes);
        if (d_scoreMatrix) cudaFree(d_scoreMatrix);
        if (d_columnState) cudaFree(d_columnState);

        finishedKernelEvents = nullptr;
        h_M = nullptr;
        h_maxScores = nullptr;
        h_maxScoreIdxs = nullptr;
        d_textBytes = nullptr;
        d_patternBytes = nullptr;
        d_scoreMatrix = nullptr;
        d_columnState = nullptr;
    };

    // Deal with memory, and find out the number of possible threads per block, given the input size.
    // Make the first compute stream deal with memory, to ensure all is ready when it starts
    // (operations are ordered within streams).
    const uint64_t NUM_THREADS_PER_BLOCK = initMemory(request, response, d_textBytes, d_patternBytes,
                                                      d_scoreMatrix, d_columnState, d_M, d_maxScore,
                                                      d_maxScoreIdx, h_M, h_maxScores, h_maxScoreIdxs,
                                                      finishedKernelEvents, currComputeStream, numComputeStreams);
    if (NUM_THREADS_PER_BLOCK == 0)
    {
        std::cout << MEM_ERROR;
        cleanUp();
        return -1;
    }
    numKernels = (numRows/NUM_THREADS_PER_BLOCK) + 1;
    /** End Allocate and transfer memory */

    // Three score buffers and the score matrix are kept in shared memory.
    const unsigned int sharedMemSize = 3 * std::min(NUM_THREADS_PER_BLOCK, numRows) * sizeof(int) +
                                       request.alphabetSize * request.alphabetSize * sizeof(int);

    #ifdef BENCHMARK
        auto begin = std::chrono::steady_clock::now();
    #endif

    // First row and first col follow from inititalisation, start from second.
    int startRow = 1;
    auto curr_h_M = h_M + numCols;
    for (int i_kernel=0; i_kernel < numKernels; ++i_kernel)
    {
        const int numThreads = std::min(NUM_THREADS_PER_BLOCK, numRows - startRow);
        const int endRow = startRow + numThreads - 1;

        // Round-robin scheduling for compute streams.
        auto i_stream = i_kernel % (numComputeStreams);
        currComputeStream = computeStreams[i_stream];
        cudaEventCreateWithFlags(&finishedKernelEvents[i_kernel], cudaEventDisableTiming);

        if (request.alignmentType == programArgs::GLOBAL)
        {
            cuda_fillMatrixNW<<<1, numThreads, sharedMemSize, currComputeStream>>>
                (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
                startRow, endRow, numCols, i_kernel, d_columnState, d_M[i_stream]);
        }
        else if (request.alignmentType == programArgs::LOCAL)
        {
            cuda_fillMatrixSW<<<1, numThreads, sharedMemSize, currComputeStream>>>
                (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
                startRow, endRow, numCols, i_kernel, d_columnState, d_maxScore[i_stream],
                d_maxScoreIdx[i_stream], d_M[i_stream]);
        }

        // Record kernel completion.
        cudaEventRecord(finishedKernelEvents[i_kernel], currComputeStream);
        // memoryStream waits for kernelCompletion before transferring d_M into h_M
        cudaStreamWaitEvent(memoryStream, finishedKernelEvents[i_kernel], 0);

        if (cudaMemcpyAsync(curr_h_M, d_M[i_stream], numThreads*numCols, cudaMemcpyDeviceToHost, memoryStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory\n";
            cudaDeviceSynchronize();
            cleanUp();
            return -1;
        }
        curr_h_M += numThreads*numCols;

        // Local alignment: keep track of maximum score and its index.
        if (request.alignmentType == programArgs::LOCAL)
        {
            // Get maxScore and its index from each local alignment kernel.
            if (cudaMemcpyAsync(h_maxScores + i_kernel, d_maxScore[i_stream], sizeof(int), cudaMemcpyDeviceToHost, memoryStream) != cudaSuccess ||
                cudaMemcpyAsync(h_maxScoreIdxs + i_kernel, d_maxScoreIdx[i_stream], sizeof(int), cudaMemcpyDeviceToHost, memoryStream) != cudaSuccess)
            {
                std::cout << "error: could not copy from device memory\n";
                cudaDeviceSynchronize();
                cleanUp();
                return -1;
            }
        }

        startRow = endRow + 1;
    }

    #ifdef BENCHMARK
        // If benchmarking, return the time taken instead of error code.
        // Measure with data transfer back to the host.
        cudaStreamSynchronize(memoryStream);
        auto end = std::chrono::steady_clock::now();
        cleanUp();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    #endif

    if (request.alignmentType == programArgs::GLOBAL)
    {
        // At the end, pull out the score from the columnState.
        if (cudaMemcpyAsync(h_maxScores, &(d_columnState[numCols - 1].score), sizeof(int), cudaMemcpyDeviceToHost, memoryStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory\n";
            cleanUp();
            return -1;
        }

        cudaStreamSynchronize(memoryStream);
        response->score = h_maxScores[0];
        traceBackNW(h_M, numRows, numCols, request, response);
    }
    else if (request.alignmentType == programArgs::LOCAL)
    {
        // Find the maximum kernel score out of every kernel invocation.
        cudaStreamSynchronize(memoryStream);
        auto iMax = std::distance(h_maxScores, std::max_element(h_maxScores, (h_maxScores + numKernels)));
        response->score = h_maxScores[iMax];
        traceBackSW(h_M, h_maxScoreIdxs[iMax], numRows, numCols, request, response);
    }

    cleanUp();
    return 0;
}
