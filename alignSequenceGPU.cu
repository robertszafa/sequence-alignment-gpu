#include <cuda.h>

#include "SequenceAlignment.hpp"

constexpr uint64_t MAX_THREADS_PER_BLOCK = 1024;

using SequenceAlignment::DIR;

struct columnState { int score; int kernelId; };

__device__ __forceinline__ void set_done(columnState* volatile colState, const unsigned int col,
                                         const int score, const int kernelId)
{
    colState[col] = {score, kernelId+1};
}

__device__ __forceinline__ void busy_wait(columnState* volatile colState, const unsigned int col,
                                          const int kernelId)
{
    if (threadIdx.x == 0)
    {
        volatile int currKernelId = colState[col].kernelId;
        while (currKernelId != kernelId)
            currKernelId = colState[col].kernelId;
    }
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
    *directionStore = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;

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

    const char patternByte = patternBytes[max(0, (tid+startRow) - 1)];

    // First half of matrix filling
    int diag_size = 0;
    for (int i_text = 0; i_text < numCols; ++i_text)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        diag_size = min(diag_size+1, numRows);
        const int threadInRowIdx = i_text - tid;

        if (startRow > 0)
            busy_wait(colState, i_text, kernelId);

        __syncthreads();

        if (tid == 0 && startRow == 0) // First row.
        {
            _thisScores[tid] = -(i_text * gapPenalty);
            M[i_text] = DIR::LEFT;
        }
        else if (tid == (diag_size-1) && i_text < numRows) // First column.
        {
            _thisScores[tid] = -((tid + startRow) * gapPenalty);
            M[tid * numCols] = DIR::TOP;
        }
        else if (tid == 0) // Not first row of M, but first row in this kernel.
        {
            const char textByte = textBytes[threadInRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            _thisScores[tid] = choose_direction(_prevScores[tid],
                                                colState[i_text].score,
                                                colState[i_text - 1].score,
                                                gapPenalty, _scoreMatrix[scoreMatrixIdx],
                                                (M + tid*numCols + threadInRowIdx));
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

        if ((tid + startRow) == endRow && tid < diag_size)
            set_done(colState, threadInRowIdx, _thisScores[tid], kernelId);
    }

    // Second half of matrix filling.
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        ping_pong_buffers(_thisScores, _prevScores, _prevPrevScores);

        const int threadInRowIdx = numCols-1 - tid + i_pattern;

        __syncthreads();

        if (tid >= i_pattern)
        {
            const char textByte = textBytes[threadInRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            _thisScores[tid] = choose_direction(_prevScores[tid],
                                                _prevScores[tid - 1],
                                                _prevPrevScores[tid - 1],
                                                gapPenalty, _scoreMatrix[scoreMatrixIdx],
                                                (M + tid*numCols + threadInRowIdx));

            if ((tid + startRow) == endRow)
                set_done(colState, threadInRowIdx, _thisScores[tid], kernelId);
        }
    }
}


int SequenceAlignment::alignSequenceGlobalGPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{
    const uint64_t numCols = request.textNumBytes + 1;
    const uint64_t numRows = request.patternNumBytes + 1;

    /** Memory allocation and transfer. */
    cudaStream_t stream0, stream1;
    cudaStreamCreate(&stream0);
    cudaStreamCreate(&stream1);

    // OS managed memory (can be swapped to disk).
    char *os_M;
    // CUDA managed host memory, pinned to physical mem address and not swappable.
    char *h_M0, *h_M1;
    int *h_score;
    // Device memory.
    char *d_textBytes, *d_patternBytes, *d_M0, *d_M1;
    int *d_scoreMatrix;
    columnState *d_columnState;

    auto cleanUp = [&]()
    {
        cudaStreamDestroy(stream0);
        cudaStreamDestroy(stream1);

        if (d_textBytes) cudaFree(d_textBytes);
        if (d_patternBytes) cudaFree(d_patternBytes);
        if (d_M0) cudaFree(d_M0);
        if (d_M1) cudaFree(d_M1);
        if (d_scoreMatrix) cudaFree(d_scoreMatrix);
        if (d_columnState) cudaFree(d_columnState);

        if (h_M0) cudaFreeHost(h_M0);
        if (h_M1) cudaFreeHost(h_M1);
        if (h_score) cudaFreeHost(h_score);

        if (os_M) delete [] os_M;
    };

    try
    {
        os_M = new char[numRows * numCols];
        response->alignedTextBytes = new char[2 * request.textNumBytes];
        response->alignedPatternBytes = new char[2 * request.textNumBytes];
    }
    catch(const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        cleanUp();
        return -1;
    }

    if (cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M0, MAX_THREADS_PER_BLOCK * numCols) != cudaSuccess ||
        cudaMalloc(&d_M1, MAX_THREADS_PER_BLOCK * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_columnState, numCols * sizeof(columnState)) != cudaSuccess ||
        cudaMallocHost(&h_M0, MAX_THREADS_PER_BLOCK * numCols) != cudaSuccess ||
        cudaMallocHost(&h_M1, MAX_THREADS_PER_BLOCK * numCols) != cudaSuccess ||
        cudaMallocHost(&h_score, sizeof(int)) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cleanUp();
        return -1;
    }
    if (cudaMemcpyAsync(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice, stream0) != cudaSuccess ||
        cudaMemcpyAsync(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice, stream0) != cudaSuccess ||
        cudaMemcpyAsync(d_scoreMatrix, request.scoreMatrix, sizeof(int) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice, stream0) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cleanUp();
        return -1;
    }

    cudaMemsetAsync(d_columnState, 0, sizeof(columnState) * numCols, stream0);
    /** End Memory allocation and transfer. */


    #ifdef BENCHMARK
        auto begin = std::chrono::steady_clock::now();
    #endif

    const unsigned int sharedMemSize = 3 * std::min(MAX_THREADS_PER_BLOCK, numRows) * sizeof(int) +
                                       request.alphabetSize * request.alphabetSize * sizeof(int);

    int startRow = 0;
    cudaStream_t currStream;
    auto curr_os_M = os_M;
    for (int i_kernel=0; i_kernel < (numRows/MAX_THREADS_PER_BLOCK + 1); ++i_kernel)
    {
        const int numThreads = std::min(MAX_THREADS_PER_BLOCK, numRows - startRow);
        const int endRow = startRow + numThreads - 1;

        currStream = (i_kernel % 2 == 0) ? stream0 : stream1;
        auto curr_d_M = (i_kernel % 2 == 0) ? d_M0 : d_M1;
        auto curr_h_M = (i_kernel % 2 == 0) ? h_M0 : h_M1;

        cuda_fillMatrixNW<<<1, numThreads, sharedMemSize, currStream>>>
            (d_textBytes, d_patternBytes, d_scoreMatrix, request.alphabetSize, request.gapPenalty,
             startRow, endRow, numCols, i_kernel, d_columnState, curr_d_M);

        // Get the filled out part of M matrix in this iteration.
        // From device memory -> CUDA managed pinned memory -> OS managed swappable memory.
        if (cudaMemcpyAsync(curr_h_M, curr_d_M, numThreads*numCols, cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
        {
            std::cout << "error: could not copy from device memory 1\n";
            cudaDeviceSynchronize();
            cleanUp();
            return -1;
        }

        cudaStreamSynchronize(currStream);
        std::copy_n(curr_h_M, numThreads*numCols, curr_os_M);
        curr_os_M += numThreads*numCols;

        startRow = endRow + 1;
    }

    if (cudaMemcpyAsync(h_score, &(d_columnState[numCols - 1].score), sizeof(int), cudaMemcpyDeviceToHost, currStream) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        cleanUp();
        return -1;
    }

    cudaStreamSynchronize(currStream);
    response->score = *h_score;

    #ifdef BENCHMARK
        auto end = std::chrono::steady_clock::now();
        cleanUp();
        return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    #endif

    traceBack(os_M, numRows, numCols, request, response);

    cleanUp();

    return 0;
}
