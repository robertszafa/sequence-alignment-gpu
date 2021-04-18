/// Old GPU version of the Needleman-Wunsch DP matrix filling algorithm.
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
    using SequenceAlignment::DIRECTION;

    extern __shared__ int _shared[];
    int *_thisScores = _shared;
    int *_prevScores = _shared + numCols;

    const int tid = threadIdx.x;

    // Each thread copies one text letter.
    const char textByte = (tid > 0) ? textBytes[tid - 1] : alphabetSize;
    // Init first row.
    _thisScores[tid] = tid * -gapPenalty;
    M[tid] = DIRECTION::LEFT;

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
            thisRowM[0] = DIRECTION::TOP;
            continue;
        }

        const char patternByte = patternBytes[i_pattern - 1];
        const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

        // We are accessing the previous row - wait for all columns to finish.
        const int fromTopScore = _prevScores[tid] - gapPenalty;
        const int fromDiagScore = _prevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

        const bool isDiagGreaterThanTop = (fromDiagScore > fromTopScore);
        const int maxFromPrev = isDiagGreaterThanTop ? fromDiagScore : fromTopScore;
        const auto tmpDir = isDiagGreaterThanTop ? DIRECTION::DIAG : DIRECTION::TOP;

        for (int i_text = 1; i_text < numCols; ++i_text)
        {
            // We are accessing the previous column within a row.
            if (tid == i_text)
            {
                const int fromLeftScore = _thisScores[tid - 1] - gapPenalty;
                const bool isPrevGreater = (maxFromPrev > fromLeftScore);

                _thisScores[tid] = isPrevGreater ? maxFromPrev : fromLeftScore;
                thisRowM[i_text] = isPrevGreater ? tmpDir : DIRECTION::LEFT;
            }
            __syncthreads();
        }

        thisRowM += numCols;
    }

    if (tid == (numCols - 1))
        finalScore[0] = _thisScores[tid];
}

int old_alignSequenceGPU(const SequenceAlignment::Request &request,
                               SequenceAlignment::Response *response)
{
    const uint64_t numRows = request.textNumBytes + 1;
    const uint64_t numCols = request.patternNumBytes + 1;
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
        return -1;
    }
    /** End Allocate host memory */

    char *d_textBytes, *d_patternBytes, *d_M;
    int *d_scoreMatrix, *d_finalScore;

    auto freeMemory = [&]()
    {
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        cudaFree(d_M);
        cudaFree(d_scoreMatrix);
        cudaFree(d_finalScore);
    };

    /** Allocate and transfer memory to device */
    if (cudaMalloc(&d_scoreMatrix, sizeof(int) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_finalScore, sizeof(int)) != cudaSuccess)
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

    #ifdef BENCHMARK
        auto begin = omp_get_wtime();
    #endif

    const unsigned int sharedMemSize = 2 * numCols * sizeof(int);
    cuda_fillMatrixNW_horizontal<<<1, numCols, sharedMemSize>>>(d_patternBytes, d_textBytes,
                                                        d_scoreMatrix, request.alphabetSize,
                                                        request.gapPenalty, numRows, numCols,
                                                        d_M, d_finalScore);


    if (cudaMemcpy(M, d_M, numRows*numCols, cudaMemcpyDeviceToHost) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        freeMemory();
        return -1;
    }
    freeMemory();

    #ifdef BENCHMARK
        auto end = omp_get_wtime();
        return (end - begin) * 1000.0;
    #endif


    return 0;
}
