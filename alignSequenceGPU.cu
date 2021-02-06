#include <cuda.h>

#include "SequenceAlignment.hpp"


enum DIR { LEFT, DIAG, TOP};


__global__ void cuda_fillMatrixNW(const char *textBytes, const uint64_t textNumBytes,
                                  const char *patternBytes, const uint64_t patternNumBytes,
                                  const short *scoreMatrix, const int alphabetSize,
                                  const short gapPenalty, const int numRows, const int numCols,
                                  char *M, int *finalScore)
{
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

__global__ void cuda_fillMatrixDiagonalNW(const char *textBytes, const uint64_t textNumBytes,
                                          const char *patternBytes, const uint64_t patternNumBytes,
                                          const short *scoreMatrix, const int alphabetSize,
                                          const short gapPenalty, const int numRows, const int numCols,
                                          char *M, int *finalScore)
{
    extern __shared__ int _shared[];
    int *_thisScores = _shared;
    int *_prevScores = _shared + numRows;
    int *_prevPrevScores = _shared + numRows*2;

    const int tid = threadIdx.x;

    const char patternByte = (tid > 0 && tid < numRows)
                             ? patternBytes[max(numRows - 1, tid - 1)]
                             : alphabetSize;

    // First Half of matrix filling
    int diag_size = 0;
    for (int i_text = 0; i_text < numCols; ++i_text)
    {
        // Advance one diag.
        auto tmp = _prevScores;
        _prevScores = _thisScores;
        _thisScores = _prevPrevScores;
        _prevPrevScores = tmp;

        if (tid == 0)
        {
            _thisScores[tid] = -(i_text * gapPenalty);
            M[i_text] = DIR::LEFT;
        }
        else if (tid == (diag_size-1))
        {
            _thisScores[tid] = -(tid * gapPenalty);
            M[tid * numCols] = DIR::TOP;
        }
        else if (tid < diag_size)
        {
            const char textByte = textBytes[i_text - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            const int fromLeftScore = _prevScores[tid] - gapPenalty;
            const int fromTopScore = _prevScores[tid - 1] - gapPenalty;
            const int fromDiagScore = _prevPrevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

            const int maxWithGap = max(fromLeftScore, fromTopScore);
            const int maxScore = max(maxWithGap, fromDiagScore);
            _thisScores[tid] = maxScore;

            const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
            M[tid*numCols + i_text] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
        }

        diag_size = min(diag_size+1, numRows);

        __syncthreads();
    }

    // Second half of matrix filling
    const char textByte = (tid > 0 && tid < numCols)
                           ? textBytes[min(numCols - tid, 0)]
                           : alphabetSize;
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Advance one diag.
        auto tmp = _prevScores;
        _prevScores = _thisScores;
        _thisScores = _prevPrevScores;
        _prevPrevScores = tmp;

        if (tid > 1 && tid < diag_size)
        {
            const char patternByte = patternBytes[tid + i_pattern - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            const int fromLeftScore = _prevScores[tid] - gapPenalty;
            const int fromTopScore = _prevScores[tid - 1] - gapPenalty;
            const int fromDiagScore = _prevPrevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

            const int maxWithGap = max(fromLeftScore, fromTopScore);
            const int maxScore = max(maxWithGap, fromDiagScore);
            _thisScores[tid - 1] = maxScore;

            const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
            M[(tid+i_pattern)*numCols - tid] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
        }

        --diag_size;

        __syncthreads();
    }

    if (tid == 0)
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

    int *d_finalScore;
    short *d_scoreMatrix;
    char *d_textBytes, *d_patternBytes, *d_M;

    auto freeMemory = [&]()
    {
        cudaFree(d_finalScore);
        cudaFree(d_scoreMatrix);
        cudaFree(d_M);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        delete [] M;
    };

    /** Allocate and transfer memory to device */
    if (cudaMalloc(&d_finalScore, sizeof(int)) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(short) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        freeMemory();
        return;
    }
    if (cudaMemcpy(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_scoreMatrix, request.scoreMatrix, sizeof(short) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        freeMemory();
        return;
    }
    /** End Allocate and transfer memory */

    const unsigned int sharedMemSize = 2 * numCols * sizeof(int);
    // std::cout << "Num col: " << numCols << "\n";
    // std::cout << "Num bytes in shared: " << sharedMemSize << "\n";
    // std::cout << "Num cols: " << numCols << "\n";

    cuda_fillMatrixNW<<<1, numCols, sharedMemSize>>>(d_textBytes, request.textNumBytes,
                                                     d_patternBytes, request.patternNumBytes,
                                                     d_scoreMatrix, request.alphabetSize,
                                                     request.gapPenalty, numRows, numCols,
                                                     d_M, d_finalScore);

    if (cudaMemcpy(&(response->score), (d_finalScore), sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess ||
        cudaMemcpy(M, d_M, numRows*numCols, cudaMemcpyDeviceToHost) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        return;
    }

    traceBack(M, numRows, numCols, request, response);

    freeMemory();

    // std::cout << "# Score: \t" << response->score << "\n";

}
