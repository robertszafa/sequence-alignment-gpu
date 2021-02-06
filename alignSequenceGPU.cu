#include <cuda.h>

#include "SequenceAlignment.hpp"



__global__ void cuda_fillMatrixNW(const char *textBytes, const char *patternBytes,
                                          const short *scoreMatrix, const int alphabetSize,
                                          const short gapPenalty, const int numRows, const int numCols,
                                          char *M, int *finalScore)
{
    using SequenceAlignment::DIR;

    extern __shared__ int _shared[];
    int *_thisScores = _shared;
    int *_prevScores = _shared + numRows;
    int *_prevPrevScores = _shared + numRows*2;

    const int tid = threadIdx.x;

    const char patternByte = patternBytes[max(0, tid - 1)];

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
            auto inRowIdx = i_text - tid;

            const char textByte = textBytes[inRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            const int fromLeftScore = _prevScores[tid] - gapPenalty;
            const int fromTopScore = _prevScores[tid - 1] - gapPenalty;
            const int fromDiagScore = _prevPrevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

            const int maxWithGap = max(fromLeftScore, fromTopScore);
            _thisScores[tid] = max(maxWithGap, fromDiagScore);

            const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
            M[tid*numCols + inRowIdx] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
        }

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

        if (tid >= i_pattern)
        {
            auto inRowIdx = numCols-1 - tid + i_pattern;

            const char textByte = textBytes[inRowIdx - 1];
            const int scoreMatrixIdx = ((int) textByte) * alphabetSize + ((int) patternByte);

            const int fromLeftScore = _prevScores[tid] - gapPenalty;
            const int fromTopScore = _prevScores[tid - 1] - gapPenalty;
            const int fromDiagScore = _prevPrevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

            const int maxWithGap = max(fromLeftScore, fromTopScore);
            _thisScores[tid] = max(maxWithGap, fromDiagScore);

            const auto dirWithGap = (fromTopScore > fromLeftScore) ? DIR::TOP : DIR::LEFT;
            M[tid*numCols + inRowIdx] = (fromDiagScore > maxWithGap) ? DIR::DIAG : dirWithGap;
        }

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


    const unsigned int sharedMemSize = 3 * numRows * sizeof(int);
    cuda_fillMatrixNW<<<1, numRows, sharedMemSize>>>(d_textBytes, d_patternBytes,
                                                             d_scoreMatrix, request.alphabetSize,
                                                             request.gapPenalty, numRows, numCols,
                                                             d_M, d_finalScore);

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
