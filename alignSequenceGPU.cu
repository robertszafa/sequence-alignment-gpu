#include <cuda.h>

#include "SequenceAlignment.hpp"


enum DIR { LEFT, DIAG, TOP};


__global__ void alignSequenceGlobalCUDA(const char *textBytes, const uint64_t textNumBytes,
                                        const char *patternBytes, const uint64_t patternNumBytes,
                                        const char *alphabet, const int alphabetSize,
                                        const short *scoreMatrix, const short gapPenalty,
                                        const int numRows, const int numCols,
                                        int *thisScores, int *prevScores, char *M)
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
        __syncthreads();
        const int fromTopScore = _prevScores[tid] - gapPenalty;
        const int fromDiagScore = _prevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

        const int maxFromPrev = max(fromDiagScore, fromTopScore);

        for (int i_text = 1; i_text < numCols; ++i_text)
        {
            // We are accessing the previous column within a row.
            if (tid == i_text)
            {
                const int fromLeftScore = _thisScores[tid - 1] - gapPenalty;
                _thisScores[tid] = max(maxFromPrev, fromLeftScore);

                const bool isFromLeft = (fromLeftScore >= maxFromPrev);
                const bool isFromTop = (fromTopScore > fromLeftScore) && (fromTopScore >= fromDiagScore);
                const bool isFromDiag = (fromDiagScore > fromLeftScore) && (fromDiagScore > fromTopScore) ;

                thisRowM[i_text] = (isFromLeft * DIR::LEFT + isFromDiag * DIR::DIAG + isFromTop * DIR::TOP);
            }
            __syncthreads();
        }

        thisRowM += numCols;
    }

    if (tid == (numCols - 1))
        thisScores[tid] = _thisScores[tid];
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

    int *d_thisScores, *d_prevScores;
    short *d_scoreMatrix;
    char *d_textBytes, *d_patternBytes, *d_alphabet;
    char *d_M;

    /** Allocate and transfer memory */
    if (cudaMalloc(&d_thisScores, sizeof(int) * numCols) != cudaSuccess ||
        cudaMalloc(&d_prevScores, sizeof(int) * numCols) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(short) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_M, numRows * numCols) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_alphabet, request.alphabetSize) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cudaFree(d_thisScores);
        cudaFree(d_prevScores);
        cudaFree(d_scoreMatrix);
        cudaFree(d_M);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return;
    }

    if (cudaMemcpy(d_textBytes, request.textBytes, request.textNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_patternBytes, request.patternBytes, request.patternNumBytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_scoreMatrix, request.scoreMatrix, sizeof(short) * (request.alphabetSize * request.alphabetSize), cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_alphabet, request.alphabet, request.alphabetSize, cudaMemcpyHostToDevice) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cudaFree(d_thisScores);
        cudaFree(d_prevScores);
        cudaFree(d_scoreMatrix);
        cudaFree(d_M);
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return;
    }
    /** End Allocate and transfer memory */

    // this and prev row scores
    const unsigned int sharedMemSize = 2 * sizeof(int) * numCols;
    std::cout << "Num bytes in col: " << sharedMemSize << "\n";
    std::cout << "Num cols: " << numCols << "\n";

    alignSequenceGlobalCUDA<<<1, numCols, sharedMemSize>>>(d_textBytes, request.textNumBytes,
                                                           d_patternBytes, request.patternNumBytes,
                                                           d_alphabet, request.alphabetSize,
                                                           d_scoreMatrix, request.gapPenalty,
                                                           numRows, numCols,
                                                           d_thisScores, d_prevScores, d_M);

    if (cudaMemcpy(&(response->score), (d_thisScores + numCols - 1), sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess ||
        cudaMemcpy(M, d_M, numRows*numCols, cudaMemcpyDeviceToHost) != cudaSuccess)
    {
        std::cout << "Could not copy back to host memory" << std::endl;
        return;
    }

    traceBack(M, numRows, numCols, request, response);

    cudaFree(d_thisScores);
    cudaFree(d_prevScores);
    cudaFree(d_scoreMatrix);
    cudaFree(d_M);
    cudaFree(d_textBytes);
    cudaFree(d_patternBytes);
    delete [] M;

    std::cout << "# Score: \t" << response->score << "\n";

}
