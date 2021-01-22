#include <cuda.h>

#include "SequenceAlignment.hpp"



__global__ void alignSequenceGlobalCUDA(const char *textBytes, const uint64_t textNumBytes,
                                        const char *patternBytes, const uint64_t patternNumBytes,
                                        const char *alphabet, const int alphabetSize,
                                        const short *scoreMatrix, const short gapPenalty,
                                        const int numRows, const int numCols,
                                        int *thisScores, int *prevScores)
{
    extern __shared__ int _shared[];
    int *_text = _shared;
    int *_thisScores = _shared + textNumBytes;
    int *_prevScores = _shared + textNumBytes + numCols;

    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Each thread copies one text letter, i.e. one column.
    _text[tid] = textBytes[tid];
    const char textByte = _text[tid];
    // Init first row.
    _thisScores[tid] = tid * -gapPenalty;
    __syncthreads();

    // Dynamic programming loop.
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Advance one row.
        int *tmp = _thisScores;
        _thisScores = _prevScores;
        _prevScores = tmp;

        if (tid == 0)
            _thisScores[tid] -= gapPenalty;

        const char patternByte = patternBytes[i_pattern];
        const int scoreMatrixIdx = ((int) patternByte) * alphabetSize + ((int) textByte);

        const int fromTopScore = _prevScores[tid] - gapPenalty;
        const int fromDiagScore = _prevScores[tid - 1] + scoreMatrix[scoreMatrixIdx];

        const int maxFromPrev = max(fromDiagScore, fromTopScore);

        for (int i_pattern2 = 1; i_pattern2 < numRows; ++i_pattern2)
        {
            if (tid == i_pattern2)
            {
                const int fromLeftScore = _thisScores[tid] - gapPenalty;
                _thisScores[tid] = max(maxFromPrev, fromLeftScore);
            }
            __syncthreads();
        }
    }


    thisScores[tid] = _thisScores[tid];
    __syncthreads();
}


void SequenceAlignment::alignSequenceGlobalGPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{
    const unsigned int numCols = request.textNumBytes + 1;
    const unsigned int numRows = request.patternNumBytes + 1;

    int *d_thisScores, *d_prevScores;
    short *d_scoreMatrix;
    char *d_textBytes, *d_patternBytes, *d_alphabet;

    if (cudaMalloc(&d_thisScores, sizeof(int) * numCols) != cudaSuccess ||
        cudaMalloc(&d_prevScores, sizeof(int) * numCols) != cudaSuccess ||
        cudaMalloc(&d_scoreMatrix, sizeof(short) * request.alphabetSize * request.alphabetSize) != cudaSuccess ||
        cudaMalloc(&d_textBytes, request.textNumBytes) != cudaSuccess ||
        cudaMalloc(&d_patternBytes, request.patternNumBytes) != cudaSuccess ||
        cudaMalloc(&d_alphabet, request.alphabetSize) != cudaSuccess)
    {
        std::cout << MEM_ERROR << std::endl;
        cudaFree(d_thisScores);
        cudaFree(d_prevScores);
        cudaFree(d_scoreMatrix);
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
        cudaFree(d_textBytes);
        cudaFree(d_patternBytes);
        return;
    }

    const unsigned int sharedMemSize = request.textNumBytes +   // text bytes buffer
                                       sizeof(int) * numCols +  // this row scores
                                       sizeof(int) * numCols;   // prev row scores

    alignSequenceGlobalCUDA<<<1, numCols, sharedMemSize>>>(d_textBytes, request.textNumBytes,
                                                           d_patternBytes, request.patternNumBytes,
                                                           d_alphabet, request.alphabetSize,
                                                           d_scoreMatrix, request.gapPenalty,
                                                           numRows, numCols,
                                                           d_thisScores, d_prevScores);

    if (cudaMemcpy(&(response->score), (d_thisScores + numCols - 1), sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
        std::cout << "Could not copy back to host memory" << std::endl;

    cudaFree(d_thisScores);
    cudaFree(d_prevScores);
    cudaFree(d_scoreMatrix);
    cudaFree(d_textBytes);
    cudaFree(d_patternBytes);

    std::cout << "# Score: \t" << response->score << "\n";

}
