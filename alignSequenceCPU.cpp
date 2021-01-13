#include "SequenceAlignment.hpp"

#include <iostream>


// void traceBack()
// {
//     alignmentNumBytes = 0;
//     int traceIdx = patternNumBytes * textNumBytes;

//     while (traceIdx > 1)
//     {
//         // Get back actual letter.
//         textBytes[alignmentNumBytes] = alphabet[alignMatrix[traceIdx].letter];
//         traceIdx = alignMatrix[traceIdx].prev;
//         ++alignmentNumBytes;
//     }

//     std::reverse(textBytes, textBytes + alignmentNumBytes);
// }

void SequenceAlignment::alignSequenceCPU(const SequenceAlignment::Request &request,
                                         SequenceAlignment::Response *response)
{

        std::cout << "test2\n";
    enum direction {LEFT, DIAG, TOP};
    struct alignPoint {int val; unsigned int prev; unsigned int gapLen; direction dir;};

    /// Buffer holding the values of the alignment matrix and a trace. Zeroed at start.
    const unsigned int matrixSize = MAX_SEQUENCE_LEN * MAX_SEQUENCE_LEN + 2 * MAX_SEQUENCE_LEN;
    alignPoint alignMatrix[matrixSize] = {};

    std::cout << "test\n";

    for (int i_text = 1; i_text <= request.textNumBytes; ++i_text)
    {
        for (int i_pattern = 1; i_pattern <= request.patternNumBytes; ++i_pattern)
        {
            // Compute this and neigbour indexes.
            const unsigned int thisAlignIdx = i_text * request.patternNumBytes + i_pattern;
            const unsigned int leftAlignIdx = thisAlignIdx - 1;
            const unsigned int topAlignIdx = thisAlignIdx - request.textNumBytes;
            const unsigned int diagonalAlignIdx = thisAlignIdx - request.textNumBytes - 1;

            // Get score for this letter combination.
            const char textByte = request.textBytes[i_text - 1];
            const char patternByte = request.patternBytes[i_pattern - 1];
            const int scoreMatrixIdx = ((int)textByte) * request.alphabetSize + ((int)patternByte);

            // Calculate all neoghbour alignment scores.
            const int fromDiagonalScore = alignMatrix[diagonalAlignIdx].val + request.scoreMatrix[scoreMatrixIdx];
            const int fromTopScore = alignMatrix[topAlignIdx].val + request.gapOpen +
                                       request.gapExtend * (alignMatrix[topAlignIdx].gapLen + 1);
            const int fromLeftScore = alignMatrix[leftAlignIdx].val + request.gapOpen +
                                      request.gapExtend * (alignMatrix[leftAlignIdx].gapLen + 1);

            // Find out the best alignment.
            const bool isFromDiagonal = fromDiagonalScore >= std::max(fromLeftScore, fromTopScore);
            const bool isFromLeft = !(isFromDiagonal) && fromLeftScore >= fromTopScore;
            const bool isFromTop = !(isFromDiagonal) && !(isFromLeft);

            // Populate this alignPoint with the best alignment.
            alignMatrix[thisAlignIdx].val = std::max(fromDiagonalScore, std::max(fromLeftScore, fromTopScore));
            alignMatrix[thisAlignIdx].prev = isFromDiagonal * diagonalAlignIdx +
                                             isFromLeft * leftAlignIdx +
                                             isFromTop * leftAlignIdx;
            // gapLen will be 0 if we align from the diagonal.
            alignMatrix[thisAlignIdx].gapLen = isFromLeft * (alignMatrix[fromLeftScore].gapLen + 1) +
                                               isFromTop * (alignMatrix[fromTopScore].gapLen + 1);

            alignMatrix[thisAlignIdx].dir = (direction) (isFromLeft * direction::LEFT +
                                                         isFromDiagonal * direction::DIAG +
                                                         isFromTop * direction::TOP);
        }
    }

    // // traceBack();
}
