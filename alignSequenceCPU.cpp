#include "SequenceAlignment.hpp"

#include <iostream>


enum direction {LEFT, DIAG, TOP};
struct alignPoint
{
    int val; unsigned short gapLenTop; unsigned short gapLenLeft; direction dir;
    /// Order alignment points based on their score.
    bool operator <(const alignPoint &b) const { return val < b.val;};
};

void SequenceAlignment::alignSequenceCPU(const SequenceAlignment::Request &request,
                                         SequenceAlignment::Response *response)
{

    /// Buffer holding the values of the alignment matrix and a trace. Zeroed at start.
    static alignPoint alignMatrix[(MAX_SEQUENCE_LEN+1) * (MAX_SEQUENCE_LEN+1)];
    /// Aditional row and column for the gap character.
    const unsigned int numRows = request.textNumBytes + 1;
    const unsigned int numCols = request.patternNumBytes + 1;

    // Init first row and column.
    alignMatrix[0].val = 0; alignMatrix[0].gapLenLeft = 0; alignMatrix[0].gapLenTop = 0;
    for (int i_pattern = 1; i_pattern < numCols; ++i_pattern)
    {
        alignMatrix[i_pattern].val = (i_pattern - 1) * request.gapExtendScore + request.gapOpenScore;
        alignMatrix[i_pattern].dir = direction::LEFT;
        alignMatrix[i_pattern].gapLenLeft = i_pattern;
        alignMatrix[i_pattern].gapLenTop = 0;
    }
    for (int i_text = 1; i_text < numRows; ++i_text)
    {
        alignMatrix[i_text * numCols].val = (i_text - 1) * request.gapExtendScore + request.gapOpenScore;
        alignMatrix[i_text * numCols].dir = direction::TOP;
        alignMatrix[i_text * numCols].gapLenTop = i_text;
        alignMatrix[i_text * numCols].gapLenLeft = 0;
    }

    for (int i_text = 1; i_text < numRows; ++i_text)
    {
        // Compute this and neigbour indexes.
        unsigned int leftAlignIdx = i_text * numCols;
        unsigned int thisAlignIdx = leftAlignIdx + 1;
        unsigned int topAlignIdx = thisAlignIdx - numCols;
        unsigned int diagonalAlignIdx = thisAlignIdx - numCols - 1;

        for (int i_pattern = 1; i_pattern < numCols; ++i_pattern)
        {
            // Get score for this letter combination. Note that i_text and i_pattern point one
            // beyond the actual text and pattern becuase of the gap character at the beginning.
            const char textByte = request.textBytes[i_text - 1];
            const char patternByte = request.patternBytes[i_pattern - 1];
            const int scoreMatrixIdx = ((int)textByte) * request.alphabetSize + ((int)patternByte);

            // Calculate all alignment scores.
            const int fromDiagonalScore = alignMatrix[diagonalAlignIdx].val + request.scoreMatrix[scoreMatrixIdx];
            // Are we opening a gap?
            const bool gapOpenFromTop = (alignMatrix[topAlignIdx].gapLenTop == 0);
            const bool gapOpneFromLeft = (alignMatrix[leftAlignIdx].gapLenLeft == 0);
            const int fromTopScore = alignMatrix[topAlignIdx].val +
                                     request.gapOpenScore * gapOpenFromTop +
                                     request.gapExtendScore * alignMatrix[topAlignIdx].gapLenTop;
            const int fromLeftScore = alignMatrix[leftAlignIdx].val +
                                     request.gapOpenScore * gapOpneFromLeft +
                                     request.gapExtendScore * alignMatrix[leftAlignIdx].gapLenLeft;

            // Find out the best alignment.
            const bool isFromDiagonal = fromDiagonalScore >= std::max(fromLeftScore, fromTopScore);
            const bool isFromLeft = !(isFromDiagonal) && fromLeftScore >= fromTopScore;
            const bool isFromTop = !(isFromDiagonal) && !(isFromLeft);

            // Populate this alignPoint with the best alignment.
            alignMatrix[thisAlignIdx].val = std::max(fromDiagonalScore, std::max(fromLeftScore, fromTopScore));
            alignMatrix[thisAlignIdx].gapLenLeft = isFromLeft * (alignMatrix[leftAlignIdx].gapLenLeft + 1);
            alignMatrix[thisAlignIdx].gapLenTop = isFromTop * (alignMatrix[topAlignIdx].gapLenTop + 1);
            alignMatrix[thisAlignIdx].dir = (direction) (isFromLeft * direction::LEFT +
                                                         isFromDiagonal * direction::DIAG +
                                                         isFromTop * direction::TOP);

            ++leftAlignIdx; ++thisAlignIdx; ++topAlignIdx; ++diagonalAlignIdx;
        }
    }

    // traceBack
    unsigned int curr = numRows*numCols - 1;
    std::cout << alignMatrix[curr].val << "\n";

    // for (int i=0; i<(std::max(request.textNumBytes, request.patternNumBytes)); ++i)
    // {
    //     // Get back actual letter.
    //     std::cout << curr->dir;
    //     curr = alignMatrix + curr->prev;
    // }

    // std::reverse(textBytes, textBytes + alignmentNumBytes);

}
