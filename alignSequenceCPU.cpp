#include "SequenceAlignment.hpp"

#include <iostream>


struct alignPoint
{
    int val; unsigned short gapLenTop; unsigned short gapLenLeft;
    bool fromLeft; bool fromDiagonal; bool fromTop;
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
        alignMatrix[i_pattern].fromLeft = true;
        alignMatrix[i_pattern].fromDiagonal = false;
        alignMatrix[i_pattern].fromTop = false;
        alignMatrix[i_pattern].gapLenLeft = i_pattern;
        alignMatrix[i_pattern].gapLenTop = 0;
    }
    for (int i_text = 1; i_text < numRows; ++i_text)
    {
        alignMatrix[i_text * numCols].val = (i_text - 1) * request.gapExtendScore + request.gapOpenScore;
        alignMatrix[i_text * numCols].fromLeft = false;
        alignMatrix[i_text * numCols].fromDiagonal = false;
        alignMatrix[i_text * numCols].fromTop = true;
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
            const int fromDiagonalScore = alignMatrix[diagonalAlignIdx].val +
                                          request.scoreMatrix[scoreMatrixIdx];
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
            alignMatrix[thisAlignIdx].fromLeft = isFromLeft;
            alignMatrix[thisAlignIdx].fromDiagonal = isFromDiagonal;
            alignMatrix[thisAlignIdx].fromTop = isFromTop;

            ++leftAlignIdx; ++thisAlignIdx; ++topAlignIdx; ++diagonalAlignIdx;
        }
    }

    // traceBack
    response->numAlignmentBytes = 0;

    unsigned int curr = numRows*numCols - 1;
    unsigned int textIndex = request.textNumBytes - 1;
    unsigned int patternIndex = request.patternNumBytes - 1;

    while (curr != 0)
    {
        const bool takeText = (alignMatrix[curr].fromDiagonal || alignMatrix[curr].fromTop);
        const bool takePattern = (alignMatrix[curr].fromDiagonal || alignMatrix[curr].fromLeft);
        response->alignedTextBytes[response->numAlignmentBytes] =
            takeText * request.textBytes[textIndex] + (!takeText) * request.alphabetSize;
        response->alignedPatternBytes[response->numAlignmentBytes] =
            takePattern * request.patternBytes[patternIndex] + (!takePattern) * request.alphabetSize;

        response->numAlignmentBytes += 1;
        textIndex -= takeText;
        patternIndex -= takePattern;
        curr -= (alignMatrix[curr].fromLeft) +
                (alignMatrix[curr].fromDiagonal * (numCols+1)) +
                (alignMatrix[curr].fromTop * (numCols));
    }

    std::reverse(response->alignedTextBytes, (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes, (response->alignedPatternBytes + response->numAlignmentBytes));

    for (int i=0; i<response->numAlignmentBytes; ++i)
    {
        response->alignedTextBytes[i] = request.alphabet[response->alignedTextBytes[i]];
        response->alignedPatternBytes[i] = request.alphabet[response->alignedPatternBytes[i]];
    }
}
