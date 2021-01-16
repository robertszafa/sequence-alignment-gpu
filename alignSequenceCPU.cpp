#include "SequenceAlignment.hpp"

#include <iostream>
#include <iterator>


void SequenceAlignment::traceBack(const alignPoint *alignMatrix,
                                  const unsigned int numRows, const unsigned int numCols,
                                  const SequenceAlignment::Request &request,
                                  SequenceAlignment::Response *response)
{
    int curr = numRows * numCols - 1;
    int textIndex = request.textNumBytes - 1;
    int patternIndex = request.patternNumBytes - 1;

    response->numAlignmentBytes = 0;
    response->score = alignMatrix[curr].score;

    while (curr != 0)
    {
        const bool takeText = (alignMatrix[curr].isFromDiag || alignMatrix[curr].isFromLeft);
        const bool takePattern = (alignMatrix[curr].isFromDiag || alignMatrix[curr].isFromTop);

        response->alignedTextBytes[response->numAlignmentBytes] =
            takeText * request.textBytes[textIndex] + (!takeText) * request.alphabetSize;
        response->alignedPatternBytes[response->numAlignmentBytes] =
            takePattern * request.patternBytes[patternIndex] + (!takePattern) * request.alphabetSize;

        response->numAlignmentBytes += 1;

        // Saturate at 0.
        textIndex = std::max(0, textIndex - takeText);
        patternIndex = std::max(0, patternIndex - takePattern);
        curr -= (alignMatrix[curr].isFromLeft) +
                (alignMatrix[curr].isFromDiag * (numCols+1)) +
                (alignMatrix[curr].isFromTop * (numCols));
    }

    // Reverse and tranform from alphabet indexes to alphabet elements.
    std::reverse(response->alignedTextBytes,
                 (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes,
                 (response->alignedPatternBytes + response->numAlignmentBytes));
    for (int i=0; i<response->numAlignmentBytes; ++i)
    {
        response->alignedTextBytes[i] = request.alphabet[response->alignedTextBytes[i]];
        response->alignedPatternBytes[i] = request.alphabet[response->alignedPatternBytes[i]];
    }
}


void SequenceAlignment::alignSequenceGlobalCPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{

    /// Buffer holding the values of the alignment matrix and a trace. Zeroed at start.
    // TODO: Store two rows in memory. When finished, swap pointers and write to disk.
    static alignPoint alignMatrix[(MAX_SEQUENCE_LEN+1) * (MAX_SEQUENCE_LEN+1)];
    /// Aditional row and column for the gap character.
    const unsigned int numCols = request.textNumBytes + 1;
    const unsigned int numRows = request.patternNumBytes + 1;

    // Init first row and column.
    alignMatrix[0].score = 0;
    alignMatrix[0].isFromLeft = alignMatrix[0].isFromTop = alignMatrix[0].isFromDiag = false;
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        alignMatrix[i_text].score = i_text * (-request.gapPenalty);
        alignMatrix[i_text].isFromLeft = true;
        alignMatrix[i_text].isFromDiag = false;
        alignMatrix[i_text].isFromTop = false;
    }
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        alignMatrix[i_pattern * numCols].score = i_pattern * (-request.gapPenalty);
        alignMatrix[i_pattern * numCols].isFromLeft = false;
        alignMatrix[i_pattern * numCols].isFromDiag = false;
        alignMatrix[i_pattern * numCols].isFromTop = true;
    }

    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Compute this and neigbour indexes.
        // TODO: Store two rows in memory. When finished, swap pointers and write to disk.
        unsigned int leftAlignIdx = i_pattern * numCols;
        unsigned int thisAlignIdx = leftAlignIdx + 1;
        unsigned int topAlignIdx = thisAlignIdx - numCols;
        unsigned int diagonalAlignIdx = thisAlignIdx - numCols - 1;

        for (int i_text = 1; i_text < numCols; ++i_text)
        {
            // Get score for this letter combination. Note that i_text and i_pattern point one
            // beyond the actual text and pattern becuase of the gap character at the beginning.
            const char textByte = request.textBytes[i_text - 1];
            const char patternByte = request.patternBytes[i_pattern - 1];
            const int scoreMatrixIdx = ((int) patternByte) * request.alphabetSize + ((int) textByte);

            // Calculate all alignment scores.
            const int fromLeftScore = alignMatrix[leftAlignIdx].score - request.gapPenalty;
            const int fromTopScore = alignMatrix[topAlignIdx].score - request.gapPenalty;
            const int fromDiagonalScore = alignMatrix[diagonalAlignIdx].score +
                                          request.scoreMatrix[scoreMatrixIdx];

            // Find out the best alignment.
            // Order in case of ties: left, top, diagonal.
            const int maxWithGap = std::max(fromLeftScore, fromTopScore);
            const int bestScore = std::max(fromDiagonalScore, maxWithGap);
            const bool isFromDiagonal = (fromDiagonalScore > maxWithGap);
            const bool isFromLeft = (fromDiagonalScore <= maxWithGap) && (fromLeftScore >= fromTopScore);
            const bool isFromTop = (fromDiagonalScore <= maxWithGap) && (fromLeftScore < fromTopScore);

            // Populate this alignPoint with the best alignment.
            alignMatrix[thisAlignIdx].score = bestScore;
            alignMatrix[thisAlignIdx].isFromLeft = isFromLeft;
            alignMatrix[thisAlignIdx].isFromDiag = isFromDiagonal;
            alignMatrix[thisAlignIdx].isFromTop = isFromTop;

            ++leftAlignIdx; ++thisAlignIdx; ++topAlignIdx; ++diagonalAlignIdx;
        }
    }

    traceBack(alignMatrix, numRows, numCols, request, response);

}
