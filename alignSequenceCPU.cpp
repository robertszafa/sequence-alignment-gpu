#include "SequenceAlignment.hpp"

#include <iostream>
#include <iterator>


/// A structure representing a cell in the alignment matrix.
struct alignPoint { int score; bool isFromLeft; bool isFromDiag; bool isFromTop; };


void traceBack(const alignPoint *M, const unsigned int numRows, const unsigned int numCols,
               const SequenceAlignment::Request &request, SequenceAlignment::Response *response)
{
    int curr = numRows * numCols - 1;
    int textIndex = request.textNumBytes - 1;
    int patternIndex = request.patternNumBytes - 1;

    response->numAlignmentBytes = 0;
    response->score = M[curr].score;

    while (curr != 0)
    {
        // Was it a match, gap in text or gap in pattern?
        const bool takeText = (M[curr].isFromDiag || M[curr].isFromLeft);
        const bool takePattern = (M[curr].isFromDiag || M[curr].isFromTop);

        // Translate from alphabet indexes to alphabet letters.
        const char translatedText = request.alphabet[request.textBytes[textIndex]];
        const char translatedPattern = request.alphabet[request.patternBytes[patternIndex]];
        const char GAP = request.alphabet[request.alphabetSize];

        response->alignedTextBytes[response->numAlignmentBytes] =
            takeText * translatedText + (!takeText) * GAP;
        response->alignedPatternBytes[response->numAlignmentBytes] =
            takePattern * translatedPattern + (!takePattern) * GAP;

        response->numAlignmentBytes += 1;

        // Update text and pattern indices depending which one was used.
        textIndex = std::max(0, textIndex - takeText);
        patternIndex = std::max(0, patternIndex - takePattern);
        // Go to next cell.
        curr -= (M[curr].isFromLeft) +
                (M[curr].isFromDiag * (numCols+1)) +
                (M[curr].isFromTop * (numCols));
    }

    std::reverse(response->alignedTextBytes,
                 (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes,
                 (response->alignedPatternBytes + response->numAlignmentBytes));
}


void SequenceAlignment::alignSequenceGlobalCPU(const SequenceAlignment::Request &request,
                                               SequenceAlignment::Response *response)
{
    /** Allocate memory */
    /// Buffer holding the values of the alignment matrix and a trace. Zeroed at start.
    /// Aditional row and column for the gap character.
    alignPoint *M;
    const unsigned int numCols = request.textNumBytes + 1;
    const unsigned int numRows = request.patternNumBytes + 1;
    try
    {
        M = new alignPoint[numRows * numCols];
        const int maxNumAlignmentBytes = std::max(request.textNumBytes, request.patternNumBytes) * 2;
        response->alignedTextBytes = new char[maxNumAlignmentBytes];
        response->alignedPatternBytes = new char[maxNumAlignmentBytes];

    }
    catch(const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        return;
    }
    /** End allocate memory */


    // Init first row and column.
    M[0].score = 0;
    M[0].isFromLeft = M[0].isFromTop = M[0].isFromDiag = false;
    for (int i_text = 1; i_text < numCols; ++i_text)
    {
        M[i_text].score = i_text * (-request.gapPenalty);
        M[i_text].isFromLeft = true;
        M[i_text].isFromDiag = false;
        M[i_text].isFromTop = false;
    }
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        M[i_pattern * numCols].score = i_pattern * (-request.gapPenalty);
        M[i_pattern * numCols].isFromLeft = false;
        M[i_pattern * numCols].isFromDiag = false;
        M[i_pattern * numCols].isFromTop = true;
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
            const int fromLeftScore = M[leftAlignIdx].score - request.gapPenalty;
            const int fromTopScore = M[topAlignIdx].score - request.gapPenalty;
            const int fromDiagonalScore = M[diagonalAlignIdx].score +
                                          request.scoreMatrix[scoreMatrixIdx];

            // Find out the best alignment.
            // Order in case of ties: left, top, diagonal.
            const int maxWithGap = std::max(fromLeftScore, fromTopScore);
            const int bestScore = std::max(fromDiagonalScore, maxWithGap);
            const bool isFromDiagonal = (fromDiagonalScore > maxWithGap);
            const bool isFromLeft = (fromDiagonalScore <= maxWithGap) && (fromLeftScore >= fromTopScore);
            const bool isFromTop = (fromDiagonalScore <= maxWithGap) && (fromLeftScore < fromTopScore);

            // Populate this alignPoint with the best alignment.
            M[thisAlignIdx].score = bestScore;
            M[thisAlignIdx].isFromLeft = isFromLeft;
            M[thisAlignIdx].isFromDiag = isFromDiagonal;
            M[thisAlignIdx].isFromTop = isFromTop;

            ++leftAlignIdx; ++thisAlignIdx; ++topAlignIdx; ++diagonalAlignIdx;
        }
    }

    traceBack(M, numRows, numCols, request, response);

}
