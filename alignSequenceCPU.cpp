#include "SequenceAlignment.hpp"

#include <iostream>
#include <iterator>


void SequenceAlignment::traceBack(const alignPoint *alignMatrix,
                                  const unsigned int numRows, const unsigned int numCols,
                                  const SequenceAlignment::Request &request,
                                  SequenceAlignment::Response *response)
{
    const int lastPointIdx = numRows * numCols - 1;
    int curr = lastPointIdx;

    int maxScore = alignMatrix[curr].score;
    int numEndGapsText = 0;
    int numEndGapsPattern = 0;

    // Check last row and col to see if we can get a better score by staring a gap at the end.
    for (int i = ((numRows-1) * numCols); i < (numRows*numCols); ++i)
    {
        int scoreFromHere = alignMatrix[i].score + request.gapOpenScore +
                            (request.gapExtendScore * (lastPointIdx - i - 1));
        maxScore = (scoreFromHere > maxScore) ? scoreFromHere : maxScore;
        curr = (scoreFromHere > maxScore) ? i : curr;
    }
    // Check last column.
    for (int i = 0; i < numRows; ++i)
    {
        const int lastCol = i*numCols + (numCols-1);
        int scoreFromHere = alignMatrix[lastCol].score + request.gapOpenScore +
                            (request.gapExtendScore * (numRows - i - 1));
        maxScore = (scoreFromHere > maxScore) ? scoreFromHere : maxScore;
        curr = (scoreFromHere > maxScore) ? lastCol : curr;
    }

    const int startOfLastRow = (numRows-1) * numCols;
    numEndGapsText = std::max(0, (curr - startOfLastRow) - (lastPointIdx - startOfLastRow));
    numEndGapsPattern = (curr < startOfLastRow) * ((int) std::floor((lastPointIdx - curr) / numCols));

    std::fill_n(response->alignedTextBytes, numEndGapsText, request.alphabetSize);
    std::fill_n(response->alignedPatternBytes, numEndGapsPattern, request.alphabetSize);

    auto endText = std::make_reverse_iterator(request.textBytes + request.textNumBytes);
    std::copy_n(endText, numEndGapsPattern, (response->alignedTextBytes + numEndGapsText));
    auto endPattern = std::make_reverse_iterator(request.patternBytes + request.patternNumBytes);
    std::copy_n(endPattern, numEndGapsText, (response->alignedPatternBytes + numEndGapsPattern));

    int textIndex = request.textNumBytes - 1 - numEndGapsPattern;
    int patternIndex = request.patternNumBytes - 1 - numEndGapsText;

    response->numAlignmentBytes = std::max(numEndGapsPattern, numEndGapsText);
    response->score = maxScore;

    while (curr != 0)
    {
        const bool takeText = (alignMatrix[curr].isFromDiag || alignMatrix[curr].isFromTop);
        const bool takePattern = (alignMatrix[curr].isFromDiag || alignMatrix[curr].isFromLeft);
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

    std::reverse(response->alignedTextBytes, (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes, (response->alignedPatternBytes + response->numAlignmentBytes));

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
        alignMatrix[i_text].score = ((i_text - 1) * request.gapExtendScore + request.gapOpenScore);
        alignMatrix[i_text].isFromLeft = true;
        alignMatrix[i_text].isFromDiag = false;
        alignMatrix[i_text].isFromTop = false;
    }
    for (int i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        alignMatrix[i_pattern * numCols].score = ((i_pattern - 1) * request.gapExtendScore +
                                                  request.gapOpenScore);
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
            // std::cout << "--------\n[" << i_text << ", " << i_pattern << "] " << request.alphabet[request.textBytes[i_text - 1]] << "-" << request.alphabet[request.patternBytes[i_pattern - 1]] << "\t";
            // Get score for this letter combination. Note that i_text and i_pattern point one
            // beyond the actual text and pattern becuase of the gap character at the beginning.
            const char textByte = request.textBytes[i_text - 1];
            const char patternByte = request.patternBytes[i_pattern - 1];
            const int scoreMatrixIdx = ((int)textByte) * request.alphabetSize + ((int)patternByte);

            // Calculate all alignment scores.
            const int fromDiagonalScore = alignMatrix[diagonalAlignIdx].score +
                                          request.scoreMatrix[scoreMatrixIdx];
            // When selecting the top or left path, we are either opening a new gap or
            // extending an existing path.
            const int fromTopScore = alignMatrix[topAlignIdx].score +
                                     request.gapOpenScore * !(alignMatrix[topAlignIdx].isFromTop) +
                                     request.gapExtendScore * alignMatrix[topAlignIdx].isFromTop;
            const int fromLeftScore = alignMatrix[leftAlignIdx].score +
                                      request.gapOpenScore * !(alignMatrix[leftAlignIdx].isFromLeft) +
                                      request.gapExtendScore * alignMatrix[leftAlignIdx].isFromLeft;

            // Find out the best alignment.
            // Order in case of ties: left, top, diagonal.
            const int maxWithGap = std::max(fromLeftScore, fromTopScore);
            const int bestScore = std::max(fromDiagonalScore, maxWithGap);
            const bool isFromDiagonal = (fromDiagonalScore > maxWithGap);
            const bool isFromLeft = (fromDiagonalScore <= maxWithGap) && (fromLeftScore >= fromTopScore);
            const bool isFromTop = (fromDiagonalScore <= maxWithGap) && (fromLeftScore < fromTopScore);

            // std::cout << "Left: " << fromLeftScore << "  |  " << "Diag: " << fromDiagonalScore << "  |  " << "Top: " << fromTopScore << "\t";

            // Populate this alignPoint with the best alignment.
            alignMatrix[thisAlignIdx].score = bestScore;
            alignMatrix[thisAlignIdx].isFromLeft = isFromLeft;
            alignMatrix[thisAlignIdx].isFromDiag = isFromDiagonal;
            alignMatrix[thisAlignIdx].isFromTop = isFromTop;

            // if (isFromTop)
            //     std::cout << "-> top: " << fromTopScore;
            // if (isFromLeft)
            //     std::cout << "-> left " << fromLeftScore;
            // if (isFromDiagonal)
            //     std::cout << "-> diag " << fromDiagonalScore;
            // std::cout << "\n";

            ++leftAlignIdx; ++thisAlignIdx; ++topAlignIdx; ++diagonalAlignIdx;
        }
    }

    traceBack(alignMatrix, numRows, numCols, request, response);

}
