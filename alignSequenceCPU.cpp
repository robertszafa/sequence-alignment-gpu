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
    // Check last row;
    for (int i = ((numRows-1) * numCols); i < (numRows*numCols); ++i)
    {
        bool isLarger = alignMatrix[i].score > maxScore;
        maxScore = isLarger ? alignMatrix[i].score : maxScore;
        curr = isLarger ? i : curr;
    }
    // Check last column.
    for (int i = 0; i < numRows; ++i)
    {
        const int lastCol = i*numCols + (numCols-1);
        bool isLarger = alignMatrix[lastCol].score > maxScore;
        maxScore = isLarger ? alignMatrix[lastCol].score : maxScore;
        curr = isLarger ? lastCol : curr;
    }

    const int startOfLastRow = (numRows-1) * numCols;
    const int numEndGapsText = std::max(0, (curr - startOfLastRow) - (lastPointIdx - startOfLastRow));
    const int numEndGapsPattern = (curr < startOfLastRow) * ((int) std::floor((lastPointIdx - curr) / numCols));

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

void SequenceAlignment::alignSequenceCPU(const SequenceAlignment::Request &request,
                                         SequenceAlignment::Response *response)
{

    /// Buffer holding the values of the alignment matrix and a trace. Zeroed at start.
    /// Aditional row and column for the gap character.
    const unsigned int numRows = request.textNumBytes + 1;
    const unsigned int numCols = request.patternNumBytes + 1;
    static alignPoint alignMatrix[(MAX_SEQUENCE_LEN+1) * (MAX_SEQUENCE_LEN+1)];

    // Init first row and column.
    alignMatrix[0].score = 0;
    alignMatrix[0].isFromLeft = false;
    alignMatrix[0].isFromTop = false;
    alignMatrix[0].isFromDiag = false;
    for (int i_pattern = 1; i_pattern < numCols; ++i_pattern)
    {
        alignMatrix[i_pattern].score = (i_pattern - 1) * request.gapExtendScore + request.gapOpenScore;
        alignMatrix[i_pattern].isFromLeft = true;
        alignMatrix[i_pattern].isFromDiag = false;
        alignMatrix[i_pattern].isFromTop = false;
    }
    for (int i_text = 1; i_text < numRows; ++i_text)
    {
        alignMatrix[i_text * numCols].score = (i_text - 1) * request.gapExtendScore + request.gapOpenScore;
        alignMatrix[i_text * numCols].isFromLeft = false;
        alignMatrix[i_text * numCols].isFromDiag = false;
        alignMatrix[i_text * numCols].isFromTop = true;
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
            // Also check if we are at the end of the pattern or text sequence.
            // If yes, there is no penalty for gaps at the end.
            // In bioinformatics, it is usually reasonable to assume that the sequences are
            // incomplete and there should be no penalty for failing to align the missing bases.
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
