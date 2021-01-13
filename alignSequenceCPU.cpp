#include "SequenceAlignment.hpp"

#include <iostream>

void SequenceAlignment::traceBack()
{
    alignmentNumBytes = 0;
    int traceIdx = patternNumBytes * textNumBytes;

    while (traceIdx > 1)
    {
        // Get back actual letter.
        textBytes[alignmentNumBytes] = alphabet[alignMatrix[traceIdx].letter];
        traceIdx = alignMatrix[traceIdx].prev;
        ++alignmentNumBytes;
    }

    std::reverse(textBytes, textBytes + alignmentNumBytes);
}

void SequenceAlignment::alignSequenceCPU()
{
    const char gapByte = alphabetSize - 1;
    const short gapScore = scoreMatrix[alphabetSize * alphabetSize - 1];

    for (int i_text = 1; i_text <= textNumBytes; ++i_text)
    {
        for (int i_pattern = 1; i_pattern <= patternNumBytes; ++i_pattern)
        {
            // Compute this and neigbour indexes.
            const unsigned int thisAlignIdx = i_text * patternNumBytes + i_pattern;
            const unsigned int leftAlignIdx = thisAlignIdx - 1;
            const unsigned int topAlignIdx = thisAlignIdx - textNumBytes;
            const unsigned int diagonalAlignIdx = thisAlignIdx - textNumBytes - 1;

            // Get score for this letter combination.
            const char textByte = textBytes[i_text - 1];
            const char patternByte = patternBytes[i_pattern - 1];
            const int scoreMatrixIdx = ((int)textByte) * SequenceAlignment::alphabetSize + ((int)patternByte);

            // Calculate all neoghbour alignment scores.
            const int fromDiagonalScore = alignMatrix[diagonalAlignIdx].val + scoreMatrix[scoreMatrixIdx];
            const int fromAboveScore = alignMatrix[topAlignIdx].val + gapScore +
                                       gapExtend * (alignMatrix[topAlignIdx].gapLen + 1);
            const int fromLeftScore = alignMatrix[leftAlignIdx].val + gapScore +
                                      gapExtend * (alignMatrix[leftAlignIdx].gapLen + 1);

            // Find out the best alignment.
            const bool isFromDiagonal = fromDiagonalScore >= std::max(fromLeftScore, fromAboveScore);
            const bool isFromLeft = !(isFromDiagonal) && fromLeftScore >= fromAboveScore;
            const bool isFromAbove = !(isFromDiagonal) && !(isFromLeft);

            // Populate this alignPoint with the best alignment.
            alignMatrix[thisAlignIdx].val = std::max(fromDiagonalScore, std::max(fromLeftScore, fromAboveScore));
            alignMatrix[thisAlignIdx].prev = isFromDiagonal * diagonalAlignIdx +
                                             isFromLeft * leftAlignIdx +
                                             isFromAbove * leftAlignIdx;
            // gapLen will be 0 if we align from the diagonal.
            alignMatrix[thisAlignIdx].gapLen = isFromLeft * (alignMatrix[fromLeftScore].gapLen + 1) +
                                               isFromAbove * (alignMatrix[fromAboveScore].gapLen + 1);

            alignMatrix[thisAlignIdx].letter = (!isFromDiagonal) * gapByte + isFromDiagonal * textByte;
        }
    }

    traceBack();
}
