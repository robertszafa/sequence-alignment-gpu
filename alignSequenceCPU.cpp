#include "SequenceAlignment.hpp"


void SequenceAlignment::traceBack()
{
    alignmentNumBytes = 0;



}

void SequenceAlignment::alignSequenceCPU()
{
    const char gapByte = alphabetSize * alphabetSize - 1;
    const short gapScore = scoreMatrix[gapByte];

    for (int i_text = 1; i_text <= textNumBytes; ++i_text)
    {
        for (int i_pattern = 1; i_pattern <= patternNumBytes; ++i_pattern)
        {
            const unsigned int thisAlignIdx = i_text * textNumBytes + i_pattern;
            const unsigned int leftAlignIdx = thisAlignIdx - 1;
            const unsigned int topAlignIdx = thisAlignIdx - textNumBytes;
            const unsigned int diagonalAlignIdx = thisAlignIdx - textNumBytes - 1;

            // Get score for this letter combination.
            const char textByte = textBytes[i_text - 1];
            const char patternByte = patternBytes[i_pattern - 1];
            const int scoreMatrixIdx = ((int)textByte) * SequenceAlignment::alphabetSize + ((int)patternByte);

            // Calculate 3 possible alignment scores.
            const int equalScore = alignMatrix[diagonalAlignIdx].val + scoreMatrix[scoreMatrixIdx];
            const int fromAboveScore = alignMatrix[topAlignIdx].val + gapScore +
                                       affineGapScore * (alignMatrix[topAlignIdx].gapLen + 1);
            const int fromLeftScore = alignMatrix[leftAlignIdx].val + gapScore +
                                      affineGapScore * (alignMatrix[leftAlignIdx].gapLen + 1);

            // Select the best alignment and populate this alignPoint struct.
            const bool isFromDiagonal = equalScore >= std::max(fromLeftScore, fromAboveScore);
            const bool isFromLeft = !(isFromDiagonal) && fromLeftScore >= fromAboveScore;
            const bool isFromAbove = !(isFromDiagonal) && !(isFromLeft);

            alignMatrix[thisAlignIdx].val = std::max(equalScore, std::max(fromLeftScore, fromAboveScore));
            alignMatrix[thisAlignIdx].prev = isFromDiagonal * diagonalAlignIdx +
                                             isFromLeft * leftAlignIdx +
                                             isFromAbove * leftAlignIdx;
            // gapLen will be 0 if we align from the diagonal.
            alignMatrix[thisAlignIdx].gapLen = isFromLeft * (alignMatrix[fromLeftScore].gapLen + 1) +
                                               isFromAbove * (alignMatrix[fromAboveScore].gapLen + 1);

            alignMatrix[thisAlignIdx].letter = (!isFromDiagonal) * gapByte + isFromDiagonal * textByte;
        }
    }
}
