#include "SequenceAlignment.hpp"

#include <iostream>
#include <iterator>


using SequenceAlignment::DIRECTION;


void SequenceAlignment::traceBackSW(const char *M, const uint64_t start, const uint64_t numRows,
                                    const uint64_t numCols, const Request &request, Response *response)
{
    int textIndex = (start % numCols) - 1;
    int patternIndex = (start / numCols) - 1;

    response->numAlignmentBytes = 0;

    auto curr = start;
    while (M[curr] != DIRECTION::STOP)
    {
        auto dir = M[curr];

        // Was it a match, gap in text or gap in pattern?
        const bool takeText = (dir == DIRECTION::DIAG || dir == DIRECTION::LEFT);
        const bool takePattern = (dir == DIRECTION::DIAG || dir == DIRECTION::TOP);

        // Translate from alphabet indexes to alphabet letters.
        const char translatedText = request.alphabet[request.textBytes[textIndex]];
        const char translatedPattern = request.alphabet[request.patternBytes[patternIndex]];
        const char GAP = request.alphabet[request.alphabetSize];

        response->alignedTextBytes[response->numAlignmentBytes] =
            takeText * translatedText + (!takeText) * GAP;
        response->alignedPatternBytes[response->numAlignmentBytes] =
            takePattern * translatedPattern + (!takePattern) * GAP;

        response->numAlignmentBytes += 1;

        // Go to next cell.
        curr -= (dir == DIRECTION::LEFT) +
                ((dir == DIRECTION::DIAG) * (numCols+1)) +
                ((dir == DIRECTION::TOP) * (numCols));

        // First row or col?
        if (curr % numCols == 0 || curr < numCols)
            break;

        // Update text and pattern indices depending which one was used.
        if (dir != DIRECTION::STOP)
        {
            textIndex = std::max(0, textIndex - takeText);
            patternIndex = std::max(0, patternIndex - takePattern);
        }
    }

    response->startInAlignedText = textIndex;
    response->startInAlignedPattern = patternIndex;
    std::reverse(response->alignedTextBytes,
                 (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes,
                 (response->alignedPatternBytes + response->numAlignmentBytes));
}

void SequenceAlignment::traceBackNW(const char *M, const uint64_t numRows, const uint64_t numCols,
                                    const Request &request, Response *response)
{
    uint64_t curr = numRows * numCols - 1;
    int textIndex = request.textNumBytes - 1;
    int patternIndex = request.patternNumBytes - 1;

    response->numAlignmentBytes = 0;

    while (curr > 0)
    {
        auto dir = M[curr];

        // First row or col?
        if (curr % numCols == 0)
            dir = DIRECTION::TOP;
        else if (curr < numCols)
            dir = DIRECTION::LEFT;

        // Was it a match, gap in text or gap in pattern?
        const bool takeText = (dir == DIRECTION::DIAG || dir == DIRECTION::LEFT);
        const bool takePattern = (dir == DIRECTION::DIAG || dir == DIRECTION::TOP);

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
        curr -= (dir == DIRECTION::LEFT) +
                ((dir == DIRECTION::DIAG) * (numCols+1)) +
                ((dir == DIRECTION::TOP) * (numCols));
    }

    response->startInAlignedText = textIndex;
    response->startInAlignedPattern = patternIndex;
    std::reverse(response->alignedTextBytes,
                 (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes,
                 (response->alignedPatternBytes + response->numAlignmentBytes));
}

std::pair<int, uint64_t> fillMatrixSW(char *M, const uint64_t numRows, const uint64_t numCols,
                                      const SequenceAlignment::Request &request)
{
    int *thisRowScores = nullptr;
    int *prevRowScores = nullptr;

    auto cleanUp = [&]()
    {
        if (thisRowScores) delete [] thisRowScores;
        if (prevRowScores) delete [] prevRowScores;
        thisRowScores = nullptr;
        prevRowScores = nullptr;
    };

    /** Allocate memory */
    try
    {
        thisRowScores = new int[numCols];
        prevRowScores = new int[numCols];
    }
    catch(const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        cleanUp();
        return {0, 0};
    }
    /** End Allocate memory */

    // Init first row with 0.
    for (uint64_t i_text = 0; i_text < numCols; ++i_text)
    {
        thisRowScores[i_text] = 0;
        M[i_text] = DIRECTION::STOP;
    }

    // Dynamic programming loop
    uint64_t maxIJ = 0;
    int maxScore = 0;
    auto thisRowM = M + numCols;
    for (uint64_t i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Ping-pong score buffers.
        auto tmp = prevRowScores;
        prevRowScores = thisRowScores;
        thisRowScores = tmp;

        // Init first column with 0.
        thisRowScores[0] = 0;
        thisRowM[0] = DIRECTION::STOP;

        for (uint64_t i_text = 1; i_text < numCols; ++i_text)
        {
            // Get score for this letter combination. Note that i_text and i_pattern point one
            // beyond the actual text and pattern becuase of the gap character at the beginning.
            const char textByte = request.textBytes[i_text - 1];
            const char patternByte = request.patternBytes[i_pattern - 1];
            const int scoreIdx = ((int) patternByte) * request.alphabetSize + ((int) textByte);

            // Calculate all alignment scores.
            const int fromLeftScore = thisRowScores[i_text - 1] - request.gapPenalty;
            const int fromTopScore = prevRowScores[i_text] - request.gapPenalty;
            const int fromDiagonalScore = prevRowScores[i_text - 1] + request.scoreMatrix[scoreIdx];

            // Find out the best alignment.
            // Order in case of ties: left, top, diagonal.
            const int maxWithGap = std::max(fromLeftScore, fromTopScore);
            const int bestScore = std::max(fromDiagonalScore, maxWithGap);
            const bool isFromDiag = (fromDiagonalScore > maxWithGap);
            const bool isFromLeft = (fromDiagonalScore <= maxWithGap) && (fromLeftScore >= fromTopScore);
            const bool isFromTop = (fromDiagonalScore <= maxWithGap) && (fromLeftScore < fromTopScore);

            // Populate this alignPoint with the best alignment.
            const auto dirNonZeroScore = isFromLeft * DIRECTION::LEFT + isFromDiag * DIRECTION::DIAG + isFromTop * DIRECTION::TOP;
            thisRowM[i_text] = bestScore > 0 ? dirNonZeroScore : DIRECTION::STOP;
            thisRowScores[i_text] = std::max(0, bestScore);
            maxIJ = thisRowScores[i_text] > maxScore ? (i_pattern * numCols + i_text) : maxIJ;
            maxScore = std::max(maxScore, thisRowScores[i_text]);
        }

        thisRowM += numCols;
    }

    cleanUp();

    return {maxScore, maxIJ};
}

int fillMatrixNW(char *M, const uint64_t numRows, const uint64_t numCols,
                 const SequenceAlignment::Request &request)
{
    int *thisRowScores = nullptr;
    int *prevRowScores = nullptr;

    auto cleanUp = [&]()
    {
        if (thisRowScores) delete [] thisRowScores;
        if (prevRowScores) delete [] prevRowScores;
        thisRowScores = nullptr;
        prevRowScores = nullptr;
    };

    /** Allocate memory */
    try
    {
        thisRowScores = new int[numCols];
        prevRowScores = new int[numCols];
    }
    catch(const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        cleanUp();
        return 0;
    }
    /** End Allocate memory */

    // Init first row.
    for (uint64_t i_text = 0; i_text < numCols; ++i_text)
    {
        thisRowScores[i_text] = i_text * -request.gapPenalty;
        M[i_text] = DIRECTION::LEFT;
    }

    // Dynamic programming loop
    auto thisRowM = M + numCols;
    for (uint64_t i_pattern = 1; i_pattern < numRows; ++i_pattern)
    {
        // Ping-pong score buffers.
        auto tmp = prevRowScores;
        prevRowScores = thisRowScores;
        thisRowScores = tmp;

        thisRowScores[0] = i_pattern * -request.gapPenalty;
        thisRowM[0] = DIRECTION::TOP;

        for (uint64_t i_text = 1; i_text < numCols; ++i_text)
        {
            // Get score for this letter combination. Note that i_text and i_pattern point one
            // beyond the actual text and pattern becuase of the gap character at the beginning.
            const char textByte = request.textBytes[i_text - 1];
            const char patternByte = request.patternBytes[i_pattern - 1];
            const int scoreIdx = ((int) patternByte) * request.alphabetSize + ((int) textByte);

            // Calculate all alignment scores.
            const int fromLeftScore = thisRowScores[i_text - 1] - request.gapPenalty;
            const int fromTopScore = prevRowScores[i_text] - request.gapPenalty;
            const int fromDiagonalScore = prevRowScores[i_text - 1] + request.scoreMatrix[scoreIdx];

            // Find out the best alignment.
            // Order in case of ties: left, top, diagonal.
            const int maxWithGap = std::max(fromLeftScore, fromTopScore);
            const int bestScore = std::max(fromDiagonalScore, maxWithGap);
            const bool isFromDiag = (fromDiagonalScore > maxWithGap);
            const bool isFromLeft = (fromDiagonalScore <= maxWithGap) && (fromLeftScore >= fromTopScore);
            const bool isFromTop = (fromDiagonalScore <= maxWithGap) && (fromLeftScore < fromTopScore);

            // Populate this alignPoint with the best alignment.
            thisRowScores[i_text] = bestScore;
            thisRowM[i_text] = (isFromLeft * DIRECTION::LEFT + isFromDiag * DIRECTION::DIAG + isFromTop * DIRECTION::TOP);
        }

        thisRowM += numCols;
    }

    int score = thisRowScores[numCols -1];

    cleanUp();

    return score;
}


int SequenceAlignment::alignSequenceCPU(const SequenceAlignment::Request &request,
                                        SequenceAlignment::Response *response)
{

    char *M = nullptr;
    /// Aditional row and column for the gap character.
    const uint64_t numCols = request.textNumBytes + 1;
    const uint64_t numRows = request.patternNumBytes + 1;

    auto cleanUp = [&] ()
    {
        if (M) delete [] M;
        M = nullptr;
    };

    /** Allocate memory */
    try
    {
        M = new char[numRows * numCols];
        response->alignedTextBytes = new char[2 * request.textNumBytes];
        response->alignedPatternBytes = new char[2 * request.textNumBytes];

    }
    catch(const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        cleanUp();
        return -1;
    }
    /** End Allocate memory */

    if (request.alignmentType == programArgs::GLOBAL)
    {
        response->score = fillMatrixNW(M, numRows, numCols, request);
        traceBackNW(M, numRows, numCols, request, response);
    }
    else if (request.alignmentType == programArgs::LOCAL)
    {
        auto scoreIdxPair = fillMatrixSW(M, numRows, numCols, request);
        response->score = scoreIdxPair.first;
        traceBackSW(M, scoreIdxPair.second, numRows, numCols, request, response);
    }

    cleanUp();

    return 0;
}
