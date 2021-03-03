#include "SequenceAlignment.hpp"

#include <iostream>
#include <iterator>


using SequenceAlignment::DIR;

void SequenceAlignment::traceBack(const char *M, const uint64_t numRows, const uint64_t numCols,
                                  const Request &request, Response *response)
{
    uint64_t curr = numRows * numCols - 1;
    int textIndex = request.textNumBytes - 1;
    int patternIndex = request.patternNumBytes - 1;

    response->numAlignmentBytes = 0;

    while (curr > 0)
    {
        // Was it a match, gap in text or gap in pattern?
        const bool takeText = (M[curr] == DIR::DIAG || M[curr] == DIR::LEFT);
        const bool takePattern = (M[curr] == DIR::DIAG || M[curr] == DIR::TOP);

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
        curr -= (M[curr] == DIR::LEFT) +
                ((M[curr] == DIR::DIAG) * (numCols+1)) +
                ((M[curr] == DIR::TOP) * (numCols));
    }

    std::reverse(response->alignedTextBytes,
                 (response->alignedTextBytes + response->numAlignmentBytes));
    std::reverse(response->alignedPatternBytes,
                 (response->alignedPatternBytes + response->numAlignmentBytes));
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
        M[i_text] = DIR::LEFT;
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
        thisRowM[0] = DIR::TOP;

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
            thisRowM[i_text] = (isFromLeft * DIR::LEFT + isFromDiag * DIR::DIAG + isFromTop * DIR::TOP);
        }

        thisRowM += numCols;
    }

    auto score = thisRowScores[numCols -1];

    cleanUp();

    return score;
}


int SequenceAlignment::alignSequenceGlobalCPU(const SequenceAlignment::Request &request,
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

    response->score = fillMatrixNW(M, numRows, numCols, request);

    traceBack(M, numRows, numCols, request, response);

    cleanUp();

    return 0;
}
