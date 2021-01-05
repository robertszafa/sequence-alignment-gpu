#include "SequenceAlignment.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>

int getScoreIndex(char char1, char char2)
{
    const auto predicate = (SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA);
    const auto sequenceOrderedChars = predicate
                                      ? SequenceAlignment::dnaScoreMatrixCharOrder
                                      : SequenceAlignment::proteinScoreMatrixCharOrder;
    const auto numChars = predicate
                          ? SequenceAlignment::NUM_DNA_CHARS
                          : SequenceAlignment::NUM_PROTEIN_CHARS;

    int val1 = std::find(sequenceOrderedChars, sequenceOrderedChars+numChars, char1) - sequenceOrderedChars;
    int val2 = std::find(sequenceOrderedChars, sequenceOrderedChars+numChars, char2) - sequenceOrderedChars;

    return val1 * numChars + val2;
}

int parseScoreMatrixFile(const std::string fname)
{
    std::ifstream f(fname);
    if (f.good())
    {
        const auto numChars = (SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA)
                            ? SequenceAlignment::NUM_DNA_CHARS
                            : SequenceAlignment::NUM_PROTEIN_CHARS;

        short nextScore;
        for (int i=0; i<numChars; ++i)
        {
            for (int j=0; j<numChars; ++j)
            {
                // Check if valid.
                if (!(f >> nextScore))
                    return -1;

                SequenceAlignment::scoreMatrix[i * numChars + j] = nextScore;
            }
        }
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
    }

    return 0;
}

/// Given a filename, copy the chars (bytes) from the file into a buffer.
/// If the text sequence buffer is empty, it is filled.
/// Else if the pattern sequence buffer is empty, it is filled.
/// Else the file is ignored.
int readSequenceBytes(const std::string fname)
{
    std::ifstream f(fname);
    if (f.good())
    {
        // Use string's range constructor to copy over entire file to memory.
        std::string fileString((std::istreambuf_iterator<char>(f)),
                                std::istreambuf_iterator<char>());

        if (SequenceAlignment::textNumBytes == 0)
        {
            SequenceAlignment::textNumBytes = fileString.length();
            std::copy_n(fileString.begin(),
                        SequenceAlignment::textNumBytes,
                        SequenceAlignment::textBytes);
        }
        else if (SequenceAlignment::patternNumBytes == 0)
        {
            SequenceAlignment::patternNumBytes = fileString.length();
            std::copy_n(fileString.begin(),
                        SequenceAlignment::patternNumBytes,
                        SequenceAlignment::patternBytes);
        }
        else
        {
            std::cerr << "Ignoring file " << fname << std::endl;
        }
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
        return -1;
    }

    return 0;
}

int parseArguments(int argc, const char *argv[])
{
    // We need at least the text and pattern file.
    if (argc < 3)
    {
        std::cerr << SequenceAlignment::USAGE;
        return -1;
    }

    bool nextIsScoreMatrixFile = false;
    for (int i = 1; i < argc; ++i)
    {
        // Flags
        if (argv[i][0] == '-' && strlen(argv[i]) > 1)
        {
            // Check if the flag is valid.
            if (SequenceAlignment::argumentMap.count(argv[i][1]) > 0)
            {
                auto setArg = SequenceAlignment::argumentMap.at(argv[i][1]);
                SequenceAlignment::deviceType = (setArg == SequenceAlignment::programArgs::CPU) ||
                                                (setArg == SequenceAlignment::programArgs::GPU)
                                                ? setArg
                                                : SequenceAlignment::deviceType;
                SequenceAlignment::sequenceType = (setArg == SequenceAlignment::programArgs::DNA) ||
                                                  (setArg == SequenceAlignment::programArgs::PROTEIN)
                                                  ? setArg
                                                  : SequenceAlignment::sequenceType;
                nextIsScoreMatrixFile = (setArg == SequenceAlignment::programArgs::SCORE_MATRIX);
            }
            else
            {
                std::cerr << "Ignoring \"" << argv[i] << "\"" << std::endl;
            }
        }
        // Files to read
        else if (nextIsScoreMatrixFile)
        {
            if (parseScoreMatrixFile(argv[i]) == -1)
            {
                // Read in defualt scores.
                std::cerr << SequenceAlignment::SCORE_MATRIX_NOT_READ_WARNING;
                const auto scoreMatrixFile = (SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA)
                                            ? SequenceAlignment::defaultDnaScoreMatrixFile
                                            : SequenceAlignment::defaultProteinScoreMatrixFile;
                parseScoreMatrixFile(scoreMatrixFile);
            }
            nextIsScoreMatrixFile = false;
        }
        else
        {
            if (readSequenceBytes(argv[i]) == -1)
            {
                std::cerr << SequenceAlignment::SEQ_NOT_READ_ERROR << SequenceAlignment::USAGE;
                return -1;
            }
        }
    }

    if (SequenceAlignment::textNumBytes == 0 || SequenceAlignment::patternNumBytes == 0)
    {
        std::cerr << SequenceAlignment::SEQ_NOT_READ_ERROR << SequenceAlignment::USAGE;
        return -1;
    }


    return 0;
}
