#include "SequenceAlignment.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>


/// Given 2 8-bit characters, construct a unique 16-bit key.
short twoCharsToKey(char char1, char char2)
{
    // Cast to short, shift char1 to upper 8 bits and OR char2 to lower 8 bits.
    return (((short) char1) << 8) | char2;
}

void parseScoreMatrixFile(const std::string fname)
{
    std::ifstream f(fname);
    if (f.good())
    {
        const auto sequenceOrderedChars = (SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA)
                                          ? SequenceAlignment::dnaScoreMatrixCharOrder
                                          : SequenceAlignment::proteinScoreMatrixCharOrder;
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
                {
                    std::cerr << SequenceAlignment::SCORE_MATRIX_NOT_READ_WARNING;
                    SequenceAlignment::scoreMap.clear();
                    return;
                }

                const short key = twoCharsToKey(sequenceOrderedChars[i], sequenceOrderedChars[j]);
                SequenceAlignment::scoreMap[key] = nextScore;
            }
        }
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
    }
}

/// Given a filename, copy the chars (bytes) from the file into a buffer.
/// If the text sequence buffer is empty, it is filled.
/// Else if the pattern sequence buffer is empty, it is filled.
/// Else the file is ignored.
void readSequenceBytes(const std::string fname)
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
    }
}

void parseArguments(int argc, const char *argv[])
{
    // We need at least the text and pattern file.
    if (argc < 3)
    {
        std::cerr << SequenceAlignment::USAGE;
        return;
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
            parseScoreMatrixFile(argv[i]);
            nextIsScoreMatrixFile = false;

        }
        else
        {
            readSequenceBytes(argv[i]);
        }
    }

    if (SequenceAlignment::textNumBytes == 0 || SequenceAlignment::patternNumBytes == 0)
    {
        std::cerr << SequenceAlignment::SEQ_NOT_READ_ERROR << SequenceAlignment::USAGE;
        return;
    }

    // Read in defualt scores.
    if (SequenceAlignment::scoreMap.size() == 0)
    {
        const auto scoreMatrixFile = (SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA)
                                     ? SequenceAlignment::defaultDnaScoreMatrixFile
                                     : SequenceAlignment::defaultProteinScoreMatrixFile;
        parseScoreMatrixFile(scoreMatrixFile);
    }
}
