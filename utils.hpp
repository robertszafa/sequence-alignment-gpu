#include "SequenceAlignment.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>


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
                f >> nextScore;
                short char1 = (short) sequenceOrderedChars[i];
                short char2 = (short) sequenceOrderedChars[j];

                const short key = (char1 << 8) & char2;
                SequenceAlignment::scoreMap[key] = nextScore;

                std::cout << sequenceOrderedChars[i] << ", " << sequenceOrderedChars[j] << " = " << nextScore <<"\n";
            }
        }
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
    }

    f.close();

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
        // Device and sequence type flags.
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
        else if (nextIsScoreMatrixFile)
        {
            parseScoreMatrixFile(argv[i]);
            nextIsScoreMatrixFile = false;
        }
        else // Sequence file names.
        {
            std::ifstream f(argv[i]);
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
                else
                {
                    SequenceAlignment::patternNumBytes = fileString.length();
                    std::copy_n(fileString.begin(),
                                SequenceAlignment::patternNumBytes,
                                SequenceAlignment::patternBytes);
                }

            }
            else
            {
                std::cerr << argv[i] << " file does not exist" << std::endl;
            }

            f.close();
        }
    }

    if (SequenceAlignment::textNumBytes == 0 || SequenceAlignment::patternNumBytes == 0)
    {
        std::cerr << "textSequence or patternSequence not read" << std::endl;
        std::cerr << SequenceAlignment::USAGE;
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
