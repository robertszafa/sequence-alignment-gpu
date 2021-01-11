#include "SequenceAlignment.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>

/// Given a letter from the dna or protein alphabet return the index of that letter
/// in the ordered alphabet array.
char letterToInt(char letter)
{
    auto letterItr = std::find(SequenceAlignment::alphabet,
                               SequenceAlignment::alphabet + SequenceAlignment::alphabetSize,
                               letter);
    if (letterItr == (SequenceAlignment::alphabet + SequenceAlignment::alphabetSize)) return -1;

    return std::distance(SequenceAlignment::alphabet, letterItr);
}

int getScoreIndex(char char1, char char2)
{
    return letterToInt(char1) * SequenceAlignment::alphabetSize + letterToInt(char2);
}

int validateAndTransform(const std::string &sequence, char *dstBuffer)
{
    for (int i=0; i<sequence.length(); ++i)
    {
        const char upperLetter = (sequence[i] > 90) ? sequence[i] - 32 : sequence[i];

        dstBuffer[i] = letterToInt(upperLetter);
        if (dstBuffer[i] == -1)
        {
            std::cerr << "'" << sequence[i] << "'" << " letter not in alphabet." << std::endl;
            return -1;
        }
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
        std::string fileString((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

        if (SequenceAlignment::textNumBytes == 0)
        {
            SequenceAlignment::textNumBytes = fileString.length();
            if (validateAndTransform(fileString, SequenceAlignment::textBytes) == -1) return -1;
        }
        else if (SequenceAlignment::patternNumBytes == 0)
        {
            SequenceAlignment::patternNumBytes = fileString.length();
            if (validateAndTransform(fileString, SequenceAlignment::patternBytes) == -1) return -1;
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

int parseScoreMatrixFile(const std::string fname)
{
    std::ifstream f(fname);
    if (f.good())
    {
        short nextScore;
        for (int i=0; i<SequenceAlignment::alphabetSize; ++i)
        {
            for (int j=0; j<SequenceAlignment::alphabetSize; ++j)
            {
                // Check if 16-bit number.
                if (!(f >> nextScore)) return -1;

                SequenceAlignment::scoreMatrix[i * SequenceAlignment::alphabetSize + j] = nextScore;
            }
        }
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
    }

    return 0;
}

int parseArguments(int argc, const char *argv[])
{
    // We need at least the text and pattern file.
    if (argc == 1)
    {
        std::cerr << SequenceAlignment::USAGE;
        return -1;
    }

    // Reset all program parameters to default.
    SequenceAlignment::deviceType = SequenceAlignment::programArgs::CPU;
    SequenceAlignment::sequenceType = SequenceAlignment::programArgs::DNA;
    SequenceAlignment::alphabet = SequenceAlignment::DNA_ALPHABET;
    SequenceAlignment::alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
    SequenceAlignment::textNumBytes = 0;
    SequenceAlignment::patternNumBytes = 0;

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

                bool isDna = SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA;
                SequenceAlignment::alphabet = isDna
                                              ? SequenceAlignment::DNA_ALPHABET
                                              : SequenceAlignment::PROTEIN_ALPHABET;
                SequenceAlignment::alphabetSize = isDna
                                              ? SequenceAlignment::NUM_DNA_CHARS
                                              : SequenceAlignment::NUM_PROTEIN_CHARS;

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
                std::cerr << SequenceAlignment::SCORE_MATRIX_NOT_READ_ERROR;
                return -1;
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
