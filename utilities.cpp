#include "SequenceAlignment.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>

/// Given a letter and an alphabet array,
/// return the index of that letter in the  alphabet array.
char indexOfLetter(const char letter, const char *alphabet, const int alphabetSize)
{
    auto letterItr = std::find(alphabet, alphabet + alphabetSize, letter);
    if (letterItr == (alphabet + alphabetSize)) return -1;
    return std::distance(alphabet, letterItr);
}

/// Given two characters, an alphabet and a scoring matrix,
/// return the score of the character combination.
int getScore(char char1, char char2, const char *alphabet, const int alphabetSize,
             const short *scoreMatrix)
{
    int idx = indexOfLetter(char1, alphabet, alphabetSize) * alphabetSize +
              indexOfLetter(char2, alphabet, alphabetSize);
    return scoreMatrix[idx];
}

int validateAndTransform(const std::string &sequence, const char *alphabet, const int alphabetSize,
                         char *dstBuffer)
{
    unsigned int numRead = 0;
    for (int i=0; i<sequence.length(); ++i)
    {
        const char upperLetter = (sequence[i] > 90) ? sequence[i] - 32 : sequence[i];
        // Ignore characters not on the A-Z range.
        if (upperLetter < 65 || upperLetter > 90) continue;

        dstBuffer[numRead] = indexOfLetter(upperLetter, alphabet, alphabetSize);
        if (dstBuffer[numRead] == -1)
        {
            std::cerr << "'" << sequence[i] << "'" << " letter not in alphabet." << std::endl;
            return -1;
        }

        ++numRead;
    }

    return numRead;
}

/// Given a filename, copy the chars (bytes) from the file into a buffer.
/// Return the number of bytes read or -1 if an error occured.
int readSequenceBytes(const std::string fname,  const char *alphabet, const int alphabetSize,
                      char *dstBuffer)
{
    std::ifstream f(fname);
    if (f.good())
    {
        // Use string's range constructor to copy over entire file to memory.
        std::string fileString((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        return validateAndTransform(fileString, alphabet, alphabetSize, dstBuffer);
    }

    std::cerr << fname << " file does not exist" << std::endl;
    return -1;
}

int parseScoreMatrixFile(const std::string& fname, const int alphabetSize, short *buffer)
{
    std::ifstream f(fname);
    if (f.good())
    {
        short nextScore;
        for (int i = 0; i < alphabetSize; ++i)
        {
            for (int j = 0; j < alphabetSize; ++j)
            {
                // Check if 16-bit number.
                if (!(f >> nextScore)) return -1;

                buffer[i*alphabetSize + j] = nextScore;
            }
        }
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
    }

    return 0;
}

int parseArguments(int argc, const char *argv[], SequenceAlignment::Request *request)
{
    // We need at least the text and pattern file.
    if (argc == 1)
    {
        std::cerr << SequenceAlignment::USAGE;
        return -1;
    }

    // Set all to default.
    request->deviceType = SequenceAlignment::DEFAULT_DEVICE;
    request->sequenceType = SequenceAlignment::DEFAULT_SEQUENCE;
    request->alignmentType = SequenceAlignment::DEFAULT_ALIGNMENT_TYPE;
    request->alphabet = SequenceAlignment::DEFAULT_ALPHABET;
    request->alphabetSize = SequenceAlignment::DEFAULT_ALPHABET_SIZE;
    request->gapOpenScore = SequenceAlignment::DEFAULT_GAP_OPEN_SCORE;
    request->gapExtendScore = SequenceAlignment::DEFAULT_GAP_EXTEND_SCORE;
    request->textNumBytes = 0;
    request->patternNumBytes = 0;

    enum flagState {NOT_READ, TO_READ, READ};
    flagState scoreMatrixState = flagState::NOT_READ;
    flagState gapOpenState = flagState::NOT_READ;
    flagState gapExtendState = flagState::NOT_READ;
    for (int i = 1; i < argc; ++i)
    {
        // Check if it's a flag.
        if (SequenceAlignment::argumentMap.count(argv[i]) > 0)
        {
            auto setArg = SequenceAlignment::argumentMap.at(argv[i]);

            request->deviceType = (setArg == SequenceAlignment::programArgs::CPU ||
                                   setArg == SequenceAlignment::programArgs::GPU)
                                   ? setArg
                                   : request->deviceType;
            request->sequenceType = (setArg == SequenceAlignment::programArgs::DNA ||
                                     setArg == SequenceAlignment::programArgs::PROTEIN)
                                     ? setArg
                                     : request->sequenceType;
            request->alignmentType = (setArg == SequenceAlignment::programArgs::GLOBAL ||
                                     setArg == SequenceAlignment::programArgs::LOCAL  ||
                                     setArg == SequenceAlignment::programArgs::SEMI_GLOBAL)
                                     ? setArg
                                     : request->alignmentType;

            // Update alphabet pointers.
            bool isDna = (request->sequenceType == SequenceAlignment::programArgs::DNA);
            request->alphabet = isDna ? SequenceAlignment::DNA_ALPHABET
                                      : SequenceAlignment::PROTEIN_ALPHABET;
            request->alphabetSize = isDna ? SequenceAlignment::NUM_DNA_CHARS
                                          : SequenceAlignment::NUM_PROTEIN_CHARS;

            scoreMatrixState = (setArg == SequenceAlignment::programArgs::SCORE_MATRIX)
                               ? flagState::TO_READ
                               : scoreMatrixState;
            gapOpenState = (setArg == SequenceAlignment::programArgs::GAP_OPEN)
                           ? flagState::TO_READ
                           : gapOpenState;
            gapExtendState = (setArg == SequenceAlignment::programArgs::GAP_EXTEND)
                             ? flagState::TO_READ
                             : gapExtendState;
        }
        else if (gapOpenState == flagState::TO_READ)
        {
            try
            {
                request->gapOpenScore = -std::stoi(argv[i]);
            }
            catch (...) // std::invalid_argument, std::out_of_range
            {
                std::cerr << SequenceAlignment::GAP_OPEN_NOT_READ_ERROR;
                return -1;
            }
            gapOpenState = flagState::READ;
        }
        else if (gapExtendState == flagState::TO_READ)
        {
            try
            {
                request->gapExtendScore = -std::stoi(argv[i]);
            }
            catch (...) // std::invalid_argument, std::out_of_range
            {
                std::cerr << SequenceAlignment::GAP_EXTEND_NOT_READ_ERROR;
                return -1;
            }
            gapExtendState = flagState::READ;
        }
        // Files to read
        else if (scoreMatrixState == flagState::TO_READ)
        {
            if (parseScoreMatrixFile(argv[i], request->alphabetSize, request->scoreMatrix) == -1)
            {
                std::cerr << SequenceAlignment::SCORE_MATRIX_NOT_READ_ERROR;
                return -1;
            }
            scoreMatrixState = flagState::READ;
        }
        else
        {
            if (request->textNumBytes == 0)
                request->textNumBytes = readSequenceBytes(argv[i], request->alphabet,
                                                          request->alphabetSize, request->textBytes);
            else if (request->patternNumBytes == 0)
                request->patternNumBytes = readSequenceBytes(argv[i], request->alphabet,
                                                          request->alphabetSize, request->patternBytes);

            // Check if read correctly.
            if (request->textNumBytes == -1 || request->patternNumBytes == -1)
            {
                std::cerr << SequenceAlignment::SEQ_NOT_READ_ERROR << SequenceAlignment::USAGE;
                return -1;
            }
        }
    }

    if (request->textNumBytes == 0 || request->patternNumBytes == 0)
    {
        std::cerr << SequenceAlignment::SEQ_NOT_READ_ERROR << SequenceAlignment::USAGE;
        return -1;
    }
    else if (request->textNumBytes < request->patternNumBytes)
    {
        std::cerr << SequenceAlignment::TEXT_SHORTER_THAN_PATTERN_ERROR;
        return -1;
    }


    if (scoreMatrixState != flagState::READ)
    {
        auto defaultScores = (request->sequenceType == SequenceAlignment::programArgs::DNA)
                             ? SequenceAlignment::DEFAULT_DNA_SCORE_MATRIX_FILE
                             : SequenceAlignment::DEFAULT_PROTEIN_SCORE_MATRIX_FILE;
        parseScoreMatrixFile(defaultScores, request->alphabetSize, request->scoreMatrix);
    }

    return 0;
}

