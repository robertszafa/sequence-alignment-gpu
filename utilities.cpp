#include "SequenceAlignment.hpp"

#include <iostream>
#include <iomanip>
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
             const int *scoreMatrix)
{
    int idx = indexOfLetter(char1, alphabet, alphabetSize) * alphabetSize +
              indexOfLetter(char2, alphabet, alphabetSize);
    return scoreMatrix[idx];
}

/// Given a string of letters:
///     - remove character, if not in the alphabet
///     - transorm character into i where character=alphabet[i], if not in the alphabet
int validateAndTransform(std::string &sequence, const char *alphabet, const int alphabetSize)
{
    unsigned int numRead = 0;
    for (int i=0; i<sequence.length(); ++i)
    {
        const char upperLetter = (sequence[i] > 90) ? sequence[i] - 32 : sequence[i];
        // Ignore characters not on the A-Z range.
        if (upperLetter < 65 || upperLetter > 90) continue;

        sequence[numRead] = indexOfLetter(upperLetter, alphabet, alphabetSize);
        if (sequence[numRead] == -1)
        {
            std::cerr << "'" << sequence[i] << "'" << " letter not in alphabet." << std::endl;
            return 0;
        }

        ++numRead;
    }

    return numRead;
}

int readSequenceFile(const std::string fname, SequenceAlignment::Request *request)
{
    std::ifstream f(fname);
    std::string fileContents = "";
    if (f.good())
        // Use string's range constructor to copy over entire file to memory.
        fileContents = std::string((std::istreambuf_iterator<char>(f)),
                                    std::istreambuf_iterator<char>());
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
        return -1;
    }

    int numLetters = validateAndTransform(fileContents, request->alphabet, request->alphabetSize);
    try
    {
        if (request->textNumBytes == 0 && numLetters > 0)
        {
            // Allocate memory for buffer and copy the normalised sequence.
            request->textBytes = new char[numLetters];
            request->textNumBytes = numLetters;
            std::copy_n(fileContents.begin(), numLetters, request->textBytes);
        }
        else if (request->patternNumBytes == 0 && numLetters > 0)
        {
            // Allocate memory for buffer and copy the normalised sequence.
            request->patternBytes = new char[numLetters];
            request->patternNumBytes = numLetters;
            std::copy_n(fileContents.begin(), numLetters, request->patternBytes);
        }
    }
    catch (const std::bad_alloc& e)
    {
        std::cerr << SequenceAlignment::MEM_ERROR;
        return -1;
    }

    return 0;
}

int parseScoreMatrixFile(const std::string& fname, const int alphabetSize, int *buffer)
{
    std::ifstream f(fname);
    if (f.good())
    {
        int nextScore;
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
    request->gapPenalty = SequenceAlignment::DEFAULT_GAP_PENALTY;
    request->textNumBytes = 0;
    request->patternNumBytes = 0;

    enum flagState {NOT_READ, TO_READ, READ};
    flagState scoreMatrixState = flagState::NOT_READ;
    flagState gapPenaltyState = flagState::NOT_READ;
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
            gapPenaltyState = (setArg == SequenceAlignment::programArgs::GAP_PENALTY)
                              ? flagState::TO_READ
                              : gapPenaltyState;
        }
        else if (gapPenaltyState == flagState::TO_READ)
        {
            try
            {
                request->gapPenalty = std::stoi(argv[i]);
            }
            catch (...) // std::invalid_argument, std::out_of_range
            {
                std::cerr << SequenceAlignment::GAP_PENALTY_NOT_READ_ERROR;
                return -1;
            }
            gapPenaltyState = flagState::READ;
        }
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
            if (readSequenceFile(argv[i], request) == -1)
            {
                std::cerr << SequenceAlignment::SEQ_NOT_READ_ERROR;
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


void prettyAlignmentPrint(SequenceAlignment::Response &response, std::ostream &stream)
{
    if (response.numAlignmentBytes == 0) return;

    const unsigned int CHARS_PER_LINE = 50;

    int charNumWidth = 0;
    int maxI = response.numAlignmentBytes;
    do { maxI /= 10; ++charNumWidth; } while (maxI != 0);

    int numIdentity = 0, numGaps = 0;
    for (int i = 0; i < response.numAlignmentBytes; i += CHARS_PER_LINE)
    {
        stream << std::setfill(' ') << std::setw(charNumWidth) << (i + 1) << " "
               << std::resetiosflags(std::ios::showbase);

        int j = i;
        for (j = i;j < (i+CHARS_PER_LINE) && j < response.numAlignmentBytes; j++)
        {
            stream << response.alignedTextBytes[j];
        }
        stream << "   " << j
               << " \n" << std::setfill(' ') << std::setw(charNumWidth) << " " << " "
               << std::resetiosflags(std::ios::showbase);


        for (j = i; j < (i+CHARS_PER_LINE) && j < response.numAlignmentBytes; j++)
        {
            if (response.alignedTextBytes[j] == response.alignedPatternBytes[j])
            {
                stream << '|';
                ++numIdentity;
            }
            else if (response.alignedTextBytes[j] == '-' || response.alignedPatternBytes[j] == '-')
            {
                stream << ' ';
                ++numGaps;
            }
            else
                stream << '.';
        }
        stream << "\n" << std::setfill(' ') << std::setw(charNumWidth) << (i + 1) << " "
               << std::resetiosflags(std::ios::showbase);

        for (j = i; j < (i+CHARS_PER_LINE) && j < response.numAlignmentBytes; j++)
        {
            stream << response.alignedPatternBytes[j];
        }

        stream << "   " << j << "\n\n";
    }

    stream << "# Length: \t" << response.numAlignmentBytes << "\n"
           << "# Identity: \t" <<  numIdentity << "/" << response.numAlignmentBytes << std::setprecision(3)
                             << " (" << (numIdentity/(response.numAlignmentBytes*1.0)*100) << "%)\n"
           << "# Gaps: \t" <<  numGaps << "/" << response.numAlignmentBytes << std::setprecision(3)
                             << " (" << (numGaps/(response.numAlignmentBytes*1.0)*100) << "%)\n"
           << "# Score: \t" << response.score << "\n";
}
