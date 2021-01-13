#pragma once

#include <iostream>
#include <algorithm>
#include <fstream>

/// Given a letter and an alphabet array,
/// return the index of that letter in the  alphabet array.
char indedOfLetter(const char letter, const char *alphabet, const int alphabetSize)
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
    int idx = indedOfLetter(char1, alphabet, alphabetSize) * alphabetSize +
              indedOfLetter(char2, alphabet, alphabetSize);
    return scoreMatrix[idx];
}

int validateAndTransform(const std::string &sequence, const char *alphabet, const int alphabetSize,
                         char *dstBuffer)
{
    for (int i=0; i<sequence.length(); ++i)
    {
        const char upperLetter = (sequence[i] > 90) ? sequence[i] - 32 : sequence[i];

        dstBuffer[i] = indedOfLetter(upperLetter, alphabet, alphabetSize);
        if (dstBuffer[i] == -1)
        {
            std::cerr << "'" << sequence[i] << "'" << " letter not in alphabet." << std::endl;
            return -1;
        }
    }

    return 0;
}

/// Given a filename, copy the chars (bytes) from the file into a buffer.
/// Return the number of bytes read or -1 if an error occured.
int readSequenceBytes(const std::string fname,  const char *alphabet, const int alphabetSize,
                      char *dstBuffer)
{
    std::ifstream f(fname);
    int numBytesRead = 0;
    if (f.good())
    {
        // Use string's range constructor to copy over entire file to memory.
        std::string fileString((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
        if (validateAndTransform(fileString, alphabet, alphabetSize, dstBuffer) == -1) return -1;
        numBytesRead = fileString.length();
    }
    else
    {
        std::cerr << fname << " file does not exist" << std::endl;
        return -1;
    }

    return numBytesRead;
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
