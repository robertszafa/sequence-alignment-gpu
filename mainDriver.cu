#include "SequenceAlignment.hpp"
#include "utilities.hpp"

#include <cuda.h>

#include <iostream>
#include <fstream>


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
    request->alphabet = SequenceAlignment::DEFAULT_ALPHABET;
    request->alphabetSize = SequenceAlignment::DEFAULT_ALPHABET_SIZE;
    request->gapOpen = SequenceAlignment::DEFAULT_GAP_OPEN;
    request->gapExtend = SequenceAlignment::DEFAULT_GAP_EXTEND;
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
                request->gapOpen = std::stoi(argv[i]);
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
                request->gapExtend = std::stoi(argv[i]);
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

    if (scoreMatrixState != flagState::READ)
    {
        auto defaultScores = (request->sequenceType == SequenceAlignment::programArgs::DNA)
                             ? SequenceAlignment::DEFAULT_DNA_SCORE_MATRIX_FILE
                             : SequenceAlignment::DEFAULT_PROTEIN_SCORE_MATRIX_FILE;
        parseScoreMatrixFile(defaultScores, request->alphabetSize, request->scoreMatrix);
    }

    return 0;
}


int main(int argc, const char *argv[])
{
    SequenceAlignment::Request request;
    SequenceAlignment::Response response;

    // Fill the request with user specified arguments.
    if (parseArguments(argc, argv, &request) == -1) return -1;

    if (request.deviceType == SequenceAlignment::programArgs::CPU)
        SequenceAlignment::alignSequenceCPU(request, &response);

    return 0;
}
