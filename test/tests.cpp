#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include "../SequenceAlignment.hpp"
#include "../utilities.hpp"


TEST_CASE("indexOfLetter")
{
    REQUIRE(indedOfLetter('A', SequenceAlignment::DNA_ALPHABET, SequenceAlignment::NUM_DNA_CHARS) == 0);
    REQUIRE(indedOfLetter('H', SequenceAlignment::DNA_ALPHABET, SequenceAlignment::NUM_DNA_CHARS) == -1);
    REQUIRE(indedOfLetter('H', SequenceAlignment::PROTEIN_ALPHABET, SequenceAlignment::NUM_PROTEIN_CHARS) == 8);
}

TEST_CASE("parseScoreMatrixFile")
{
    SequenceAlignment::Request request;
    request.alphabet = SequenceAlignment::DNA_ALPHABET;
    request.alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
    parseScoreMatrixFile("scoreMatrices/dna/blast.txt", request.alphabetSize, request.scoreMatrix);

    REQUIRE(getScore('A', 'A', request.alphabet, request.alphabetSize, request.scoreMatrix) == 5);
    REQUIRE(getScore('G', 'T', request.alphabet, request.alphabetSize, request.scoreMatrix) == -4);
}

TEST_CASE("readSequenceBytes")
{
    const int argc = 3;
    const char *argv[argc] = { "./alignSequence", "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
    SequenceAlignment::Request request;
    parseArguments(argc, argv, &request);

    const char expectedText[] = {0, 2, 0, 2, 3, 2, 1, 0, 3};
    const char expectedPattern[] = {2, 2, 1, 0, 1, 3, 3, 2, 1, 3};

    REQUIRE(std::equal(request.textBytes, request.textBytes + request.textNumBytes, expectedText));
    REQUIRE(std::equal(request.patternBytes, request.patternBytes + request.patternNumBytes, expectedPattern));
}

TEST_CASE("parseArguments")
{
    // Catch stderr to string to test error messages.
    std::stringstream buffer;
    std::streambuf *old = std::cerr.rdbuf(buffer.rdbuf());

    SECTION("usage")
    {
        const int argc = 1;
        const char *argv[argc] = { "./alignSequence"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string stderrString = buffer.str();
        REQUIRE(stderrString == SequenceAlignment::USAGE);
    }

    SECTION("no or empty sequence files")
    {
        const int argc = 3;
        const char *argv[argc] = { "./alignSequence", "-p", "-c"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string expectedMsg = SequenceAlignment::SEQ_NOT_READ_ERROR + SequenceAlignment::USAGE;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);

        REQUIRE(request.deviceType == SequenceAlignment::programArgs::CPU);
        REQUIRE(request.sequenceType == SequenceAlignment::programArgs::PROTEIN);
    }

    SECTION("incorrect score matrix")
    {
        const int argc = 5;
        const char *argv[argc] = { "./alignSequence", "--score-matrix", "test/corruptScoreMatrix.txt",
                                   "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        SequenceAlignment::Request request;
        parseArguments(argc, argv, &request);

        std::string expectedMsg = SequenceAlignment::SCORE_MATRIX_NOT_READ_ERROR;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);
    }

    // Restore old cerr.
    std::cerr.rdbuf(old);
}

TEST_CASE("alignSequenceCPU")
{
    // const int argc = 3;
    // const char *argv[argc] = { "./alignSequence", "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
    // parseArguments(argc, argv);

    // SequenceAlignment::alignSequenceCPU();

    // const std::string expectedAlignment = "CCGCTG";
    // auto gotAlignment = std::string(request.textBytes,
    //                                 request.textBytes + request.alignmentNumBytes);

    // REQUIRE(expectedAlignment == gotAlignment);
}