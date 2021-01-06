#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include "../SequenceAlignment.hpp"
#include "../utils.hpp"

#include <algorithm>


TEST_CASE( "parseArguments - error messages")
{
    // Catch stderr to string to test error messages.
    std::stringstream buffer;
    std::streambuf *old = std::cerr.rdbuf(buffer.rdbuf());

    SECTION("usage")
    {
        const int argc = 1;
        const char *argv[argc] = { "./alignSequence"};
        parseArguments(argc, argv);

        std::string stderrString = buffer.str();
        REQUIRE(stderrString == SequenceAlignment::USAGE);
    }

    SECTION("no or empty sequence files")
    {
        const int argc = 3;
        const char *argv[argc] = { "./alignSequence", "-p", "-c"};
        parseArguments(argc, argv);

        std::string expectedMsg = SequenceAlignment::SEQ_NOT_READ_ERROR + SequenceAlignment::USAGE;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);

        REQUIRE(SequenceAlignment::deviceType == SequenceAlignment::programArgs::CPU);
        REQUIRE(SequenceAlignment::sequenceType == SequenceAlignment::programArgs::PROTEIN);
    }

    SECTION("incorrect score matrix")
    {
        const int argc = 5;
        const char *argv[argc] = { "./alignSequence", "-m", "test/corruptScoreMatrix.txt",
                                   "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        parseArguments(argc, argv);

        std::string expectedMsg = SequenceAlignment::SCORE_MATRIX_NOT_READ_ERROR;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);
    }

    // Restore old cerr.
    std::cerr.rdbuf(old);

}

TEST_CASE("validateAndTransform")
{
    SequenceAlignment::textNumBytes = 9;
    SequenceAlignment::patternNumBytes = 10;
    SequenceAlignment::sequenceType = SequenceAlignment::programArgs::DNA;
    SequenceAlignment::alphabet = SequenceAlignment::DNA_ALPHABET;
    SequenceAlignment::alphabetSize = SequenceAlignment::NUM_DNA_CHARS;

    std::ifstream textFile("data/dna/dna_01.txt");
    std::ifstream patternFile("data/dna/dna_02.txt");
    std::string textString((std::istreambuf_iterator<char>(textFile)), std::istreambuf_iterator<char>());
    std::string patternString((std::istreambuf_iterator<char>(patternFile)), std::istreambuf_iterator<char>());

    validateAndTransform(textString, SequenceAlignment::textBytes);
    validateAndTransform(patternString, SequenceAlignment::patternBytes);

    const char expectedText[] = {0, 2, 0, 2, 3, 2, 1, 0, 3};
    const char expectedPattern[] = {2, 2, 1, 0, 1, 3, 3, 2, 1, 3};

    REQUIRE(std::equal(SequenceAlignment::textBytes,
                       SequenceAlignment::textBytes + SequenceAlignment::textNumBytes,
                       expectedText));
    REQUIRE(std::equal(SequenceAlignment::patternBytes,
                        SequenceAlignment::patternBytes + SequenceAlignment::patternNumBytes,
                        expectedPattern));
}

TEST_CASE("parseScoreMatrixFile")
{
    SequenceAlignment::alphabet = SequenceAlignment::DNA_ALPHABET;
    SequenceAlignment::alphabetSize = SequenceAlignment::NUM_DNA_CHARS;
    parseScoreMatrixFile("scoreMatrices/dna/blast.txt");

    REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('A', 'A')] == 5);
    REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('G', 'T')] == -4);
    REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('C', '*')] == -1);
    REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('*', '*')] == 1);
}

