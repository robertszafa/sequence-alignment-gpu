#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"

#include "../SequenceAlignment.hpp"
#include "../utils.hpp"

#include <algorithm>


TEST_CASE( "parseArguments")
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

        std::string expectedMsg = SequenceAlignment::SCORE_MATRIX_NOT_READ_WARNING;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);
    }

    SECTION("reading sequences from files")
    {
        const int argc = 5;
        const char *argv[argc] = { "./alignSequence", "-g", "-d",
                                   "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        parseArguments(argc, argv);

        const char* expectedText = "ACACGCTAG";
        const char* expectedPattern = "CCTATGGCTG";

        REQUIRE(SequenceAlignment::deviceType == SequenceAlignment::programArgs::GPU);
        REQUIRE(SequenceAlignment::sequenceType == SequenceAlignment::programArgs::DNA);

        REQUIRE(std::equal(SequenceAlignment::textBytes,
                           SequenceAlignment::textBytes + SequenceAlignment::textNumBytes,
                           expectedText));
        REQUIRE(std::equal(SequenceAlignment::patternBytes,
                           SequenceAlignment::patternBytes + SequenceAlignment::patternNumBytes,
                           expectedPattern));
    }

    SECTION("reading score matrix from file")
    {
        const int argc = 6;
        const char *argv[argc] = { "./alignSequence", "-d", "-m", "scoreMatrices/dna/blast.txt",
                                   "data/dna/dna_01.txt", "data/dna/dna_02.txt"};
        parseArguments(argc, argv);

        REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('A', 'A')] == 5);
        REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('G', 'T')] == -4);
        REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('C', '*')] == -1);
        REQUIRE(SequenceAlignment::scoreMatrix[getScoreIndex('*', '*')] == 1);
    }
}
