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

        std::string expectedMsg = "textSequence or patternSequence not read\n" + SequenceAlignment::USAGE;
        std::string stderrString = buffer.str();
        REQUIRE(stderrString == expectedMsg);

        REQUIRE(SequenceAlignment::deviceType == SequenceAlignment::programArgs::CPU);
        REQUIRE(SequenceAlignment::sequenceType == SequenceAlignment::programArgs::PROTEIN);
    }

    SECTION("reading data from files")
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
}
